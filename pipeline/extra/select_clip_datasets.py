# One CLIP dataset per reader and species: human ENCODE eCLIP first, else the GEO dataset of the clipflow manifest
# (clipflow: the CLIP library catalogue of a sibling project, not part of this repository). The m6A and m5C readers
# were selected and processed together; this study uses the 18 m6A reader sets (script 07)
# reads the clipflow library index and sample sheets (CLIPFLOW below), the ENCODE portal, ENA and NCBI E-utilities
# writes data/raw/clip_standard/selection.tsv and no_set.tsv
#
# GEO order: brain tissue, then an input control of the same condition, then the most IP libraries, then the most
# crosslink events; two IP replicates. Only protocols that nf-core/clipseq processes (single-end reads whose 5' end is
# the RT stop: iCLIP, eCLIP, seCLIP, irCLIP) are eligible; a reader without one has no set. The crosslink read and the
# UMI layout of each dataset are the ones measured by clipflow.
import io
import json
import re
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import utils  # noqa: E402

CLIPFLOW = Path.home() / "projects/clipflow/data"  # checkout of the sibling project
OUT = utils.DATA / "raw/clip_standard"
CORPUS = {"human": "clip_corpus", "mouse": "clip_corpus_mouse"}
SHEETS = [*(CLIPFLOW / "batches").glob("*.csv"), CLIPFLOW / "geo_input_rerun/geo_full_input_fixed.csv"]
QC_FAILED = {"GSE107767"}  # eCLIP mirror of ENCODE with undetermined UMIs
BRAIN = ("brain", "hippocamp", "cortex", "cerebell", "striat")
SUPPORTED = {"eCLIP", "seCLIP", "iCLIP", "iCLIP2", "irCLIP"}
# protocol of datasets labelled CLIP-seq or not at all in the manifest, from GEO and the reads
METHOD_OF = {"GSE115429": "irCLIP", "PRJNA961171": "eCLIP"}
# libraries of another genotype (reason in the notes)
EXCLUDE = {"GSM3242764": "Mir140 mutant chondrocytes", "GSM3242768": "Mir140 mutant chondrocytes"}
# crosslink read (1 or 2) and UMI of the IP libraries: an umi_tools extract pattern of the read start (X discarded,
# N UMI), the UMI in the read name (prefix "UMI:", suffix "_UMI", or "UMI@" after the read id), or none
READS = {
    "GSE271167": (1, "XXXXXNNNNNNNNXXXXXXXXX"),  # irCLIPv2: barcode, 8 nt UMI, G tail and one extra base
    "GSE115429": (1, "name_suffix"),
    "GSE304795": (1, "NNNNNNNNNN"),
    "GSE315406": (1, "NNNNNNNNNN"),
    "GSE224672": (1, "NNNNNNNNNN"),
    "GSE135642": (1, "name_at"),
    "PRJNA961171": (1, "NNNNNNNNNN"),
    "PRJNA719255": (1, "NNNNNNNNNN"),
    "GSE156115": (1, "NNNNNNNNNN"),
    "GSE98309": (1, "NNNNNNNNNN"),
    "GSE148039": (1, "NNNNNNNNNN"),
    "GSE188846": (1, "none"),  # barcodes removed before deposit
    "GSE196427": (1, "none"),
}
NOTES = {
    ("human", "YTHDF2"): (
        "neural progenitor and glioblastoma stem cell iCLIP (GSE158742, no input) not counted as brain tissue"
    ),
    ("mouse", "ALYREF"): "IP library of the Zc3h11a project (sample CLIP_Alyref, mESC); no input",
    ("mouse", "YBX1"): "wild-type chondrocytes only",
}
# cell type of the GEO sets whose catalogue entry is unknown or a raw tag, from their GEO and BioSample sample records
CELL = {"GSE224672": "H1 hESC", "GSE135642": "HTR8/SVneo", "PRJNA719255": "testis", "GSE156115": "CD11b+ cells",
        "GSE196427": "testis", "GSE188846": "CD4+ T cells"}
GENERIC = re.compile(
    r"^(ip|input|rep|repl|replicate|sminput|sm|size|matched|for|of|clip|iclip|eclip|seclip|irclip)\d*$|^\d+$|^rr\d+$"
)


def fetch(url):
    return urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json"}), timeout=120)


def tag(title):
    # title without replicate, library-type and size-matched words: libraries of one condition share it
    words = re.split(r"[^a-z0-9+]+", str(title).lower())
    return " ".join(w for w in words if w and not GENERIC.match(w))


reader_table = utils.read(utils.ROOT / "ref/gene_sets/rna_mod_regulators.tsv")
reader_table = reader_table[(reader_table.role == "reader") & reader_table.modification.isin(["m6A", "m5C"])]
reader_of = {
    s.upper(): h for h, m in zip(reader_table.human_symbol, reader_table.mouse_symbol, strict=True) for s in (h, m)
}
# readers with an IP library in the manifest (corpus or ENCODE) of human or mouse
pairs = []
for species, corpus in CORPUS.items():
    index = pd.read_csv(CLIPFLOW / corpus / "library_index.tsv", sep="\t", low_memory=False)
    index = index[(index.role == "ip") & (index.status == "ok") & ~index.dataset.isin(QC_FAILED)]
    pairs += [(species, reader) for reader in sorted(index.target.str.upper().map(reader_of).dropna().unique())]
print(f"{len(pairs)} reader-species pairs")


# human ENCODE eCLIP: K562, else HepG2; paired-end reads in both replicates, a size-matched input, earliest released
def encode_experiment(reader):
    url = (
        "https://www.encodeproject.org/search/?type=Experiment&assay_title=eCLIP"
        f"&target.label={reader}&status=released&field=accession&limit=all&format=json"
    )
    try:
        hits = json.load(fetch(url))["@graph"]
    except urllib.error.HTTPError:  # the portal answers 404 to a search without hits
        hits = []
    found = []
    for hit in hits:
        exp = json.load(fetch(f"https://www.encodeproject.org/experiments/{hit['accession']}/?format=json"))
        fastq = [f for f in exp["files"] if f["file_format"] == "fastq"]
        paired = len(fastq) == 4 and all(f["run_type"] == "paired-ended" for f in fastq)
        cell = exp["biosample_ontology"]["term_name"]
        if cell in ("K562", "HepG2") and paired and exp["possible_controls"]:
            found.append((cell != "K562", exp["date_released"], exp))
    return min(found, key=lambda x: x[:2])[2] if found else None


def encode_row(species, reader, exp):
    control = json.load(
        fetch(f"https://www.encodeproject.org/experiments/{exp['possible_controls'][0]['accession']}/?format=json")
    )
    library_of = {r["biological_replicate_number"]: r["library"]["accession"] for r in exp["replicates"]}
    control_library = {r["library"]["accession"] for r in control["replicates"]}
    # read 2 starts at the crosslink and is the single-end input of the pipeline
    files = [f for f in exp["files"] if f["file_format"] == "fastq" and f["paired_end"] == "2"]
    fastq = [f"{library_of[f['biological_replicates'][0]]}|https://www.encodeproject.org{f['href']}" for f in files]
    return dict(
        reader=reader,
        species=species,
        source="ENCODE",
        accession=f"{exp['accession']} (input {control['accession']})",
        cell=exp["biosample_ontology"]["term_name"],
        method="eCLIP",
        ip_libraries=";".join(sorted(library_of.values())),
        control_libraries=";".join(sorted(control_library)),
        crosslink_read=2,
        umi="name_prefix",
        fastq=";".join(sorted(fastq)),
        notes="",
    )


# GEO datasets of the clipflow manifest
sheets = (
    pd.concat([pd.read_csv(p, dtype=str) for p in SHEETS]).drop_duplicates("sample", keep="last").set_index("sample")
)


def ena_runs(library):
    # run accessions and fastq urls of a library: the run accession of the sheet, else the GEO sample alias
    accession = sheets.accession.get(library, library)
    query = (
        f'run_accession="{accession}"' if re.match(r"[SDE]RR\d+$", str(accession)) else f'experiment_alias="{library}"'
    )
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&format=tsv"
        "&fields=run_accession,library_layout,fastq_ftp&query="
    )
    runs = pd.read_csv(io.StringIO(urllib.request.urlopen(url + urllib.parse.quote(query)).read().decode()), sep="\t")
    if runs.empty:  # recent GEO samples are not aliased in ENA: resolve the runs at NCBI
        ids = json.load(
            urllib.request.urlopen(
                f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={library}&retmode=json"
            )
        )
        info = urllib.request.urlopen(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra"
            f"&id={','.join(ids['esearchresult']['idlist'])}&rettype=runinfo&retmode=text"
        )
        runs = pd.concat([ena_runs(run) for run in pd.read_csv(info).Run])
    assert len(runs) > 0, f"no run for {library}"
    return runs


def geo_choice(species, reader):
    index = pd.read_csv(CLIPFLOW / CORPUS[species] / "library_index.tsv", sep="\t", low_memory=False)
    index = index[(index.origin == "clipflow") & ~index.dataset.isin(QC_FAILED) & (index.condition == "untreated")]
    index = index[(index.target.str.upper().map(reader_of) == reader) & index.role.isin(["ip", "input"])]
    index = index[~index.library_id.isin(EXCLUDE)].assign(tag=lambda d: d.title.map(tag))
    index["method"] = index.dataset.map(METHOD_OF).fillna(index.inferred_method)
    groups = []
    for (dataset, cell, group_tag), ips in index[index.role == "ip"].groupby(["dataset", "cell", "tag"]):
        if not ips.method.isin(SUPPORTED).all():
            continue
        same = (index.dataset == dataset) & (index.cell == cell) & (index.tag == group_tag)
        inputs = index[same & (index.role == "input")]
        brain = any(word in str(cell).lower() for word in BRAIN)
        groups.append((-brain, -(len(inputs) > 0), -len(ips), -ips.crosslinks.sum(), dataset, cell, ips, inputs))
    if not groups:
        return None
    _, _, _, _, dataset, cell, ips, inputs = min(groups, key=lambda g: g[:4])
    # two replicates: the first library of each replicate named in the title, else the first two libraries
    ips = ips.sort_values("library_id")
    replicate = ips.title.str.extract(r"(?i)rep(?:licate|l)?\W*(\d+)")[0]
    ips = ips[~replicate.duplicated() | replicate.isna()].head(2)
    return dataset, cell, ips, inputs, len(groups)


def geo_row(species, reader, choice):
    dataset, cell, ips, inputs, n_groups = choice
    mate, umi = READS[dataset]
    fastq = []
    for library in ips.library_id:
        for _, run in ena_runs(library).sort_values("run_accession").iterrows():
            urls = run.fastq_ftp.split(";")
            # paired runs: the file of the mate
            url = next((u for u in urls if u.endswith(f"_{mate}.fastq.gz")), urls[0])
            fastq.append(f"{library}|https://{url}")
    print(species, reader, dataset, cell, f"{len(ips)} IP, {len(inputs)} input, {n_groups} eligible groups")
    return dict(
        reader=reader,
        species=species,
        source="GEO",
        accession=dataset,
        cell=CELL.get(dataset, cell),
        method=ips.method.iloc[0],
        ip_libraries=";".join(ips.library_id),
        control_libraries=";".join(inputs.library_id),
        crosslink_read=mate,
        umi=umi,
        fastq=";".join(fastq),
        notes=NOTES.get((species, reader), ""),
    )


rows, no_set = [], []
for species, reader in pairs:
    experiment = encode_experiment(reader) if species == "human" else None
    if experiment:
        rows.append(encode_row(species, reader, experiment))
    elif choice := geo_choice(species, reader):
        rows.append(geo_row(species, reader, choice))
    else:
        no_set.append(dict(reader=reader, species=species, reason="no single-end truncation protocol in the manifest"))
selection = pd.DataFrame(rows).sort_values(["species", "reader"])
OUT.mkdir(parents=True, exist_ok=True)
utils.write(selection, OUT / "selection.tsv")
utils.write(pd.DataFrame(no_set), OUT / "no_set.tsv")
print(selection.drop(columns="fastq").to_string())
print("no set:", *[f"{r['species']} {r['reader']}" for r in no_set])
