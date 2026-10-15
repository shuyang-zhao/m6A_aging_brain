# Gene models, RNA modification site maps (hg38, mm10), site atlas
# reads data/raw (00_download.sh); writes output/prepared/resources/gencode_<species> and site_atlas
import gzip
import io
import re
import shutil
import tarfile
import zipfile

import numpy as np
import pandas as pd
import pysam
from pyliftover import LiftOver

import utils

RAW, OUT = utils.DATA / "raw", utils.RESOURCES
KEY = ["chrom", "pos", "strand"]
SITE_BASE = {"m6A": "A", "Psi": "T", "m5C": "C", "m7G": "G"}  # base of the modified nucleotide (Psi from U)
genomes = {}
to_hg38, to_mm10 = (LiftOver(str(RAW / f"ucsc_liftover/{name}.over.chain.gz")) for name in ("hg19ToHg38", "mm39ToMm10"))


def fetch(species, chrom, start, end, strand):
    # interval bases, reverse-complemented on the minus strand
    genome = genomes[species]
    if chrom not in genome.references:
        return ""
    bases = genome.fetch(chrom, max(start, 1) - 1, min(end, genome.get_reference_length(chrom))).upper()
    return bases.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1] if strand == "-" else bases


def site_table(frame, **extra):
    # site map columns chrom, pos, strand, stoich, plus extras
    table = pd.DataFrame(
        {
            "chrom": frame.chrom.astype(str).str.replace(r"^(?:chr)?", "chr", regex=True),
            "pos": frame.pos.astype("int64"),
            "strand": frame.strand,
            "stoich": frame.get("stoich", np.nan),
        }
    )
    extra = {name: values.to_numpy() if hasattr(values, "to_numpy") else values for name, values in extra.items()}
    return table.reset_index(drop=True).assign(**extra)


def best(frame):
    # highest stoichiometry per site
    return frame.groupby(KEY, as_index=False).stoich.max()


def read_xls(path, **kwargs):
    # gzip-packed xls; two header rows joined as upper|lower
    raw = pd.read_excel(io.BytesIO(gzip.open(path).read()), engine="xlrd", header=None, **kwargs)
    upper, lower = raw.iloc[0].ffill(), raw.iloc[1]
    raw.columns = [f"{head}|{name}" if isinstance(head, str) else name for head, name in zip(upper, lower, strict=True)]
    return raw.iloc[2:]


def tx_offset(tx_id, pos, strand, exons):
    # offset along the transcript of a position
    # from the first covering exon on the same strand
    sites = pd.DataFrame({"tx_id": tx_id, "pos": pos, "strand": strand}).dropna(subset=["tx_id"])
    hit = sites.reset_index(names="row").merge(exons, on="tx_id", suffixes=("", "_exon"))
    hit = hit[(hit.start <= hit.pos) & (hit.end >= hit.pos) & (hit.strand == hit.strand_exon)].drop_duplicates("row")
    offset = hit.tx_offset + np.where(hit.strand == "+", hit.pos - hit.start, hit.end - hit.pos)
    return pd.Series(offset.to_numpy(), index=hit.row)


# Canonical transcript per coding gene: exons, spliced sequence
# human v44: MANE, Ensembl canonical, APPRIS, longest CDS
# mouse vM25 with the vM33 canonical ids
gtfs, models = {}, {}
for species, gtf_file, genome_file in (
    ("human", "gencode_v44/gencode.v44.basic.annotation.gtf.gz", "gencode_v44/GRCh38.primary_assembly.genome.fa.gz"),
    ("mouse", "gencode_vM25/gencode.vM25.basic.annotation.gtf.gz", "gencode_vM25/GRCm38.primary_assembly.genome.fa.gz"),
):
    gtf = gtfs[species] = utils.read_gtf(RAW / gtf_file, ("gene", "transcript", "exon", "CDS", "stop_codon", "UTR"))
    # coding genes without the chrY PAR copies
    genes = gtf[(gtf.feat == "gene") & (gtf.gene_type == "protein_coding")]
    genes = genes[~genes.gene_id.str.endswith("_PAR_Y")].set_index("gene_id")
    # their coding transcripts, with APPRIS and lengths for the ranking
    tx = gtf[(gtf.feat == "transcript") & (gtf.tx_type == "protein_coding") & gtf.gene_id.isin(genes.index)].copy()
    exons = gtf[(gtf.feat == "exon") & gtf.tx_id.isin(tx.tx_id)]
    cds = gtf[(gtf.feat == "CDS") & gtf.tx_id.isin(tx.tx_id)]
    tx["appris"] = tx.tags.str.extract(r"appris_principal_(\d)")[0].astype(float).fillna(9)
    tx["cds_len"] = tx.tx_id.map((cds.end - cds.start + 1).groupby(cds.tx_id).sum()).fillna(0)
    tx["tx_len"] = tx.tx_id.map((exons.end - exons.start + 1).groupby(exons.tx_id).sum())
    if species == "human":
        mane, canonical = tx.tags.str.contains("MANE_Select"), tx.tags.str.contains("Ensembl_canonical")
        tx["rank"] = np.select([mane, canonical], [0, 1], 2)
        source = {0: "MANE_Select", 1: "Ensembl_canonical", 2: "APPRIS/longest CDS"}
    else:
        newer = utils.read_gtf(RAW / "gencode_vM33/gencode.vM33.basic.annotation.gtf.gz", ("transcript",))
        canonical = set(newer.tx_id[newer.tags.str.contains("Ensembl_canonical")].str.split(".").str[0])
        tx["rank"] = np.where(tx.tx_id.str.split(".").str[0].isin(canonical), 0, 2)
        source = {0: "Ensembl_canonical (vM33 id)", 2: "APPRIS/longest CDS"}
    # best transcript first: canonical tag, APPRIS, longest CDS and transcript
    order = ["gene_id", "rank", "appris", "cds_len", "tx_len"]
    tx = tx.sort_values(order, ascending=[True, True, True, False, False]).drop_duplicates("gene_id")

    # exons 5' to 3' with transcript offsets
    exons = exons[exons.tx_id.isin(tx.tx_id)]
    exons = exons.assign(key=np.where(exons.strand == "+", exons.start, -exons.end)).sort_values(["tx_id", "key"])
    exons["exon_no"] = exons.groupby("tx_id").cumcount() + 1
    exons["len"] = exons.end - exons.start + 1
    exons["tx_offset"] = exons.groupby("tx_id").len.cumsum() - exons.len
    exons = exons[["tx_id", "exon_no", "chrom", "start", "end", "strand", "len", "tx_offset"]]

    # CDS ends (stop codon included) as transcript offsets
    coding = gtf[gtf.feat.isin(["CDS", "stop_codon"]) & gtf.tx_id.isin(tx.tx_id)]
    span = coding.groupby("tx_id").agg(low=("start", "min"), high=("end", "max"))
    strand = tx.set_index("tx_id").strand.reindex(span.index)
    five = span.low.where(strand == "+", span.high)
    three = span.high.where(strand == "+", span.low)
    cds_start = tx_offset(span.index.to_series(), five, strand, exons)
    cds_end = tx_offset(span.index.to_series(), three, strand, exons) + 1

    # transcript table with CDS limits and UTR lengths
    transcripts = tx[["gene_id", "gene_name", "tx_id", "chrom", "strand", "start", "end"]]
    transcripts = transcripts.rename(columns={"start": "tx_start", "end": "tx_end"})
    transcripts["canonical_source"] = tx["rank"].map(source)
    transcripts["gene_start"] = transcripts.gene_id.map(genes.start)
    transcripts["gene_end"] = transcripts.gene_id.map(genes.end)
    transcripts["cds_start"] = transcripts.tx_id.map(cds_start)
    transcripts["cds_end"] = transcripts.tx_id.map(cds_end)
    transcripts["n_exons"] = transcripts.tx_id.map(exons.groupby("tx_id").size())
    transcripts["tx_len"] = transcripts.tx_id.map(exons.groupby("tx_id").len.sum())
    transcripts = transcripts.dropna(subset=["cds_start", "cds_end"]).astype({"cds_start": int, "cds_end": int})
    transcripts["utr5_len"] = transcripts.cds_start
    transcripts["cds_len"] = transcripts.cds_end - transcripts.cds_start
    transcripts["utr3_len"] = transcripts.tx_len - transcripts.cds_end
    transcripts["gene_id"] = transcripts.gene_id.str.split(".").str[0]
    exons = exons[exons.tx_id.isin(transcripts.tx_id)]
    # exons of the protein-coding transcripts of coding genes, for the gene lookup
    gene_exons = gtf[(gtf.feat == "exon") & (gtf.gene_type == "protein_coding") & (gtf.tx_type == "protein_coding")]
    gene_exons = gene_exons.assign(gene_id=gene_exons.gene_id.str.split(".").str[0])
    gene_exons = gene_exons[["chrom", "start", "end", "strand", "gene_id"]].drop_duplicates()

    # plain FASTA with index, then spliced sequences
    folder = OUT / f"gencode_{species}"
    folder.mkdir(parents=True, exist_ok=True)
    with gzip.open(RAW / genome_file, "rb") as packed, open(folder / "genome.fa", "wb") as plain:
        shutil.copyfileobj(packed, plain, 1 << 24)
    pysam.faidx(str(folder / "genome.fa"))
    genomes[species] = pysam.FastaFile(str(folder / "genome.fa"))
    pieces = dict(list(exons.reset_index(drop=True).groupby("tx_id")))
    spliced = [
        "".join(fetch(species, e.chrom, e.start, e.end, e.strand) for e in pieces[tx_id].itertuples())
        for tx_id in transcripts.tx_id
    ]
    sequences = pd.DataFrame(
        {"gene_id": transcripts.gene_id.to_numpy(), "tx_id": transcripts.tx_id.to_numpy(), "seq": spliced}
    )
    transcripts, exons, gene_exons = (table.reset_index(drop=True) for table in (transcripts, exons, gene_exons))
    assert transcripts.gene_id.is_unique and len(transcripts) > 19000
    models[species] = transcripts, exons, gene_exons
    utils.write(transcripts, folder / "transcripts.parquet")
    utils.write(exons, folder / "exons.parquet")
    utils.write(gene_exons, folder / "gene_exons.parquet")
    utils.write(sequences, folder / "sequences.parquet")

maps = {}  # map name -> table; the atlas list at the end gives each map its label


def merge_ids(first, second):
    # ids of two replicates joined by ";" where they differ
    return [";".join(sorted({a, b} - {""})) for a, b in zip(first, second, strict=True)]


# m6A in HEK293T and HeLa (GLORI, Liu 2023)
# sites in both replicates, RefSeq mRNAs only
files = sorted((RAW / "liu2023_gse210563").glob("*.csv.gz"))
lines = []
for pair in zip(files[0::2], files[1::2], strict=True):
    reps = []
    for file in pair:
        calls = pd.read_csv(file, sep="\t", dtype={"Gene": str, "Transcript": str}, keep_default_na=False)
        calls = calls[(calls.P_adjust <= 0.01) & (calls.AGcov >= 20) & (calls.NormeRatio >= 0.05)]
        order = ["P_adjust", "AGcov", "NormeRatio", "Transcript"]
        calls = calls.sort_values(order, ascending=[True, False, False, True], kind="stable")
        reps.append(calls.drop_duplicates(["Chr", "Sites", "Strand"]).set_index(["Chr", "Sites", "Strand"]))
    shared = reps[0].index.intersection(reps[1].index)
    first, second = reps[0].loc[shared], reps[1].loc[shared]
    merged = {
        "transcript": merge_ids(first.Transcript, second.Transcript),
        "gene": merge_ids(first.Gene, second.Gene),
        "stoich": (first.NormeRatio + second.NormeRatio) / 2,
    }
    lines.append(pd.DataFrame(merged, index=shared).reset_index(names=KEY))
# RefSeq mRNA (NM_) sites, one row per site
sites = pd.concat(lines)
sites = sites[sites.transcript.str.startswith("NM_")].groupby([*KEY, "gene"], as_index=False).stoich.mean()
maps["m6a_hek293t_hela"] = site_table(sites.drop_duplicates(KEY))

# m5C in HeLa and HEK293T (TAC-seq, Lu 2024, hg19)
# plus the NSUN2 and NSUN6 dependent sites
book = RAW / "lu2024_gse242724/1-s2.0-S1097276524005288-mmc5.xlsx"
sheets = ("HeLa_cutoff-1", "HEK293T_cutoff-1")
table = pd.concat([pd.read_excel(book, sheet_name=sheet) for sheet in sheets], ignore_index=True)
sites = utils.lift(to_hg38, table.Chromosome, table.Position, table.Strand)
# mean TAC-seq ratio, percent to fraction
sites["stoich"] = table.filter(like="C-to-T_ratio_m5C-TAC").mean(axis=1) / 100
sites["writer"] = table["Enzyme-dependent"].astype(str)
sites = sites.dropna(subset=["chrom"]).astype({"pos": "int64"})
maps["m5c_tac_union"] = site_table(best(sites))
# writer-dependent subsets
for writer in ("NSUN2", "NSUN6"):
    maps[f"m5c_writer_{writer}"] = site_table(best(sites[sites.writer == writer]), writer_label=writer)

# m5C in HEK293T, HeLa (UBS-seq, Dai 2024), ratio >= 5%
calls = pd.concat(
    [
        pd.read_csv(RAW / f"dai2024_gse225614/GSE225614_{line}-WT_sites.tsv.gz", sep="\t", na_values=["."]).dropna(
            subset=["chromosome", "position"]
        )
        for line in ("HEK293T", "HeLa")
    ]
)
calls = calls[calls.ratio >= 0.05]
contig = calls.chromosome.astype(str).str.replace(r"^chr", "", regex=True).replace("MT", "M")
ubs = pd.DataFrame({"chrom": contig, "pos": calls.position, "strand": calls.strand, "stoich": calls.ratio})
maps["m5c_ubs_union"] = site_table(best(ubs))

# m5C in frontal cortex (RNA-BisSeq, Huang 2019, hg19)
table = read_xls(RAW / "modmaps/GSE122254/GSE122254_human_sites.xls.gz", sheet_name="Human tissue")
table = table[table["Frontal cortex|if m5C site"].astype(str).str.lower() == "true"]
sites = utils.lift(to_hg38, "chr" + table.Chromosome.astype(str), table.Position.astype(int), table.Strand)
sites["stoich"] = pd.to_numeric(table["Frontal cortex|m5C level"], errors="coerce").to_numpy()
sites = sites.dropna(subset=["chrom"]).drop_duplicates(KEY).astype({"pos": "int64"})
maps["m5c_frontal_cortex_huang2019"] = site_table(sites)

# m7G in HeLa and HepG2 (m7G-seq, Zhang 2019, hg38)
columns = {"usecols": [0, 1, 5, 10], "names": ["chrom", "start", "strand", "freq"]}
calls = [
    pd.read_csv(RAW / f"zhang2019_m7g_scimodom/{line}_m7G.bedrmod", sep="\t", comment="#", header=None, **columns)
    for line in ("hela", "hepg2")
]
sites = pd.concat(calls).groupby(["chrom", "start", "strand"], as_index=False).freq.max()
maps["m7g_hela_hepg2"] = site_table(sites.assign(pos=sites.start + 1, stoich=sites.freq / 100))

# METTL1-dependent m7G regions from MeRIP peaks (hg19 lifted)
# overlapping peaks form a region
# region = intersection of its peaks (narrowest if empty)
# kept when a peak is lost on METTL1 knock-down
zhang = ("HeLa_HighFC", "HepG2_HighFC", "HEK293T_HighFC", "HeLa_OverallFC", "HepG2_OverallFC")
zhang += ("HeLa_siMETTL1_ConsistentDiff_peak", "HepG2_shMETTL1_ConsistentDiff_peak")
zhao = ("94_r1", "95_r2", "98_m7g.1", "99_m7g.2")
peak_files = [(f"zhang2019_gse112276/GSE112276_{name}", "hg19", "Diff" in name) for name in zhang]
peak_files += [(f"zhao2023_gse193039/GSM57667{name}.mrna.peak", "hg38", False) for name in zhao]
parts = []
for path, build, differential in peak_files:
    table = pd.read_csv(RAW / f"{path}.txt.gz", sep="\t")
    part = pd.DataFrame({"chrom": table.chr, "start": table.chromStart, "end": table.chromEnd, "strand": table.strand})
    part["lost"] = table["diff.log2.fc"] < 0 if differential else False
    if build == "hg19":
        # both ends lifted to the same chromosome and strand, length within two-fold
        kept = {}
        for row, (chrom, start, end) in enumerate(zip(part.chrom, part.start, part.end, strict=True)):
            left = to_hg38.convert_coordinate(chrom, start)
            right = to_hg38.convert_coordinate(chrom, max(end - 1, 0))
            if left and right and left[0][0] == right[0][0] and left[0][2] == right[0][2]:
                low, high = sorted((left[0][1], right[0][1] + 1))
                if 0.5 <= (high - low) / max(end - start, 1) <= 2:
                    kept[row] = (left[0][0], low, high)
        lifted = pd.DataFrame.from_dict(kept, orient="index", columns=["chrom", "start", "end"])
        part = part.loc[lifted.index].assign(**lifted)
    parts.append(part)
peaks = pd.concat(parts).sort_values(["chrom", "strand", "start", "end"], ignore_index=True)
same = ((peaks.chrom == peaks.chrom.shift()) & (peaks.strand == peaks.strand.shift())).to_numpy()
region, reach, number = np.zeros(len(peaks), int), -1, 0
for row, (follows, start, stop) in enumerate(zip(same, peaks.start, peaks.end, strict=True)):
    number, reach = (number, max(reach, stop)) if follows and start < reach else (number + 1, stop)
    region[row] = number
grouped = peaks.groupby(region)
regions = pd.DataFrame(
    {
        "chrom": grouped.chrom.first(),
        "strand": grouped.strand.first(),
        "start": grouped.start.max(),
        "end": grouped.end.min(),
    }
)
narrowest = peaks.loc[(peaks.end - peaks.start).groupby(region).idxmin()].set_index(regions.index)
empty = regions.start >= regions.end
regions.loc[empty, ["start", "end"]] = narrowest.loc[empty, ["start", "end"]]
regions = regions[regions.index.isin(region[peaks.lost.to_numpy()])]
maps["m7g_mettl1_dependent"] = site_table(
    regions.assign(pos=(regions.start + regions.end) // 2), start=regions.start + 1, end=regions.end, interval=True
)

# Psi (hg38): BID-seq, BACS, PRAISE, PRAISE-seq, nanopore
# PRAISE writer labels, copied to single maps by position
WRITERS = ("TRUB1", "PUS7", "DKC1", "PUS1")
human_genes = gtfs["human"][gtfs["human"].feat == "gene"].groupby("chrom")


def gene_strand(chrom, pos):
    # strand of the single gene covering a position
    rows = human_genes.get_group(chrom) if chrom in human_genes.groups else human_genes.head(0)
    found = set(rows.strand[(rows.start <= pos) & (rows.end >= pos)])
    return found.pop() if len(found) == 1 else None


def motif_centre(chrom, pos, strand, motif, reach=6):
    # lowest position of the motif centre within reach, else pos
    half = len(motif) // 2
    start = pos - reach - half - 1
    bases = fetch("human", chrom, start + 1, pos + reach + half, strand)
    hits = [at for at in range(len(bases) - len(motif) + 1) if bases[at : at + len(motif)] == motif]
    return min([start + at + half + 1 if strand == "+" else start + len(bases) - at - half for at in hits], default=pos)


def first_u(window):
    # first U of a PRAISE window
    # single bases take the strand of their gene
    chrom, _, where = str(window).strip().rpartition("_")
    found = re.fullmatch(r"(\d+)(?:-(\d+))?", where)
    if not found or not re.fullmatch(r"chr(?:[1-9]|1\d|2[0-2]|X|Y|M)", chrom):
        return None
    first, last = int(found[1]), int(found[2] or found[1])
    strand = "+" if last > first else "-" if last < first else gene_strand(chrom, first)
    low, high = min(first, last), max(first, last)
    if strand is None or high - low > 9:  # strand unresolved, or window over 10 nt
        return None
    for pos in range(low, high + 1) if strand == "+" else range(high, low - 1, -1):
        if fetch("human", chrom, pos, pos, strand) == "T":
            return chrom, pos, strand


def praise_sheet(name):
    # one PRAISE supplementary sheet, from its header row on
    raw = pd.read_excel(RAW / "zhang2023_praise/41589_2023_1304_MOESM3_ESM.xlsx", sheet_name=name, header=None)
    top = raw.index[raw.iloc[:, 0].astype(str).str.strip() == "chr_name"][0]
    table = raw.iloc[top + 1 :].set_axis([str(column).strip() for column in raw.iloc[top]], axis=1)
    return table[table.chr_name.notna()]


# BID-seq (Dai 2023): positions moved to the motif centre
book = RAW / "dai2023_bidseq/41587_2022_1505_MOESM3_ESM.xlsx"
bid = pd.concat(
    [pd.read_excel(book, sheet_name=f"Supplementary Table {number}", header=3) for number in (5, 6, 7)],
    ignore_index=True,
)
bid = bid.rename(columns={"chr": "chrom"}).assign(published=lambda frame: frame.pos)
bid["pos"] = [motif_centre(*site) for site in zip(bid.chrom, bid.published, bid.strand, bid.Motif_1, strict=True)]
# BACS (Xu 2024)
bacs = pd.read_excel(RAW / "xu2024_bacs/41592_2024_2439_MOESM3_ESM.xlsx", sheet_name="Supplementary Table 8", header=2)
bacs = bacs.rename(columns={"chr": "chrom"})
# PRAISE sites (Zhang 2023): first U of each window
table = praise_sheet("Supplementary Dataset 2")
praise = pd.DataFrame([site for site in map(first_u, table.chr_site) if site], columns=KEY).drop_duplicates(KEY)
# PRAISE writer labels; a site with two writers is multi
labels = {}
table = praise_sheet("Supplementary Dataset 3")
for window, writer in zip(table.chr_site, table.enzyme_dependency.astype(str).str.strip(), strict=True):
    site = first_u(window)
    if site and writer in WRITERS:
        labels[site] = "multi" if labels.get(site, writer) != writer else writer
# nanopore calls (McCormick 2026): strand from the reference k-mer
with zipfile.ZipFile(RAW / "mccormick2026_nar/gkag353_supplemental_files.zip") as bundle:
    table = pd.read_excel(io.BytesIO(bundle.read("SupplementaryTableS6.xlsx")), sheet_name="HepG2")
strands = []
for chrom, pos, kmer in zip(table.chr, table.position, table.kmer, strict=True):
    kmer, half = kmer.upper().replace("U", "T"), len(kmer) // 2
    forward = fetch("human", chrom, pos - half, pos + half, "+") == kmer
    reverse = fetch("human", chrom, pos - half, pos + half, "-") == kmer
    strands.append("+" if forward and not reverse else "-" if reverse and not forward else None)
nanopore = pd.DataFrame({"chrom": table.chr, "pos": table.position, "strand": strands})
nanopore = nanopore.dropna(subset=["strand"]).drop_duplicates(KEY)
# union of the maps with labels and moved positions
union = pd.concat([frame[KEY] for frame in (nanopore, bacs, bid.assign(pos=bid.published), praise)]).drop_duplicates()
keys = list(union[KEY].itertuples(index=False, name=None))
moved = {
    (chrom, published, strand): pos
    for chrom, published, strand, pos in zip(bid.chrom, bid.published, bid.strand, bid.pos, strict=True)
}
union["writer_label"] = [labels.get(key, "") for key in keys]
union["pos"] = [moved.get(key, key[1]) for key in keys]
union = union.sort_values("writer_label", ascending=False, kind="stable").drop_duplicates(KEY)
labelled = union[union.writer_label != ""]
maps["psi_union"] = site_table(union, writer_label=union.writer_label.to_numpy())
for writer in WRITERS:
    maps[f"psi_writer_{writer}"] = site_table(labelled[labelled.writer_label == writer], writer_label=writer)
# labelled sites by position
known = {
    (chrom, pos, strand): writer
    for chrom, pos, strand, writer in labelled[[*KEY, "writer_label"]].itertuples(index=False, name=None)
}
# PRAISE-seq (HEK293T): coding genes, no intron or intergenic sites
praiseseq = pd.read_excel(
    RAW / "praiseseq2026_gse326442/GSE326442_PRAISE-seq-HEK293T.xlsx", sheet_name="pu_premRNA_anno"
)
praiseseq = praiseseq[(praiseseq.Gene_Type == "protein_coding") & ~praiseseq.annotation.isin(["Intron", "Intergenic"])]
best_bid = bid.groupby(KEY, as_index=False)["Frac_Ave %"].max()
singles = {
    "psi_hek293t": praiseseq.rename(columns={"chr": "chrom", "ratio": "stoich"}),
    "psi_bidseq_3cells": best_bid.assign(stoich=best_bid["Frac_Ave %"] / 100),
    "psi_bacs_hela": bacs.assign(stoich=bacs["Ψ level"]),
}
# singles take the label of a labelled site
# up to 2 bases away along a run of U
for name, single in singles.items():
    single, label = site_table(single), []
    for chrom, pos, strand in zip(single.chrom, single.pos, single.strand, strict=True):
        label.append(known.get((chrom, pos, strand), ""))
        for shift in (1, -1, 2, -2):
            near = known.get((chrom, pos + shift, strand))
            if not label[-1] and near:
                run = fetch("human", chrom, min(pos, pos + shift), max(pos, pos + shift), strand)
                label[-1] = near if set(run) == {"T"} else ""
    maps[name] = single.assign(writer_label=label)

# Mouse maps (mm10)
# GLORI and BID-seq sites kept if A or U
# in coding CDS or UTR, off exon starts
# m5C and m7G keep their published sites
gtf = gtfs["mouse"]
gtf["first"] = np.where(gtf.strand == "+", gtf.start, gtf.end)
coding = gtf[gtf.gene_type == "protein_coding"]
coding_genes = coding[coding.feat == "gene"]
symbols = set(coding_genes.gene_id.str.split(".").str[0]) | set(coding_genes.gene_name)


def blocks(rows, start, end):
    # sorted block starts with running maximum of ends
    out = {}
    for key, group in rows.groupby(["chrom", "strand"]):
        order = np.argsort(group[start].to_numpy())
        out[key] = (group[start].to_numpy()[order], np.maximum.accumulate(group[end].to_numpy()[order]))
    return out


def inside(sites, block):
    # inside: the running maximum end reaches the site
    position, mask = sites.pos.to_numpy(), np.zeros(len(sites), bool)
    for key, rows in sites.groupby(["chrom", "strand"]).indices.items():
        if key in block:
            at = np.searchsorted(block[key][0], position[rows], side="right") - 1
            mask[rows] = (at >= 0) & (position[rows] <= block[key][1][at])
    return mask


# exons of coding genes, CDS and UTR blocks, exon starts
exon_blocks = blocks(gtf[(gtf.feat == "exon") & (gtf.gene_type == "protein_coding")], "start", "end")
coding_blocks = blocks(coding[coding.feat.isin(["CDS", "UTR"])], "start", "end")
start_blocks = blocks(gtf[gtf.feat == "exon"], "first", "first")


def mouse_sites(table, chrom, pos, strand, stoich=None, per=1.0, keep=()):
    # one row per site, highest stoichiometry kept
    frame = pd.DataFrame({"chrom": table[chrom], "pos": table[pos], "strand": table[strand]})
    if stoich:
        frame["stoich"] = table[stoich].astype(float) / per
    sites = site_table(frame, **{name: table[name] for name in keep})
    sites = sites.sort_values("stoich", ascending=False, kind="stable")
    return sites.drop_duplicates(KEY).sort_index().reset_index(drop=True)


def coding_sites(sites, base):
    # standard chromosomes, CDS or UTR, off exon starts, base matches
    sites = sites[sites.chrom.str.fullmatch(r"chr(\d+|X|Y)")]
    sites = sites[inside(sites, coding_blocks) & ~inside(sites, start_blocks)].reset_index(drop=True)
    read = [
        fetch("mouse", chrom, pos, pos, strand)
        for chrom, pos, strand in zip(sites.chrom, sites.pos, sites.strand, strict=True)
    ]
    return sites[pd.Series(read, index=sites.index) == base].reset_index(drop=True)


def consensus(reps):
    # sites seen in at least two replicates, mean stoichiometry
    calls = pd.concat(reps).groupby(KEY, as_index=False).stoich.agg(["size", "mean"]).rename(columns={"mean": "stoich"})
    return calls[calls["size"] >= 2].drop(columns="size")


# m6A in mouse hippocampus (GLORI 3.0, Sun 2025)
# cytoplasm and synaptosome, sites in two replicates
with tarfile.open(RAW / "modmaps/GSE270643/GSE270643_RAW.tar") as tar:
    for name, label in (("m6a_hippocampus_cyto", "GLORI_3_0_Cyto_rep"), ("m6a_hippocampus_syn", "GLORI_3_0_Syn_rep")):
        reps = []
        for member in sorted(tar.getnames()):
            if label in member:
                columns = ["Chr", "Sites", "Strand", "Ratio"]
                calls = pd.read_csv(tar.extractfile(member), sep="\t", compression="gzip", usecols=columns)
                reps.append(mouse_sites(calls, "Chr", "Sites", "Strand", "Ratio"))
        maps[name] = coding_sites(consensus(reps), "A")

# Psi in mouse cortex and cerebellum (BID-seq, Dai 2023)
# stem cell sites (Zhang 2024): GRCm39, GEO says mm10
# (U only at chance in mm10), so they are lifted
book = RAW / "modmaps/GSE179798/GSE179798_Mouse_mRNA_12-Tissues_BID-seq.xlsx"
for name, sheet in (("psi_cortex", "Cerebral cortex"), ("psi_cerebellum", "Cerebellum")):
    table = pd.read_excel(book, sheet_name=sheet, header=3)
    maps[name] = coding_sites(mouse_sites(table[table.chr.notna()], "chr", "pos", "strand", "Frac_Ave %", per=100), "T")
table = pd.read_excel(RAW / "modmaps/GSE238245/GSE238245_Pseudouridine_mESC-WT_polyA-RNA.xlsx", header=1)
table = table[table.chr.notna() & (table["Ref Type"] == "genome")]
sites = mouse_sites(table, "chr", "position", "strand", "mESCWT_fraction")
lifted = utils.lift(to_mm10, sites.chrom, sites.pos, sites.strand)
mapped = lifted.chrom.notna().to_numpy()
sites = sites[mapped].assign(pos=lifted.pos[mapped].astype(int), strand=lifted.strand[mapped])
maps["psi_mesc"] = coding_sites(sites, "T")
assert len(maps["psi_mesc"]) > 1000  # most sites are U only after the lift to mm10

# Psi in neural stem cells (Uli-epic, He 2025), both replicates
reps = []
for rep in (1, 2):
    table = pd.read_csv(RAW / f"modmaps/GSE277358/GSM85209{15 + rep}_WT-NSC-Uli-epic-BIDSeq-rep{rep}.xls.gz", sep="\t")
    passed = table.filter(regex="_treat_passed$").iloc[:, 0]
    table = table.assign(fraction=table.filter(regex="_treat_fraction$").iloc[:, 0])
    reps.append(mouse_sites(table[passed == 1], "chr", "pos", "strand", "fraction"))
maps["psi_nsc_uliepic"] = coding_sites(consensus(reps), "T")

# m5C in mouse brain (RNA-BisSeq): Yang 2017, Amort 2017
# Yang: mRNA sites; Amort: sites inside exons
table = read_xls(
    RAW / "modmaps/GSE93749/GSE93749_Mouse_m5C_sites_information.xls.gz", sheet_name="Mouse_6_tissues_total_m5C_sites"
)
level = "Brain_8W-BS|Final Level"
genic = table["Gene ID"].isin(symbols) | table["Gene Symbol"].isin(symbols)
measured = pd.to_numeric(table[level], errors="coerce").notna()
maps["m5c_brain_yang2017"] = mouse_sites(
    table[genic & measured & (table["Gene Type"] == "mRNA")], "Chromosome", "Position", "Strand", level
)
table = pd.read_excel(RAW / "modmaps/GSE83432/GSE83432_Brain_total.xlsx", sheet_name="Brain total poly(A)RNA", header=3)
sites = mouse_sites(
    table.set_axis(table.columns.str.lstrip("#"), axis=1), "SeqID", "refPos", "refStrand", "avg/methRate"
)
maps["m5c_brain_amort2017"] = sites[inside(sites, exon_blocks)].reset_index(drop=True)

# m6A in mouse brain (DART-seq, Tegowski 2024)
# the site is the 5' neighbour of the edited C
for name, path, columns in (
    (
        "m6a_cortex_scdart_tegowski2024",
        "GSE266876_youngmousecortex_singlecellsites",
        ["chrom", "s", "e", "n_cells", "gene", "id", "strand"],
    ),
    ("m6a_neuron_dart_tegowski2024", "GSE267234_NeuronNucleisites", ["chrom", "s", "e", "gene", "rate", "strand"]),
):
    bed = pd.read_csv(RAW / f"modmaps/tegowski2024/{path}.bed.gz", sep="\t", header=None, names=columns)
    bed["site"] = bed.s + np.where(bed.strand == "+", 0, 2)
    maps[name] = mouse_sites(bed, "chrom", "site", "strand")

# m7G regions in mouse brain and stem cells (miCLIP, m7GHub)
columns = ["seqnames", "start", "end", "strand", "source", "species", "assembly"]
table = pd.read_csv(RAW / "modmaps/m7ghub2/m7G-miCLIP-seq.csv", usecols=columns)
table = table[(table.species == "mouse") & (table.assembly == "mm10")]
for name, label in (("m7g_brain_miclip", "Mouse Brain"), ("m7g_mesc_miclip", "mESC")):
    rows = table[table.source.str.endswith(label)]
    rows = rows.assign(mid=(rows.start + rows.end) // 2)
    maps[name] = mouse_sites(rows, "seqnames", "mid", "strand", keep=("start", "end")).assign(interval=True)

# m6A in five human brain regions (m6A-seq peaks of 25 donors, Shafik 2026; hg38, 1-based starts)
peaks = pd.read_excel(RAW / "modmaps/shafik2026/41593_2025_2112_MOESM2_ESM.xlsx", sheet_name="Table 2", usecols="A:C,F")
peaks.columns = ["chrom", "start", "end", "strand"]
peaks = peaks.assign(pos=(peaks.start + peaks.end) // 2)
maps["m6a_brain_shafik2026"] = site_table(peaks, start=peaks.start, end=peaks.end, interval=True)

# Site atlas: maps placed on the canonical transcripts
# reference base, gene, offset, region, codon, motif, writer class
atlas = {
    "human": {
        "m6a_hek293t_hela": "HEK293T + HeLa (GLORI consensus)",
        "m6a_brain_shafik2026": "human brain, five regions (m6A-seq, Shafik 2026)",
        "m5c_tac_union": "HeLa + HEK293T (TAC-seq)",
        "m5c_writer_NSUN2": "HeLa + HEK293T (TAC-seq, NSUN2-dependent)",
        "m5c_writer_NSUN6": "HeLa + HEK293T (TAC-seq, NSUN6-dependent)",
        "m5c_ubs_union": "human cell lines (UBS-seq)",
        "m5c_frontal_cortex_huang2019": "human frontal cortex (RNA-BisSeq, Huang 2019)",
        "m7g_hela_hepg2": "HeLa + HepG2 (base-resolution m7G)",
        "m7g_mettl1_dependent": "human cell lines (METTL1-dependent regions)",
        "psi_hek293t": "HEK293T (PRAISE-seq)",
        "psi_union": "human cell lines (union)",
        **{f"psi_writer_{writer}": f"human cell lines (union, {writer}-dependent subset)" for writer in WRITERS},
        "psi_bidseq_3cells": "HeLa, HEK293T, A549 (BID-seq)",
        "psi_bacs_hela": "HeLa (BACS)",
    },
    "mouse": {
        "m6a_hippocampus_cyto": "mouse hippocampus, cytoplasmic fraction (GLORI 3.0)",
        "m6a_hippocampus_syn": "mouse hippocampus, synaptosome fraction (GLORI 3.0)",
        "psi_cortex": "mouse cerebral cortex (BID-seq, 12-tissue)",
        "psi_cerebellum": "mouse cerebellum (BID-seq, 12-tissue)",
        "psi_mesc": "mouse ESC (polyA Psi)",
        "psi_nsc_uliepic": "mouse NSC (Uli-Epic BID-seq)",
        "m5c_brain_yang2017": "mouse brain (RNA-BisSeq, Yang 2017)",
        "m5c_brain_amort2017": "mouse brain (RNA-BisSeq, Amort 2017)",
        "m6a_cortex_scdart_tegowski2024": "mouse cortex, young (scDART-seq single-cell sites, GSE266876)",
        "m6a_neuron_dart_tegowski2024": "mouse neuronal nuclei (DART-seq, GSE267234)",
        "m7g_brain_miclip": "mouse brain (m7G-miCLIP, m7GHub compilation)",
        "m7g_mesc_miclip": "mouse ESC (m7G-miCLIP)",
    },
}
MODIFICATION = {name.lower(): name for name in utils.MODIFICATIONS}
# motif classes of the 15-nt flank (site at base 7)
# as (class, offset, pattern)
motifs = {
    "Psi": [("PUS7_motif", 5, "T..A[AG]"), ("TRUB1_motif", 5, "GT.C.A..C"), ("TRUB1_motif", 6, "G.TC.A..C")],
    "m5C": [("NSUN2_motif", 9, "GGG"), ("NSUN6_motif", 4, "CTC.A")],
    "m6A": [("DRACH", 5, "[AGT][AG].C[ACT]")],
}
# fmt: off
COLUMNS = [
    "map", "modification", "species", "chrom", "start", "end", "strand", "resolution", "stoich", "gene_id", "tx_id",
    "region", "source", "pos", "ref_ok", "tx_offset", "motif_class",
]
# fmt: on


def motif_class(flank, modification):
    if len(flank) != 15 or modification not in motifs:
        return ""
    return next((name for name, start, pattern in motifs[modification] if re.match(pattern, flank[start:])), "other")


def covering_gene(sites, gene_exons):
    # gene whose exons cover a site on its strand
    # several genes make the site ambiguous
    found = np.full(len(sites), None, dtype=object)
    for (chrom, strand), rows in sites.groupby(["chrom", "strand"]):
        exons = gene_exons[(gene_exons.chrom == chrom) & (gene_exons.strand == strand)].sort_values("start")
        starts, ends, genes = exons.start.to_numpy(), exons.end.to_numpy(), exons.gene_id.to_numpy()
        reach = np.maximum.accumulate(ends)
        for row, pos, at in zip(rows.index, rows.pos, np.searchsorted(starts, rows.pos, side="right") - 1, strict=True):
            hits = set()
            while at >= 0 and reach[at] >= pos:
                if ends[at] >= pos:
                    hits.add(genes[at])
                at -= 1
            if hits:
                found[row] = hits.pop() if len(hits) == 1 else "ambiguous"
    return pd.Series(found, index=sites.index)


assert set(maps) == {name for entries in atlas.values() for name in entries}
tables = []
for species, entries in atlas.items():
    transcripts, exons, gene_exons = models[species]
    by_gene = transcripts.set_index("gene_id")
    for name, tissue in entries.items():
        modification = MODIFICATION[name[:3]]
        resolution = "interval" if "interval" in maps[name] else "single_base"
        sites = maps[name].drop_duplicates(KEY).reset_index(drop=True)
        if resolution == "single_base":
            sites["flank"] = [
                fetch(species, chrom, pos - 7, pos + 7, strand)
                for chrom, pos, strand in zip(sites.chrom, sites.pos, sites.strand, strict=True)
            ]
            sites["ref_ok"] = sites.flank.str[7] == SITE_BASE[modification]
        else:
            sites["flank"], sites["ref_ok"] = "", True
        gene = covering_gene(sites, gene_exons)
        sites["gene_id"] = gene.where(gene != "ambiguous")
        sites["tx_id"] = sites.gene_id.map(by_gene.tx_id)
        offset = tx_offset(sites.tx_id, sites.pos, sites.strand, exons).reindex(sites.index)
        cds_start, cds_end = sites.gene_id.map(by_gene.cds_start), sites.gene_id.map(by_gene.cds_end)
        stop = (offset >= cds_end - 100) & (offset < cds_end + 100)  # 100 nt either side of the CDS end
        region = np.select([offset < cds_start, stop, offset < cds_end], ["5UTR", "stop", "CDS"], "3UTR")
        sites["region"] = np.where(offset.isna(), np.where(sites.gene_id.notna(), "noncanonical", None), region)
        sites["tx_offset"] = offset
        sites["motif_class"] = [motif_class(flank, modification) for flank in sites.flank]
        for column in ("start", "end"):
            sites[column] = sites.get(column, sites.pos).fillna(sites.pos).astype("int64")
        sites["map"] = name
        sites["modification"] = modification
        sites["resolution"] = resolution
        sites["species"] = {"human": "Homo sapiens", "mouse": "Mus musculus"}[species]
        sites["source"] = tissue
        tables.append(sites[COLUMNS])
sites = pd.concat(tables, ignore_index=True)
# reference base as expected except in DART-seq maps
# a gene found for most sites
assert sites[~sites["map"].str.contains("dart")].groupby("map").ref_ok.mean().gt(0.99).all()
assert sites.gene_id.notna().groupby(sites["map"]).mean().gt(0.6).all()
utils.write(sites, OUT / "site_atlas/sites.parquet")
