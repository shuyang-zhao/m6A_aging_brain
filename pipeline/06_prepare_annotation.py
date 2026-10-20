# ID maps, GO programmes, lifetimes, half-lives, cell-type origin, GPU predictions, weights of published clocks
# reads data/raw (00_download.sh) and the pseudobulk of 03; writes output/prepared/ref and output/prepared/resources/*
import gzip
import io
import re
import zipfile
from itertools import pairwise

import numpy as np
import pandas as pd

import utils

RAW, OUT = utils.DATA / "raw", utils.RESOURCES

# gene id maps: Ensembl and Swiss-Prot to symbol, orthologs
hgnc = pd.read_csv(RAW / "hgnc/hgnc_complete_set.txt", sep="\t", dtype=str, keep_default_na=False)
mgi_ensembl = pd.read_csv(RAW / "mgi/MRK_ENSEMBL.rpt", sep="\t", header=None, dtype=str, keep_default_na=False)
approved = hgnc[(hgnc.status == "Approved") & (hgnc.ensembl_gene_id != "")]
gene_tables = {"human": approved[["ensembl_gene_id", "symbol"]], "mouse": mgi_ensembl[mgi_ensembl[5] != ""][[5, 1]]}
for species, table in gene_tables.items():
    table = table.set_axis(["gene_id", "symbol"], axis=1).drop_duplicates("gene_id").reset_index(drop=True)
    utils.write(table, utils.PREPARED / f"ref/gene_symbol_{species}.tsv")
    uniprot = pd.read_csv(RAW / f"uniprot/{species}_swissprot_keywords.tsv", sep="\t", dtype=str)
    uniprot = uniprot.rename(columns={"Entry": "uniprot", "Gene Names (primary)": "symbol"})[["uniprot", "symbol"]]
    utils.write(
        uniprot.dropna().drop_duplicates("uniprot").reset_index(drop=True),
        utils.PREPARED / f"ref/uniprot_gene_{species}.tsv",
    )
homology = pd.read_csv(RAW / "mgi/HOM_MouseHumanSequence.rpt", sep="\t", usecols=[0, 1, 3]).set_axis(
    ["key", "organism", "symbol"], axis=1
)
pairs = homology[homology.organism.str.startswith("mouse")].merge(
    homology[homology.organism == "human"], on="key", suffixes=("_mouse", "_human")
)
pairs = (
    pairs[["symbol_mouse", "symbol_human"]]
    .set_axis(["mouse", "human"], axis=1)
    .apply(lambda column: column.str.upper())
    .drop_duplicates()
)
pairs = pairs[~pairs.mouse.duplicated(keep=False) & ~pairs.human.duplicated(keep=False)]
assert len(pairs) > 15000
utils.write(pairs.reset_index(drop=True), utils.PREPARED / "ref/orthologs.tsv")

# Symbol maps for old names: tables of (step, query, value)
# steps: ensembl, official, official_upper (human), synonym, synonym_upper
# to_ensembl takes a name to its gene id (human: also unambiguous previous symbols and aliases; mouse: the
# current MGI symbol)
mgi = pd.read_csv(RAW / "mgi/MRK_List2.rpt", sep="\t", dtype=str, keep_default_na=False)
mgi = mgi[(mgi.Status == "O") & (mgi["Marker Type"] == "Gene")]
mouse = mgi[["MGI Accession ID", "Marker Symbol"]].set_axis(["mgi", "symbol"], axis=1)
mouse = mouse.merge(mgi_ensembl[[0, 5]].set_axis(["mgi", "ensembl"], axis=1), on="mgi", how="left")
mouse_old = mgi[["Marker Symbol", "Marker Synonyms (pipe-separated)"]].set_axis(["symbol", "old"], axis=1)
mouse_names = mgi_ensembl[[1, 5]].set_axis(["name", "ensembl"], axis=1).replace("", np.nan)
human = hgnc[hgnc.status == "Approved"].rename(columns={"ensembl_gene_id": "ensembl"})
human_old = [
    human[["symbol", column]].set_axis(["symbol", "old"], axis=1) for column in ("prev_symbol", "alias_symbol")
]
# human names: unambiguous alias and previous symbols, then official
with_gene = hgnc.replace("", np.nan).dropna(subset=["ensembl_gene_id"])
human_names = []
for column in ("alias_symbol", "prev_symbol"):
    part = with_gene[["ensembl_gene_id", column]].dropna()
    part = part.assign(name=part[column].str.split("|")).explode("name")
    unique = part.name.map(part.groupby("name").ensembl_gene_id.nunique()) == 1
    human_names.append(part[unique][["name", "ensembl_gene_id"]])
official = with_gene[["symbol", "ensembl_gene_id"]].set_axis(["name", "ensembl_gene_id"], axis=1)
human_names = pd.concat([*human_names, official]).set_axis(["name", "ensembl"], axis=1)
for species, current, synonyms, names in (
    ("human", human, pd.concat(human_old), human_names),
    ("mouse", mouse, mouse_old, mouse_names),
):
    synonyms = synonyms.assign(old=synonyms.old.str.split("|")).explode("old")
    synonyms = synonyms[(synonyms.old.fillna("") != "") & ~synonyms.old.isin(set(current.symbol))]
    upper = synonyms.old.str.upper()
    shared = upper.map(synonyms.assign(key=upper).groupby("key").symbol.nunique()) > 1
    synonyms = synonyms[~shared]
    by_ensembl = current[current.ensembl.fillna("") != ""].drop_duplicates("ensembl")
    steps = {
        "ensembl": (by_ensembl.ensembl, by_ensembl.symbol),
        "official": (current.symbol, current.symbol),
        "synonym": (synonyms.old, synonyms.symbol),
        "synonym_upper": (synonyms.old.str.upper(), synonyms.symbol),
        "to_ensembl": (names.name, names.ensembl),
    }
    if species == "human":
        steps["official_upper"] = (current.symbol.str.upper(), current.symbol)
    tables = [
        pd.DataFrame({"step": step, "query": query, "value": value}).drop_duplicates("query", keep="last")
        for step, (query, value) in steps.items()
    ]
    utils.write(pd.concat(tables).dropna().reset_index(drop=True), utils.PREPARED / f"ref/symbols_{species}.tsv")

# GO programmes: a term and its is_a, part_of descendants
# genes annotated to them, without NOT; myelin also takes 46 curated myelin and oligodendrocyte markers
roots = {"myelin": "GO:0043209", "lysosome": "GO:0005764", "mitochondrion": "GO:0005739", "synapse": "GO:0045202"}
myelin_markers = (
    "Mbp Plp1 Mog Mag Mobp Cnp Cldn11 Mal Ermn Gjb1 Opalin Tspan2 Aspa Ugt8a Gsn Qdpr Trf Apod Pllp Gpr37 Enpp6 Fa2h "
    "Sox10 Olig1 Olig2 Myrf Nkx2-2 Tppp Qk Hapln2 Anln Ndrg1 Sirt2 Cryab Bcas1 Cd9 Mal2 Gjc2 Gjb2 Plekhb1 Tmem10 "
    "Stmn4 Kif5a Tf Slc44a1 Ctps"
).split()
parents, term = {}, None
for line in open(RAW / "go_release/go-basic.obo"):
    if line.startswith(("[Term]", "is_obsolete: true")):
        term = None
    elif line.startswith("id: GO:"):
        term = line[4:].strip()
    elif term and line.startswith(("is_a: ", "relationship: part_of ")):
        parents.setdefault(term, set()).add(line.split()[1 if line.startswith("is_a") else 2])
children = {}
for child, group in parents.items():
    for parent in group:
        children.setdefault(parent, set()).add(child)
terms = {}
for name, root in roots.items():
    terms[name], stack = {root}, [root]
    while stack:
        for child in children.get(stack.pop(), ()):
            if child not in terms[name]:
                terms[name].add(child)
                stack.append(child)
rows = []
for gaf, species in (("mgi.gaf.gz", "Mus musculus"), ("goa_human.gaf.gz", "Homo sapiens")):
    symbols = {name: set() for name in roots}
    for line in gzip.open(RAW / "go_release" / gaf, "rt"):
        field = line.split("\t")
        if line[0] != "!" and "NOT" not in field[3]:
            for name, ids in terms.items():
                if field[4] in ids:
                    symbols[name].add(field[2])
    symbols["myelin"] |= {name.upper() if species == "Homo sapiens" else name for name in myelin_markers}
    assert all(len(group) > 50 for group in symbols.values())
    rows += [(species, name, symbol) for name, group in symbols.items() for symbol in sorted(group)]
utils.write(pd.DataFrame(rows, columns=["species", "programme", "symbol"]), OUT / "go_release/programmes.parquet")

# Protein half-lives at 5 and 21 months (Kluever 2022)
zipped = zipfile.ZipFile(RAW / "kluever2022_sciadv/sciadv.abn4437_tables_s1_to_s7.zip")
raw = pd.read_excel(
    io.BytesIO(zipped.read("sciadv.abn4437_table_s3.xlsx")), sheet_name="lifetimes and scaling", header=None
)
header, body = raw.iloc[0], raw.iloc[1:]
starts = [i for i, name in enumerate(header) if name == "protein ID"] + [raw.shape[1]]
pattern = r"(?P<scaled>SCALED_)?(?P<age>5m|21m)_(?P<tissue>Cer|Ctx)_(?P<fraction>\w+) (?:(?P<bound>c1|c2)|[Tt][Hh]alf)"
bound, blocks = {None: "t", "c1": "lo", "c2": "hi"}, []
for first, last in pairwise(starts):
    columns = {}
    for column in range(first + 1, last):
        found = re.match(pattern, str(header.iloc[column]))
        if found and found["fraction"] == "Hom":
            age = "21m_scaled" if found["scaled"] else found["age"]
            columns[f"{bound[found['bound']]}_{age}"] = pd.to_numeric(body.iloc[:, column])
            tissue = {"Cer": "cerebellum", "Ctx": "cortex"}[found["tissue"]]
    if columns:
        blocks.append(
            pd.DataFrame(columns).assign(accession=body.iloc[:, first].astype(str).str.split("-").str[0], tissue=tissue)
        )
life = pd.concat(blocks, ignore_index=True)
life = life[life[["t_5m", "t_21m", "t_21m_scaled"]].gt(0).all(axis=1)].reset_index(drop=True)
assert len(life) > 4000


def log2_se(half_life, low, high):
    # log2 half-life error from its 95% bounds; wide or missing bounds give 1, the result stays within 0.05-1
    half = (high - low) / 2
    se = half / utils.Z95 / (half_life * np.log(2))
    return np.clip(np.where((half >= 0.5 * half_life) | ~np.isfinite(se), 1.0, se), 0.05, 1.0)


se_young = log2_se(life.t_5m, life.lo_5m, life.hi_5m)
life["log2_t_young"] = np.log2(life.t_5m)
life["se_log2_t_young"] = se_young
for age, name in (("21m_scaled", "lfc_scaled"), ("21m", "lfc_unscaled")):
    life[name] = np.log2(life[f"t_{age}"] / life.t_5m)
    life[f"se_{name}"] = np.clip(
        np.hypot(se_young, log2_se(life[f"t_{age}"], life[f"lo_{age}"], life[f"hi_{age}"])), 0.05, 1.0
    )
utils.write(life, OUT / "kluever2022_sciadv/kluever_lifetimes.parquet")

# mRNA half-life PC1 (Saluki, Agarwal and Kelley 2022)
half = []
for species, name in (("Mus musculus", "mouse"), ("Homo sapiens", "human")):
    table = pd.read_csv(RAW / f"saluki_halflife/all_HLs_{name}_PC1.txt", sep="\t")
    half.append(
        pd.DataFrame(
            {
                "species": species,
                "gene_id": table.ensembl_gene_id,
                "symbol": table.external_gene_name,
                "hl_pc1": table.iloc[:, -1],
            }
        )
    )
utils.write(pd.concat(half, ignore_index=True), OUT / "saluki_halflife/half_life.parquet")

# Cell-type origin: share of pooled CPM per cell class
# adult Allen mouse atlas, Sloan human cortex
to_ensembl = utils.symbol_steps("human")["to_ensembl"]
classes = {
    "jin2025_zengagingmouse10xv3": ("mouse", "cell_class", ["Glut|GABA|Dopa|Sero", "^30 ", "^31 ", "^34 ", "^33 "]),
    "sloan2026_snrna": ("human", "cell_type", ["^(?:Exc|Inh)$", "^Ast$", "^(?:Oli|OPC)$", "^Mic$", "^(?:End|Peri)$"]),
}
origin = []
for dataset, (species, column, patterns) in classes.items():
    counts = utils.read(utils.PREPARED / dataset / "pseudobulk.parquet")
    samples = utils.read(utils.PREPARED / dataset / "pseudobulk_samples.tsv")
    features = utils.read(utils.PREPARED / dataset / "pseudobulk_features.tsv", index_col=0)
    gene = (features.ensembl if species == "mouse" else features.symbol.map(to_ensembl)).to_numpy()
    use = (samples.age_group == "adult") if species == "mouse" else pd.Series(True, index=samples.index)
    cpm = {}
    for name, regex in zip(["neuron", "astro", "oligo", "micro", "vascular"], patterns, strict=True):
        pooled = counts.loc[:, (use & samples[column].astype(str).str.contains(regex)).to_numpy()].sum(axis=1)
        cpm[name] = (pooled / pooled.sum() * 1e6).groupby(gene).sum()
    cpm = pd.DataFrame(cpm)
    total = cpm[["neuron", "astro", "oligo", "micro"]].sum(axis=1)
    share = cpm[["neuron", "astro", "oligo", "micro"]].div(total.replace(0, np.nan), axis=0).add_prefix("ct_frac_")
    origin.append(
        share.assign(ct_vascular_cpm=cpm.vascular, ct_total_cpm=total, species=species).rename_axis("gene_id")
    )
origin = pd.concat(origin)
utils.write(origin, OUT / "cell_type_origin/shares.parquet")

# GPU predictions of the release: Orthrus, RiboNN
for species in ("mouse", "human"):
    orthrus = np.load(RAW / f"gpu_precomputed/orthrus_{species}.npz")
    ids = pd.Index(orthrus["tx_id"].astype(str), name="tx_id")
    embedding = pd.DataFrame(orthrus["emb"], index=ids).add_prefix("emb").assign(gene_id=orthrus["gene_id"].astype(str))
    utils.write(embedding, OUT / f"orthrus_{species}/embeddings.parquet")
    te = pd.read_csv(RAW / f"gpu_precomputed/ribonn_{species}.tsv", sep="\t")
    utils.write(te, OUT / f"ribonn_{species}/predictions.parquet")

# Weights of published linear proteomic ageing clocks (benchmarks of script 52): the coefficient of each protein on its
# standardised level, summed over the assays of a gene; assays of several genes are dropped. Oh 2023: the mean over
# the 500 bootstrap LASSO models with zeros (ST7 averages only the bootstraps where the protein is selected);
# Goeminne 2025: the fold-1 coefficient on Olink NPX times the NPX standard deviation of training (Table S3), as the
# authors apply the models to other platforms; Melendez 2024: the elastic-net weights on min-max scaled levels.
# SomaScan symbols that are not current (old, ambiguous or mangled by Excel such as 3-Sep) go by their Entrez id
approved = set(hgnc.symbol[hgnc.status == "Approved"])
entrez = hgnc[(hgnc.entrez_id != "") & (hgnc.ensembl_gene_id != "")].drop_duplicates("entrez_id", keep=False)
entrez = entrez.set_index("entrez_id").ensembl_gene_id


def clock(label, source, fluid, assays, weights, entrez_ids=None):
    rows = pd.DataFrame({"assay": np.asarray(assays, dtype=object), "weight": np.asarray(weights, dtype=float)})
    rows["ensembl"] = np.nan if entrez_ids is None else pd.Series(np.asarray(entrez_ids, dtype=object)).map(entrez)
    rows = rows[rows.weight.fillna(0) != 0].fillna({"assay": ""})
    rows = rows[~rows.assay.str.contains(r"[ ,|._]")]
    rows["symbol"] = utils.current_symbol(rows.assay.to_numpy(), "human", np.full(len(rows), np.nan))
    by_entrez = utils.current_symbol(rows.assay.to_numpy(), "human", rows.ensembl.to_numpy())
    rows["symbol"] = rows.symbol.where(rows.symbol.isin(approved), by_entrez)
    rows = rows[rows.symbol.isin(approved)]
    genes = rows.groupby(rows.symbol.str.upper()).weight.agg(weight="sum", n_assays="size").reset_index()
    return genes.assign(clock=label, source=source, fluid=fluid)[["clock", "source", "fluid", *genes]]


clocks = []
st = RAW / "lehallier2019_clock/41591_2019_673_MOESM3_ESM.xlsx"
annotation = pd.read_excel(st, sheet_name="ST1 Nomenclature 2,925 proteins", header=2, dtype=str).set_index("ID")
lasso = pd.read_excel(st, sheet_name="ST7 Plasma proteomic clock", header=2)
lasso = lasso[lasso.Variable.isin(annotation.index)]
# symbols turned into dates by Excel (SEPT10) or left empty (APOE, EPHA10) are read from the assay name
symbols = annotation.EntrezGeneSymbol.mask(annotation.EntrezGeneSymbol.str.contains(":", na=False))
symbols = lasso.Variable.map(symbols).fillna(lasso.Variable.str.replace(r"(\.\d+){3}$", "", regex=True))
entrez_ids = lasso.Variable.map(annotation.EntrezGeneID)
clocks.append(clock("Lehallier 2019, plasma", "lehallier2019", "plasma", symbols, lasso.iloc[:, 1], entrez_ids))
st = RAW / "oh2023_clock/41586_2023_6802_MOESM5_ESM.xlsx"
annotation = pd.read_excel(st, sheet_name="ST3", dtype=str).set_index("Somamer")
mean = pd.read_excel(st, sheet_name="ST7")
for organ in ("Brain", "Organismal"):
    rows = mean[mean.organ == organ]
    symbols, entrez_ids = rows.Somamer.map(annotation.EntrezGeneSymbol), rows.Somamer.map(annotation.EntrezGeneID)
    weight = rows.mean_coef * rows.percboot_nonzero / 100
    clocks.append(clock(f"Oh 2023 {organ.lower()}, plasma", "oh2023", "plasma", symbols, weight, entrez_ids))
st = RAW / "melendez2024_clock/ACEL-23-e14230-s002.xlsx"
net = pd.read_excel(st, sheet_name="Full Model Protein Weights", dtype={"EntrezGeneID": str})
symbols, entrez_ids = net.EntrezGeneSymbol, net.EntrezGeneID
clocks.append(clock("Melendez 2024, CSF", "melendez2024", "CSF", symbols, net["Model Weight"], entrez_ids))
st = RAW / "oh2025_clock/41591_2025_3798_MOESM3_ESM.xlsx"
lasso = pd.read_excel(st, sheet_name="ST3", index_col="organ").drop(columns="intercept")
assays = lasso.columns.map(pd.read_excel(st, sheet_name="ST2").set_index("prot").Assay)
for organ in ("Brain", "Organismal"):
    clocks.append(clock(f"Oh 2025 {organ.lower()}, plasma", "oh2025", "plasma", assays, lasso.loc[organ]))
st = RAW / "goeminne2025_clock/1-s2.0-S1550413124004017-mmc"
net = pd.read_excel(f"{st}2.xlsx", sheet_name="Table S1A").query("Model == 'Brain' and Fold == 1")
net = net.drop(columns=["Model", "Fold", "Intercept"]).iloc[0]
sd = pd.read_excel(f"{st}4.xlsx", sheet_name="Table S3", index_col="Type").loc["Full models"]
clocks.append(clock("Goeminne 2025 brain, plasma", "goeminne2025", "plasma", net.index, net * sd.reindex(net.index)))
utils.write(pd.concat(clocks, ignore_index=True), utils.PREPARED / "ref/published_clocks.tsv")
