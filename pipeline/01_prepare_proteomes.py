# Protein-group matrices and sample tables of the proteome cohorts
# Paired cohorts also get RNA-seq sample tables, read by 02
# reads data/raw (00_download.sh); writes output/prepared/<dataset>/protein* and rna_samples.tsv
import gzip
import io
import json
import re

import numpy as np
import pandas as pd

import utils

RAW = utils.DATA / "raw"

# UniProt gene symbols; reviewed entries pick the lead accession
uniprot = pd.concat(
    pd.read_csv(RAW / f"uniprot/{s}_genes.tsv", sep="\t", dtype=str, keep_default_na=False) for s in ("mouse", "human")
)
uniprot = uniprot.drop_duplicates("Entry").set_index("Entry").rename(columns={"Gene Names (primary)": "gene"})
reviewed = uniprot.gene[uniprot.Reviewed == "reviewed"]
gene_of = uniprot.gene[uniprot.gene != ""]


def write_samples(dataset, layer, samples):
    # sample table with dates as text
    for column in samples.select_dtypes("datetime"):
        samples[column] = samples[column].dt.strftime("%Y-%m-%d")
    utils.write(samples, utils.PREPARED / dataset / f"{layer}_samples.tsv")


def save_protein(dataset, values, samples, groups, listed, names):
    # float32 matrix and features of detected groups
    features = pd.DataFrame(
        {
            "uniprot": groups.lead.str.split("-").str[0].to_numpy(),
            "gene_symbol": groups.gene.replace("", np.nan).to_numpy(),
            "genes": pd.Series(np.asarray(listed)).replace("", np.nan).to_numpy(),
        },
        index=pd.Index(names, name="protein_id"),
    )
    features["gene_symbol"] = features.gene_symbol.fillna(features.uniprot.map(gene_of))
    features["genes"] = features.genes.fillna(features.gene_symbol)
    matrix = values.set_axis(features.index).set_axis(samples.sample_id.to_numpy(), axis=1)
    seen = matrix.notna().any(axis=1).to_numpy()
    assert matrix.columns.is_unique and features.index.is_unique
    utils.write(matrix[seen].astype("float32"), utils.PREPARED / dataset / "protein.parquet")
    utils.write(features[seen], utils.PREPARED / dataset / "protein_features.tsv")
    write_samples(dataset, "protein", samples)


def save_groups(dataset, values, samples, accessions, genes, lead="first", sep=";", missing="", split=None):
    # lead accession and gene of each protein group
    # lead: first accession, or the first reviewed one whose gene is
    # the first (reviewed_first_gene) or any (reviewed_any_gene) listed gene
    # missing: name part of a group without a gene, None for no part
    # split: accessions and genes in one text, kept as the group name
    texts = accessions
    if split:
        accessions, genes = (texts.str.split(split, n=1).str[i] for i in (0, 1))
        genes = genes.str.rstrip(";")
    groups = []
    for accession, text in zip(accessions, genes, strict=True):
        parts = str(accession).split(sep)
        names = [g for g in str(text).split(";") if g and g != "nan"]
        chosen = parts[0]
        if lead != "first":
            listed = list(dict.fromkeys(a.split("-")[0] for a in parts))
            known = [a for a in listed if a in reviewed.index]
            allowed = names if lead == "reviewed_any_gene" else names[:1]
            chosen = next((a for a in known if reviewed[a] and reviewed[a] in allowed), None) or (known or listed)[0]
        own = reviewed.get(chosen.split("-")[0], "")
        groups.append((chosen, own if lead == "reviewed_any_gene" and own in names else (names or [""])[0]))
    groups = pd.DataFrame(groups, columns=["lead", "gene"])
    # group names <accession>_<gene>, repeats numbered by detection
    if missing is None:
        base = groups.lead.where(groups.gene == "", groups.lead + "_" + groups.gene)
    else:
        base = groups.lead + "_" + groups.gene.where(groups.gene != "", missing)
    rank = pd.Series(values.notna().sum(axis=1).rank(ascending=False, method="first").to_numpy())
    repeat = rank.groupby(base).rank(method="first").astype(int)
    names = texts if split else base.where(repeat == 1, base + "~" + repeat.astype(str))
    save_protein(dataset, values, samples, groups, (genes if split else genes.fillna("")).to_numpy(), names.to_numpy())


def match(columns, pattern):
    # samples from the columns matching a pattern
    columns = pd.Series(list(columns))
    samples = pd.DataFrame({"column": columns[columns.str.match(pattern)].to_numpy()})
    return samples.join(samples.column.str.extract(pattern)) if re.compile(pattern).groups else samples


# Takasugi 2024 cerebrum: TMT proteome and RNA-seq
dataset = "takasugi2024_pxd040722"
table = pd.read_excel(RAW / dataset / "supp/41467_2024_52845_MOESM5_ESM.xlsb", sheet_name="Brain_WTL", engine="pyxlsb")
table = table.drop_duplicates("Gene name").reset_index(drop=True)
# RNA columns of the GEO table
rna_columns = pd.read_csv(RAW / dataset / "geo/GSE225576_Brain.txt.gz", sep="\t", nrows=0).columns
# protein and RNA samples from the column names
for layer, columns, pattern, batch in (
    ("protein", table.columns, r"^(?P<age_months>\d+)mo-(?P<replicate>\d+)$", "TMT02_Brain_WTL"),
    ("rna", rna_columns, r"^(?P<age_months>\d+)M(?P<replicate>\d+)$", "single"),
):
    samples = match(columns, pattern).assign(sheet="Brain_WTL", tissue="Cerebrum", sex="male", batch=batch)
    samples["sample_id"] = "Brain_WTL_" + samples.age_months + "mo_" + samples.replicate
    samples["donor_id"] = samples.age_months + "mo_rep" + samples.replicate
    samples["age_months"] = samples.age_months.astype(int)
    assert len(samples) == 16
    if layer == "rna":
        write_samples(dataset, layer, samples.drop(columns="sheet"))
    else:
        save_groups(
            dataset,
            np.log2(table[samples.column].astype(float)),
            samples,
            table["Protein ID"],
            table["Gene name"],
            sep="|",
        )

# Langfelder 2016 striatum: MaxQuant workbook (log10)
# flagged groups (reverse, contaminant, site-only) come from proteinGroups
dataset = "langfelder2016_pxd006302"
GENOTYPE = {
    "Q0WT": "WT",
    "Q020": "Q20",
    "Q050": "Q50",
    "Q080": "Q80",
    "Q092": "Q92",
    "Wild Type": "WT",
    "Het (Q20) Knock-In": "Q20",
    "Het (Q80) Knock-In": "Q80",
    "Het (Q92) Knock-In": "Q92",
    "Het (Q111) Knock-In": "Q111",
    "Het (Q140) Knock-In": "Q140",
    "Het (Q175) Knock-In": "Q175",
}
AGE = {"02M": 2, "06M": 6, "10M": 10, "2M": 2, "6M": 6, "2 month": 2, "6 month": 6, "10 month": 10}
book = pd.read_excel(
    gzip.open(RAW / dataset / "hdinhd/KIQ_144062_Evotec_Proteome_Part2_WP10_reprocessed_Striatum_Results.xlsx.gz"),
    header=1,
)
meta = pd.read_excel(
    gzip.open(RAW / dataset / "hdinhd/KIQ-Brain-Samples-160725-Animal-to-Sample-mapping.xlsx.gz"), header=1
)
# samples from the workbook columns, animals from the sheet
samples = match(book.columns, r"^(?:LFQIntensity)?(?P<experiment>(?:Q\d+|Q0WT)_[FM]_\d+M_[A-Za-z]+_R\d_\d+)_?$")
samples = samples.merge(meta, left_on="experiment", right_on="Experimental Design Name", how="left")
samples["tissue"] = "Striatum"
samples["donor_id"] = samples.PGIDNumber
samples["genotype"] = samples["Factor Genotype"].replace(GENOTYPE)
samples["age_months"] = pd.to_numeric(samples["Factor Age"].replace(AGE), errors="coerce")
samples["sex"] = samples["Factor Sex"].replace(utils.SEX)
samples["batch"] = samples["Work-Package"]
samples["arm"] = np.where(samples.genotype == "WT", "aging", "disease")
samples["sample_id"] = "Str_" + samples.genotype + "_" + samples.donor_id.astype(str)
# intensities, zero means not detected
values = book[samples.column].apply(pd.to_numeric, errors="coerce")
values = values.where(values > 0)
info = pd.DataFrame(
    {"id": pd.to_numeric(book["Id"], errors="coerce"), "proteins": book["Uniprot Ids"].fillna("").astype(str)}
)
# groups without ids and intensities are dropped
blank = info.proteins.isin(["", "nan"]) & values.isna().all(axis=1)
info, values = info[~blank], values[~blank]
# proteinGroups gives flags, majority ids and gene names
flags = ["Reverse", "Potential contaminant", "Only identified by site"]
columns = ["id", "Protein IDs", "Majority protein IDs", "Gene names", "Razor + unique peptides", *flags]
run = pd.read_csv(
    RAW / dataset / "maxquant/striatum_PXD006302/proteinGroups.txt.gz", sep="\t", usecols=columns, low_memory=False
).set_index("id")
# groups without an id are found by their protein ids
lost = info.id.isna()
info.loc[lost, "id"] = info.proteins[lost].map(pd.Series(run.index, index=run["Protein IDs"]))
run = run.loc[info.id.astype(int)]
# flagged and undetected groups leave; log10 to log2
keep = ~run[flags].eq("+").any(axis=1).to_numpy() & values.notna().any(axis=1).to_numpy()
run, values = run[keep], values[keep] * np.log2(10)
# name groups by detection (most detected first); a repeated name gets a ~2, ~3 suffix
taken, named = set(), {}
order = np.lexsort((-run["Razor + unique peptides"].astype(float).to_numpy(), -values.notna().sum(axis=1).to_numpy()))
for row in order:
    majority = [accession for accession in str(run["Majority protein IDs"].iloc[row]).split(";") if accession]
    majority = majority or run["Protein IDs"].iloc[row].split(";")
    gene = run["Gene names"].iloc[row]
    gene = str(gene).split(";")[0] if pd.notna(gene) else ""
    known = [accession for accession in (name.split("-")[0] for name in majority) if accession in reviewed.index]
    lead = next((a for a in known if gene and reviewed[a].lower() == gene.lower()), None) or (known or majority)[0]
    base = f"{lead}_{gene}" if gene else lead
    name, repeat = base, 1
    while name in taken:
        repeat += 1
        name = f"{base}~{repeat}"
    taken.add(name)
    named[row] = (name, lead, gene)
names, leads, genes = zip(*(named[row] for row in range(len(run))), strict=True)
groups = pd.DataFrame({"lead": leads, "gene": genes})
save_protein(dataset, values, samples, groups, groups.gene.to_numpy(), list(names))

# RNA samples: decoder ring maps columns to animals
decoder = pd.read_excel(
    gzip.open(RAW / dataset / "hdinhd/Allelic_Series_Decoder_Ring-1.5.xlsx.gz"),
    sheet_name="Allelic Series Decoder Ring ",
)
rna = pd.DataFrame(
    {
        "column": pd.read_csv(
            RAW / dataset / "geo/GSE65774_Striatum_mRNA_counts_processedData.txt.gz", sep="\t", nrows=0
        ).columns[1:]
    }
)
rna = rna.merge(decoder, left_on="column", right_on="Sample Name Column Header in Count file", how="left")
rna["tissue"] = "Striatum"
rna["donor_id"] = rna["Animal ID"]
rna["genotype"] = rna.Genotype.replace(GENOTYPE)
rna["age_months"] = pd.to_numeric(rna.Age.replace(AGE), errors="coerce")
rna["sex"] = rna.Sex.replace(utils.SEX)
rna["batch"] = rna.age_months.replace({2: "lib_2m10m", 6: "lib_6m", 10: "lib_2m10m"})
rna["arm"] = np.where(rna.genotype == "WT", "aging", "disease")
rna["sample_id"] = "Str_" + rna.genotype + "_" + rna.donor_id.astype(str)
assert samples.donor_id.notna().all() and rna.donor_id.notna().all()
write_samples(dataset, "rna", rna)

# Luo 2024 frontal cortex: TMT proteome (S2), RNA-seq (S3)
# 2016-957 is 2015-4 in S1B, 2016-1532 is 2015-1532 in S1D
dataset = "luo2024_syn32136022"
book = RAW / dataset / "supp/41380_2024_2576_MOESM2_ESM.xlsx"
sheets = pd.read_excel(book, sheet_name=["Table S2", "Table S1B", "Table S1D"], header=2)
rna_columns = pd.read_excel(book, sheet_name="Table S3", header=2, nrows=0, index_col="Ensembl ID").columns
table = sheets["Table S2"]
for layer, columns, meta, key, rename, batch in (
    ("protein", table.columns, sheets["Table S1B"], "Samplel ID", {"2016-957": "2015-4"}, "Batch"),
    ("rna", rna_columns, sheets["Table S1D"], "Sample ID", {"2016-1532": "2015-1532"}, "RNA library"),
):
    samples = match(columns, r"^\d{4}-\d+$")
    # sample sheet by donor id
    samples = samples.merge(meta, left_on=samples.column.replace(rename), right_on=key, how="left")
    samples = samples.assign(
        **({"sheet": "Table S2"} if layer == "protein" else {}),
        tissue="Frontal cortex",
        donor_id=samples.column,
        sample_id=samples.column,
    )
    samples["age_years"] = pd.to_numeric(samples["Age of death"].replace("90+", 90), errors="coerce")
    samples["sex"] = samples["Sex M=1, F=2"].map({1: "male", 2: "female"})
    samples["batch"] = samples[batch]
    samples["brain_bank"] = samples["Brain bank"].replace({"American": None})
    samples["pmi_h"] = pd.to_numeric(samples.PMI, errors="coerce")
    if layer == "rna":
        samples["rin"] = pd.to_numeric(samples.RIN, errors="coerce")
    samples["arm"] = np.where(samples.Group == "CTRL", "aging", "disease")
    assert samples.Group.notna().all()
    if layer == "rna":
        write_samples(dataset, layer, samples)
    else:
        save_groups(
            dataset,
            table[samples.column].astype(float),
            samples,
            table["Protein accession"],
            table["Gene symbol"],
            missing="NA",
        )

# Koopmans 2026 cortex layers (DIA-NN), runs matched to sheet
dataset = "koopmans2026_pxd058441"
table = pd.read_csv(RAW / dataset / "diann/20220120_rklaas_SUN-schizo-cortex-layer1-3+4-6_MBR.pg_matrix.tsv", sep="\t")
samples = match(table.columns, r"^.*\\(?P<run>[^\\]+)$")
samples = samples.merge(
    pd.read_excel(RAW / dataset / "supp/MOESM3_supplementary_data_1_samples.xlsx", sheet_name="data"),
    left_on="run",
    right_on="sample_id",
    how="left",
)
samples["sheet"] = None
samples["sample_id"] = samples.shortname
samples["donor_id"] = samples.case_id
samples["age_years"] = samples.age
assert len(samples) == 192 and samples.shortname.notna().all()
samples["batch"] = "gel_day" + samples.gel_preparation_day.astype(str)
samples["tissue"] = "DPFC " + samples.cortical_layers
samples["cohort"] = samples.diagnosis
samples["pmd_h"] = samples.postmortem_delay
samples["brain_bank"] = samples.origin
samples["sex"] = samples.sex.replace(utils.SEX)
save_groups(
    dataset,
    np.log2(table[samples.column].astype(float)),
    samples,
    table["Protein.Group"],
    table["Genes"],
    lead="reviewed_first_gene",
)

# Eggers 2025 substantia nigra (Spectronaut, log2)
dataset = "eggers2025_pxd051145"
table = pd.read_excel(RAW / dataset / "pride/norm_proteins.xlsx")
samples = match(table.columns, r"^(?:young|middle|old)_\d+$")
samples = samples.merge(
    pd.read_excel(RAW / dataset / "pride/sample_information.xlsx"),
    left_on="column",
    right_on="Sample Name ",
    how="left",
)
samples = samples.assign(
    sheet=None,
    tissue="Substantia nigra",
    batch="spectronaut16_directdia",
    donor_id=samples.column,
    sample_id=samples.column,
)
samples["age_years"] = pd.to_numeric(samples.AGE, errors="coerce")
samples["pmi_h"] = pd.to_numeric(samples.PMI, errors="coerce")
samples["ms_id"] = samples["MS ID"]
samples["sex"] = samples.SEX.replace(utils.SEX)
samples["ms_run_order"] = samples.ms_id.rank(method="first").astype(int)
save_groups(
    dataset,
    table[samples.column].astype(float),
    samples,
    table["PG.ProteinAccessions"],
    table["PG.Genes"],
    lead="reviewed_any_gene",
    missing=None,
)

# Anderton 2026 hippocampus (Spectronaut), unidentified or floor values missing
dataset = "anderton2026_biorxiv714902"
table = pd.read_excel(
    RAW / dataset / "github/data_raw/EA6_lysate_SN17_2pep_Analysis_V02.xlsx",
    sheet_name="Identified Protein Groups",
    header=3,
)
pattern = (
    r"^X\.\d+\.\.(?P<yy>\d{2})_(?P<month>\d{2})(?P<day>\d{2})_(?P<run_no>\d{3})_(?P<ms_batch>EA6|EA5C)"
    r"_(?P<position>\d{2})_(?P<author_group>[YMOG])_(?P<run_donor>[A-Za-z0-9]+)_TH_SW_01\.htrms\.PG\.Quantity$"
)
samples = match(table.columns, pattern)
key = ("S" + samples.run_donor).replace({"SS12676": "SS12376", "S3823": "Sm3823"})
meta = pd.read_csv(RAW / dataset / "github/data_raw/final_metadata.csv", encoding="utf-8-sig", na_values=["#N/A"])
samples = samples.merge(meta, left_on=key, right_on="individual", how="left")
samples = samples.assign(
    sheet="Identified Protein Groups",
    tissue="Hippocampus",
    sex="female",
    sample_id=samples.individual,
    donor_id=samples.individual,
)
samples["batch"] = samples.ms_batch + "_20" + samples.yy + "-" + samples.month
samples["age_years"] = pd.to_numeric(samples.age, errors="coerce")
samples["pmi_h"] = pd.to_numeric(samples.PMI, errors="coerce")
counts = table[[c.replace("PG.Quantity", "PG.NrOfPrecursorsIdentified") for c in samples.column]].set_axis(
    samples.column, axis=1
)
values = table[samples.column].astype(float)
save_groups(
    dataset, np.log2(values.where((counts != 0) & (values > 2.0))), samples, table["PG.UniProtIds"], table["PG.Genes"]
)

# Lu 2025 mouse brain regions (PulseDIA, log2), aging and diet arms (the units use the aging arm)
dataset = "lu2025_pxd053154"
table = pd.read_csv(RAW / dataset / "gdrive/Aging_mouse_organ_PulseDIA_Matrix.tsv.xz", sep="\t")
samples = match(table.columns, r"^(?P<run>.+_(?P<plate>b\d+))_\d+$")
samples = samples.merge(
    pd.read_excel(RAW / dataset / "gdrive/Aging_mouse_organ_SampleInformation.xlsx"),
    left_on="column",
    right_on="Sample_name",
    how="left",
)
samples["sheet"] = None
samples["animal"] = samples.id.astype(str).str.extract(r"^([^_]+)")[0]
samples["sample_id"] = samples.column
samples["donor_id"] = samples.Group + "_" + samples.animal  # group A1: aging arm; G1-G4: diet groups
samples["tissue"] = samples.Organ
samples["age_months"] = pd.to_numeric(samples.Month, errors="coerce")
samples["sex"] = samples.Sex.replace(utils.SEX)
samples["arm"] = np.where(samples.Type == "Aging", "aging", "diet")
samples["batch"] = samples.run
samples = samples[samples.Organ.isin(["Brainstem", "Cerebellum", "Hippocampus", "Olfactory bulb"])]
save_groups(dataset, table[samples.column].astype(float), samples, table["...1"], table["...1"], split="_")

# Tsumagari 2023 cortex, hippocampus (TMT), joined on protein id
dataset = "tsumagari2023_pxd041485"
parts = [
    pd.read_excel(RAW / dataset / f"supp/41598_2023_45570_MOESM{n}_ESM.xlsx", sheet_name="expression")
    .drop_duplicates("id")
    .set_index("id", drop=False)
    for n in (3, 4)
]
samples = pd.concat(match(p.columns, r"^(?P<code>Cx|Hipp)_(?P<age_months>\d+)mo_(?P<replicate>\d+)$") for p in parts)
samples = samples.assign(
    sheet="expression", sex="male", sample_id=samples.column, donor_id=samples.age_months + "mo_" + samples.replicate
)
samples["tissue"] = samples.code.replace({"Cx": "Cortex", "Hipp": "Hippocampus"})
samples["plex"] = samples.replicate.replace({"1": 1, "2": 1, "3": 1, "4": 2, "5": 2, "6": 2})
samples["batch"] = samples.code + "-" + samples.plex.astype(str)
samples["age_months"] = samples.age_months.astype(int)
values = pd.concat([p[[c for c in samples.column if c in p.columns]].astype(float) for p in parts], axis=1)
annotation = (
    pd.concat([p[["UniProt accession", "Gene name"]] for p in parts])
    .groupby(level=0, sort=False)
    .first()
    .loc[values.index]
)
save_groups(dataset, values, samples.reset_index(drop=True), annotation["UniProt accession"], annotation["Gene name"])

# Wang 2025 brain (DIA-NN), runs 81-120 as in Table S1
dataset = "wang2025_pxd047296"
table = pd.read_csv(RAW / dataset / "pride/BrainAging_report.pg_matrix.tsv", sep="\t", quoting=3)
samples = match(table.columns, r"^.*\\(?P<run>\d+)\.mzML$")
order = pd.MultiIndex.from_product(
    [[4, 8, 12, 20], ["M", "F"], range(1, 6)], names=["age_months", "sex", "replicate"]
).to_frame(index=False)
samples = samples.join(order).assign(sheet=None, tissue="Brain", batch="diann_brain_search")
samples["sample_id"] = "brain" + samples.age_months.astype(str) + samples.sex + samples.replicate.astype(str)
samples["donor_id"] = samples.sample_id.str[5:]
assert len(samples) == 40
samples["sex"] = samples.sex.replace(utils.SEX)
save_groups(
    dataset,
    np.log2(table[samples.column].astype(float)),
    samples,
    table["Protein.Group"],
    table["Genes"],
    lead="reviewed_first_gene",
)

# Keele 2025 cortex, hippocampus, striatum (TMT), joined on UniProt
dataset = "keele2025_pxd058160"
parts, frames = [], []
for sheet, short in (("Cortex", "Ctx"), ("Hippocampus", "Hip"), ("Striatum", "Str")):
    part = (
        pd.read_excel(RAW / dataset / "figshare/table_s3.xlsx", sheet_name=sheet)
        .drop_duplicates("uniprot_id")
        .set_index("uniprot_id", drop=False)
    )
    frame = match(part.columns, r"^tag_\w+$").merge(
        pd.read_excel(RAW / dataset / "figshare/table_s1.xlsx", sheet_name=sheet, dtype={"tmt_tag": str}),
        left_on="column",
        right_on="tmt_tag",
        how="left",
    )
    frame = frame.assign(sheet=sheet, tissue=sheet, short=short)
    frame["sample_id"] = frame.short + "_" + frame.mouse_id
    frame["donor_id"] = frame.mouse_id
    frame["age_months"] = pd.to_numeric(frame.age, errors="coerce")
    frame["sex"] = frame.sex.replace(utils.SEX)
    values = part[frame.column].astype(float)
    parts.append(np.log2(values.where(values != 0)))
    frames.append((frame, part[["uniprot_id", "symbol"]]))
samples = pd.concat(f for f, _ in frames)
annotation = pd.concat(a for _, a in frames).groupby(level=0, sort=False).first()
values = pd.concat(parts, axis=1)
save_groups(
    dataset,
    values,
    samples.reset_index(drop=True),
    annotation.loc[values.index, "uniprot_id"],
    annotation.loc[values.index, "symbol"],
    missing=None,
)

# Marino 2025 whole brain (ubiquitin-enriched TMT10 PSMs), protein medians
dataset = "marino2025_pxd047951"
channels = ["126", "127_N", "127_C", "128_N", "128_C", "129_N", "129_C", "130_N", "130_C", "131"]
psm = pd.read_csv(
    RAW / dataset / "massive/200904_AM_Brain_Mouse_UbEnr_TMT10_fr1_TargetPeptideSpectrumMatch.txt",
    sep="\t",
    usecols=["Annotated Sequence", "Master Protein Accessions", "Quan Info", "Ions Score", *channels],
    low_memory=False,
)
master = psm["Master Protein Accessions"].fillna("")
intensity = psm[channels].astype(float)
intensity = intensity.where(intensity > 0)
# unique PSMs of one protein, no contaminants, ion score of at least 15
keep = (
    ~master.str.contains(";")
    & ~master.str.contains("CON_")
    & (master != "")
    & (psm["Quan Info"] == "Unique")
    & (psm["Ions Score"] >= 15)
)
# a channel of 1000 or more, at most one channel missing
keep &= (intensity >= 1e3).any(axis=1) & (intensity.isna().sum(axis=1) <= 1)
psm, intensity = psm[keep], np.log2(intensity[keep])
# channel medians aligned to their common median
centre = intensity.median()
intensity = intensity - centre + centre.median()
# protein medians over PSMs, proteins with at least two peptides
protein = psm["Master Protein Accessions"]
peptides = psm["Annotated Sequence"].str.upper().groupby(protein).nunique()
values = intensity.groupby(protein).median()[peptides >= 2]
accessions = values.index.to_series()
genes = accessions.str.split("-").str[0].map(gene_of).fillna("")
mice = pd.DataFrame(
    {
        "channel": channels,
        "sample": [f"Young{i}" for i in range(1, 6)] + [f"Old{i}" for i in range(1, 6)],
        "animal": ["27", "28", "137", "139", "140", "303", "305", "306", "307", "308"],
    }
)
samples = pd.DataFrame({"column": channels}).merge(
    mice.assign(group=["young"] * 5 + ["old"] * 5), left_on="column", right_on="channel", how="left"
)
samples = samples.assign(
    tissue="Whole brain", sex="male", batch="tmt10_plex1", sample_id=samples["sample"], donor_id=samples.animal
)
samples["age_months"] = samples.group.map({"young": 3.5, "old": 33.0})
samples["tmt_channel"] = samples.channel
save_protein(
    dataset,
    values[channels],
    samples,
    pd.DataFrame({"lead": accessions.to_numpy(), "gene": genes.to_numpy()}),
    genes.to_numpy(),
    (accessions + "_" + genes).to_numpy(),
)

# Lee 2023 hippocampus CA1 (Spectronaut), neurologically normal donors
# quantities: Spectronaut protein-group quantities of the PRIDE report (authors' normalised abundances)
# authors dropped donors 12 (Braak III) and 57 (outlier), kept here and marked
dataset = "lee2023_pxd043880"
table = pd.read_csv(RAW / dataset / "pride/JL_EXP3_300819_Upload.xls", sep="\t")
# donor runs, age, sex, PMI and Braak stage
donors = pd.read_excel(
    RAW / dataset / "supp/13024_2023_650_MOESM1_ESM.xlsx",
    sheet_name="Source Data Proteins",
    header=1,
    usecols=["Sample ID", "FileName (Proteomics)", "Sex", "Age", "PMI", "Braak stage"],
    nrows=76,
)
# brain bank and pH from table S1
bank = pd.read_excel(
    RAW / dataset / "supp/13024_2023_650_MOESM2_ESM.xlsx",
    sheet_name="S1_Sample_Information",
    header=1,
    usecols=["Sample ID", "Brain Bank", "Brain pH", "Age"],
    nrows=76,
)
for sheet in (donors, bank):
    sheet["number"] = sheet["Sample ID"].astype(str).str.extract(r"^(\d+)")[0]
donors = donors.merge(bank, on="number", suffixes=("", "_s1"), validate="1:1")
assert (donors.Age == donors.Age_s1).all()
donors = pd.DataFrame(
    {
        "run": donors["FileName (Proteomics)"].str.extract(r"_EXP3_(\d+)\.raw$")[0],
        "sex": donors.Sex.replace(utils.SEX),
        "age_years": donors.Age,
        "pmi_h": donors.PMI,
        "braak": donors["Braak stage"],
        "brain_bank": donors["Brain Bank"],
        "brain_ph": donors["Brain pH"],
        "author_excluded": donors["Sample ID"].astype(str).str.extract(r"\((\w+)\)$")[0].str.lower(),
    }
)
# donor runs only, others are replicates or failed
samples = match(table.columns, r"^\[\d+\] .*_EXP3_(?P<run>\d+)\.raw\.PG\.Quantity$").merge(donors, on="run")
samples = samples.assign(tissue="Hippocampus CA1", arm="aging", diagnosis="neurologically normal", batch="EXP3")
samples["sample_id"] = "EXP3_" + samples.run
samples["donor_id"] = samples.sample_id
samples["ms_run_order"] = samples.pop("run").astype(int)
assert len(samples) == 76 and samples.donor_id.is_unique and samples.age_years.notna().all()
assert samples.author_excluded.value_counts().to_dict() == {"removed": 1, "outlier": 1}
assert samples.brain_bank.value_counts().to_dict() == {"NSWBTRC": 38, "QBB": 38}
# quantities are linear, missing when not quantified
save_groups(
    dataset,
    np.log2(table[samples.column].astype(float)),
    samples,
    table["PG.ProteinAccessions"],
    table["PG.Genes"],
    lead="reviewed_any_gene",
    missing=None,
)

# Teng 2023 six brain regions (TMTpro PSMs)
# quantities: log2 ratios to the 126 pool, protein medians (Lee 2020 pipeline), rebuilt from PRIDE PSMs
dataset = "teng2023_pxd040884"
CHANNELS = "126 127N 127C 128N 128C 129N 129C 130N 130C 131N 131C 132N 132C 133N 133C 134N".split()
POOLS = ["126", "127N"]
PSM_COLUMNS = [
    "# Proteins",
    "Master Protein Accessions",
    "Master Protein Descriptions",
    "Protein Accessions",
    "Isolation Interference [%]",
]
CHANNEL_MAP = {
    "Protein Sample / Subject ID": "individual",
    "Subcollection": "tissue",
    "Sample Type": "group",
    "Type of Replicate": "replicate",
    "TMT Sample Plex Number": "plex",
    "TMT Reporter Ion": "tmt_channel",
    "Associated Search Result File": "psm_file",
}
CONTAMINANT = r"^(?:P00761|IRT\d+)$"
MAX_INTERFERENCE = 50  # isolation interference, percent
MIN_PSMS = 2  # PSMs per protein
REGION = {
    "Naccumbens": "Nucleus accumbens",
    "Prefrontal Cortex": "Prefrontal cortex",
    "Ventral Tegmental Area": "Ventral tegmental area",
}
# channel map and PSM table of each plex
reference = pd.read_excel(RAW / dataset / "pride/Brain_AUD_Sample_Reference_File.xlsx", usecols=list(CHANNEL_MAP))
reference = reference.rename(columns=CHANNEL_MAP)
reference["batch"] = reference.plex.map("P{:02d}".format)
reference["tmt_channel"] = reference.tmt_channel.str.replace("Abundance: ", "")
reference["tissue"] = reference.tissue.replace(REGION)
# age, PMI and pH of the 27 individuals
people = pd.read_excel(
    RAW / dataset / "supp/41398_2023_2605_MOESM3_ESM.xlsx",
    header=None,
    skiprows=3,
    nrows=27,
    usecols=range(5),
    names=["individual", "diagnosis", "age_years", "pmi_h", "brain_ph"],
)
# individual channels without the two pools
samples = reference[~reference.tmt_channel.isin(POOLS)].astype({"individual": int})
samples = samples.merge(people, on="individual", how="left")
samples["column"] = samples.batch + "_" + samples.tmt_channel
samples["donor_id"] = "S" + samples.individual.astype(str)
samples["sample_id"] = samples.column + "_" + samples.donor_id
samples["sex"] = "male"
samples["arm"] = np.where(samples.diagnosis == "Control", "aging", "disease")
# technical replicates are individuals in both plexes
samples["technical_replicate"] = samples.replicate == "Technical"
assert len(samples) == 168
assert (samples.group == samples.diagnosis).all()
assert (samples.technical_replicate == samples.duplicated(["tissue", "individual"], keep=False)).all()
samples = samples.drop(columns=["individual", "group", "replicate", "plex", "psm_file"])


def mode(ratios):
    # peak of the smoothed histogram of one channel
    counts, edges = np.histogram(ratios.dropna(), bins=np.arange(-4, 4.01, 0.02))
    return edges[np.convolve(counts, np.hanning(21), mode="same").argmax()] + 0.01


# PSM ratios to the pool, proteins from their PSMs
reporters = [f"Abundance: {channel}" for channel in CHANNELS]
parts, levels, descriptions = [], [], []
for batch, file in reference.dropna(subset="psm_file").set_index("batch").psm_file.items():
    psm = pd.read_csv(RAW / dataset / f"pride/{file}.txt", sep="\t", usecols=PSM_COLUMNS + reporters, low_memory=False)
    intensity = psm[reporters].set_axis(CHANNELS, axis=1)
    intensity = intensity.where(intensity > 0)
    master = psm["Master Protein Accessions"].fillna("")
    # pool signal, low isolation interference, no contaminants
    keep = intensity["126"].notna() & (psm["Isolation Interference [%]"] < MAX_INTERFERENCE)
    keep &= (master != "") & ~master.str.match(CONTAMINANT)
    psm, intensity, master = psm[keep], intensity[keep], master[keep]
    # log2 ratios to the pool, channels centred on their mode
    ratio = np.log2(intensity.drop(columns=POOLS).div(intensity["126"], axis=0))
    ratio = ratio - ratio.apply(mode)
    # PSMs of several proteins join the protein with the closest ratios
    single = psm["# Proteins"] == 1
    profile = ratio[single].groupby(master[single]).median()
    shared = psm.loc[~single & (ratio.notna().mean(axis=1) >= 0.5), "Protein Accessions"].str.split("; ").explode()
    shared = shared[shared.isin(profile.index)]
    gap = np.nanmean((ratio.loc[shared.index].to_numpy() - profile.loc[shared].to_numpy()) ** 2, axis=1)
    options = pd.DataFrame({"psm": shared.index, "protein": shared.to_numpy(), "gap": gap})
    closest = options.loc[options.groupby("psm").gap.idxmin()]
    protein = master.where(single)
    protein.loc[closest.psm.to_numpy()] = closest.protein.to_numpy()
    # medians of proteins with at least two PSMs
    count = ratio.notna().any(axis=1).groupby(protein).sum()
    parts.append(ratio.groupby(protein).median()[count >= MIN_PSMS].add_prefix(f"{batch}_"))
    levels.append(np.log2(intensity["126"].groupby(protein).sum())[count >= MIN_PSMS])
    descriptions.append(psm[single].groupby(master[single])["Master Protein Descriptions"].first())
# pool intensity as the abundance level of a protein, the same in all plexes
level = pd.concat(levels, axis=1).mean(axis=1)
values = pd.concat(parts, axis=1).add(level, axis=0).sort_index()[samples.column]
assert values.notna().sum().min() > 4000
# gene names from the protein descriptions
descriptions = pd.concat(descriptions)
genes = descriptions[~descriptions.index.duplicated()].reindex(values.index).str.extract(r"GN=(\S+)")[0]
save_groups(dataset, values, samples, values.index.to_series(), genes, missing=None)

# Merrihew 2023 SMTG, hippocampus, IPL, caudate (DIA, batch-adjusted log2), joined on protein group
# quantities: authors' Level 3B protein groups, normalised and batch-adjusted (Panorama batchadj tables)
dataset = "merrihew2023_pxd034525"
parts, frames = [], []
for folder, prefix, batches, tissue in (
    ("SMTG", "SMTG", "SMTG-B1-B5", "Superior and middle temporal gyri (SMTG)"),
    ("Hippocampus", "Hipp", "Hipp-B1-B3", "Hippocampus"),
    ("IPL", "IPL", "IPL-B1-B5", "Inferior parietal lobule (IPL)"),
    ("Caudate", "Caudate", "Caudate-B1-B4", "Caudate"),
):
    table = pd.read_csv(RAW / dataset / f"panorama/{folder}/level3b/{prefix}-protGrp.batchadj.csv", index_col=0)
    meta = pd.read_csv(
        RAW / dataset / f"panorama/{folder}/metadata/HuAD_Clean_Diagnosis_Samples-{batches}-Meta.csv",
        usecols=lambda name: not name.startswith("Unnamed"),
    )
    # donor runs only, pooled references are no samples
    frame = match(table.columns, r"^[A-Z]ZR\d+$").merge(meta, left_on="column", right_on="Sample Label", how="left")
    frame["tissue"] = tissue
    frame["batch"] = prefix + "-B" + frame.Batch.astype(int).astype(str)
    # rows named by accession set, unmapped peptides leave
    table = table[~table.index.str.endswith("@ Library Peptides")]
    table.index = table.index.str.findall(r"\|([^|]+)\|").str.join(";")
    assert table.index.is_unique
    # zero means below the detection limit
    values = table[frame.column]
    parts.append(values.where(values != 0))
    frames.append(frame)
samples = pd.concat(frames, ignore_index=True)
values = pd.concat(parts, axis=1)
samples["sample_id"] = samples.column
samples["age_years"] = samples.Age
samples["sex"] = samples.Sex.replace(utils.SEX)
samples["pmi_h"] = samples["PMI (hrs)"]
samples["diagnosis"] = samples.Condition
samples["arm"] = np.where(samples.Condition == "HCF/Low ADNC", "aging", "disease")
# no donor id is deposited; one donor has the same clinical record in every region it was sampled in
record = "Age,Sex,PMI (hrs),APOE Alleles,Study Name,Braak Stage,CERAD Score,Condition,Mutation Status".split(",")
number = samples[record].groupby(record, dropna=False).ngroup()
samples["donor_id"] = "donor" + (number + 1).astype(str).str.zfill(2)
assert len(samples) == 225 and samples.Condition.notna().all()
assert not samples.duplicated(["donor_id", "tissue"]).any() and samples.donor_id.nunique() == 67
assert list(values.columns) == list(samples.column) and (samples.arm == "aging").sum() == 37
# UniProt gene of each accession in the group
accessions = values.index.to_series()
symbols = gene_of.str.split(";").str[0].to_dict()
genes = accessions.str.split(";").map(lambda names: ";".join(dict.fromkeys(symbols[a] for a in names if a in symbols)))
save_groups(dataset, values, samples, accessions, genes, missing=None)

# Blumenreich 2024 five brain regions: MaxQuant LFQ and supplementary table (occipital cortex)
# quantities: MaxQuant LFQ of the deposited proteinGroups; occipital cortex from Supplementary Data 3
dataset = "blumenreich2024_pxd047134"
SEARCHES = {
    "STR": ("9619_STR", "Striatum"),
    "SN": ("10401_SN", "Substantia nigra"),
    "CG": ("10017_CG", "Cingulate gyrus"),
    "MTG": ("9105_MTG", "Middle temporal gyrus"),
}
FLAGS = ["Reverse", "Potential contaminant", "Only identified by site"]
COLUMNS = ["Majority protein IDs", "Gene names", *FLAGS]
PADDING = r"(?<!\d)0+(?=\d)"  # case numbers are zero-padded differently between searches
OCC_FLOOR = 17.339  # value imputed for missing intensities of the occipital cortex
# contaminant or site-only groups of the deposited occipital search
OCC_FLAGGED = "ALB CLDN10 CLEC3B GSN HBE1 KRT1 KRT10 KRT14 KRT2 KRT222 KRT9 LUM TMSB10 TUBB1".split()
# gene names that Excel turned into month-year dates
EXCEL_GENES = {"2002-03": "MARC2", "2005-03": "MARCH5", "2015-09": "SEP15"}
EXCEL_GENES |= {f"{2000 + n}-09": f"SEPT{n}" for n in range(2, 12)}

# donors by author id, case numbers only for matching
donors = pd.read_excel(RAW / dataset / "pride/data_tissue_key.xlsx", sheet_name="Patient_data")
donors = donors.rename(
    columns={"ID": "donor_id", "Group": "diagnosis", "Age": "age_years", "Mutations": "gba_mutation"}
)
donors["sex"] = donors.Sex.replace(utils.SEX)
donors["key"] = donors["Sample Name"].str.replace("/", "-").str.replace(PADDING, "", regex=True)
donors = donors[["donor_id", "diagnosis", "sex", "age_years", "gba_mutation", "key"]]
# run set of the raw file used in MaxQuant
runs = pd.read_csv(RAW / dataset / "manifests/raw_file_map.tsv", sep="\t", dtype=str)
runs = runs[(runs.acquisition == "DDA") & (runs.sample_type == "donor")]
runs = runs.sort_values("used_in_maxquant", ascending=False).drop_duplicates(["region", "author_id"])
run_set = runs.set_index(runs.region + "_" + runs.author_id).run_set
# Supplementary Data 3, samples in rows
supp = pd.read_excel(RAW / dataset / "supp/41531_2024_645_MOESM4_ESM.xlsx")

# log2 LFQ, groups named by first protein
parts, notes, frames = [], [], []
for code, (search, tissue) in SEARCHES.items():
    path = RAW / dataset / f"maxquant/{search}/proteinGroups.txt"
    table = pd.read_csv(
        path, sep="\t", usecols=lambda name: name.startswith("LFQ intensity ") or name in COLUMNS, low_memory=False
    )
    table = table[~table[FLAGS].eq("+").any(axis=1)]
    lfq = match(table.columns, r"^LFQ intensity (?:Control|PD-GBA|PD)-(?P<name>.+)$")
    lfq["key"] = lfq.name.str.replace("_", "-").str.replace(PADDING, "", regex=True)
    samples = donors.merge(lfq, on="key").assign(tissue=tissue, search=search, source="maxquant")
    samples["sample_id"] = code + "_" + samples.donor_id
    assert len(samples) == len(lfq)
    # zero intensity means not detected
    values = table[samples.column].astype(float)
    values = np.log2(values.where(values > 0)).set_axis(samples.sample_id, axis=1)
    lead = table["Majority protein IDs"].str.split(";").str[0]
    seen = values.notna().any(axis=1).to_numpy()
    parts.append(values[seen].set_axis(lead[seen]))
    notes.append(table.loc[seen, ["Majority protein IDs", "Gene names"]].set_axis(lead[seen]))
    frames.append(samples)

# occipital cortex, one column per gene
occ = supp[supp["Brain Region"] == "Occipital cortex"].set_index("Sample ID").drop(columns=["Brain Region", "Group"])
occ = occ.dropna(axis=1, how="all")
occ.columns = [name if isinstance(name, str) else EXCEL_GENES[name.strftime("%Y-%m")] for name in occ.columns]
occ = occ.T[donors.donor_id]
assert occ.min().min().round(3) == OCC_FLOOR
occ = occ.mask(occ.round(3) == OCC_FLOOR)
# gene columns with identical values are one protein group
group = occ.groupby(list(occ.columns), dropna=False, sort=False).ngroup().to_numpy()
genes = pd.Series(occ.index).groupby(group).agg(";".join)
occ = occ.groupby(group).first().set_axis(genes.str.split(";").str[0])
keep = ~occ.index.isin(OCC_FLAGGED)
assert (~keep).sum() == len(OCC_FLAGGED)
samples = donors.assign(sample_id="OCC_" + donors.donor_id, tissue="Occipital cortex", search="8018_OCC", source="supp")
parts.append(occ[keep].set_axis(samples.sample_id, axis=1))
notes.append(pd.DataFrame({"Gene names": genes.to_numpy()[keep]}, index=occ.index[keep]))
frames.append(samples)

# donors of all regions in one table
values = pd.concat(parts, axis=1)
samples = pd.concat(frames, ignore_index=True).drop(columns=["key", "column", "name"])
samples["batch"] = samples.sample_id.map(run_set)
samples["arm"] = np.where(samples.diagnosis == "Control", "aging", "disease")
samples["pmi_h"] = np.nan
# authors excluded donors absent from their table
listed = supp["Brain Region"] + "|" + supp["Sample ID"]
samples["author_excluded"] = ~(samples.tissue + "|" + samples.donor_id).isin(listed)
notes = pd.concat(notes).groupby(level=0, sort=False).first().loc[values.index]
assert values.columns.equals(pd.Index(samples.sample_id))
assert len(samples) == 299 and samples.batch.notna().all() and samples.author_excluded.sum() == 1
lead = notes["Majority protein IDs"].str.split(";").str[0]
gene = notes["Gene names"].str.split(";").str[0]
names = (lead + "_" + gene).fillna(gene).fillna(lead)
save_protein(
    dataset,
    values,
    samples,
    pd.DataFrame({"lead": lead, "gene": gene}).fillna(""),
    notes["Gene names"].fillna(""),
    names,
)

# Ping 2020 dorsolateral prefrontal cortex: TMT11 proteome of control, AsymAD and AD donors
# quantities: raw PD master-protein abundances of the authors' search (PRIDE Proteins.txt)
# PD files F4, F5, F6 are TMT batches 1, 2, 3; channels 126 and 131C hold the pooled standard
dataset = "ping2020_pxd020296"
BATCH_OF_FILE = {"4": "b1", "5": "b2", "6": "b3"}
INFO_COLUMNS = ("Master", "Accession", "Description", "Gene Symbol")
# master proteins with raw PD abundances
table = pd.read_csv(
    RAW / dataset / "pride/totalonly_multiconsensus__Proteins.txt",
    sep="\t",
    usecols=lambda name: name in INFO_COLUMNS or name.startswith("Abundance: "),
)
table = table[table.Master == "IsMasterProtein"].reset_index(drop=True)
# traits of Table 2 (PMI: death to discovery of the body)
traits = pd.read_excel(RAW / dataset / "pride/sampletraits.xlsx", header=3, nrows=33)
traits.columns = traits.columns.str.strip()
traits = traits[["Batch", "Channel", "Diagnosis", "CERAD", "Braak", "PMI", "Age", "Sex"]]
traits = traits.rename(
    columns={
        "Batch": "batch",
        "Channel": "channel",
        "Diagnosis": "diagnosis",
        "CERAD": "cerad",
        "Braak": "braak",
        "PMI": "pmi_h",
        "Age": "age_years",
        "Sex": "sex",
    }
)
traits["batch"] = "b" + traits.batch.astype(int).astype(str)
traits["channel"] = traits.channel.astype(str)
# donor channels, the pooled standard is no sample
samples = match(table.columns, r"^Abundance: F(?P<pd_file>\d): (?P<channel>\w+), (?:Control|Sample)$")
samples["batch"] = samples.pd_file.map(BATCH_OF_FILE)
samples = samples.merge(traits, on=["batch", "channel"], how="left")
samples = samples[samples.diagnosis != "GIS"].reset_index(drop=True)
samples["sample_id"] = samples.batch + "_" + samples.channel
samples["donor_id"] = samples.sample_id
samples["tissue"] = "DLPFC"
samples["arm"] = np.where(samples.diagnosis == "control", "aging", "disease")
assert len(samples) == 27 and samples.age_years.notna().all()
assert samples.sex.isin(["male", "female"]).all() and samples.batch.value_counts().eq(9).all()
assert samples.diagnosis.value_counts().to_dict() == {"control": 10, "AD": 9, "AsymAD": 8}
assert not (table[samples.column] == 0).any().any()
# gene from the FASTA header, else the first annotated symbol
genes = table.Description.str.extract(r"\bGN=(\S+)")[0].fillna(table["Gene Symbol"].str.split(";").str[0].str.strip())
save_groups(dataset, np.log2(table[samples.column].astype(float)), samples, table.Accession, genes, missing="NA")


# ==== External validation proteomes ====
# Mouse brain ageing proteomes of the external validation set (08 reads them as the context mouse_protein_external)
# control animals with known ages only, the arms of disease models and interventions are not prepared


def gene_symbols(accessions):
    # UniProt gene of each accession, empty when unknown
    return accessions.str.split("-").str[0].map(gene_of).fillna("")


def maxquant_groups(table):
    # protein groups without reverse, contaminant and site-only hits
    # decoy and contaminant accessions leave the majority ids of the others
    table = table[~table[["Reverse", "Potential contaminant", "Only identified by site"]].eq("+").any(axis=1)]
    ids = table["Majority protein IDs"].str.replace(r"(?:^|;)(?:REV|CON)__[^;]*", "", regex=True)
    return table.assign(**{"Majority protein IDs": ids.str.lstrip(";")})


# Keele 2023 striatum, cerebellum, hippocampus (TMT10 plexes with a bridge channel), Table S1
# quantities: authors' log2 ratios to the bridge of the plex, mice and tissues in Table S1 as published
# cells without reporter signal in the mouse's channel or in the bridge are missing
dataset = "keele2023_pxd034029"
book = RAW / dataset / "figshare/supp_tables/table_s1.xlsx"
mice = pd.read_excel(book, sheet_name="Mouse annotations", dtype=str).set_index("mouse.id")
annotation = pd.read_excel(book, sheet_name="Protein annotations", dtype=str).set_index("protein.id")
# UniProt accession of an Ensembl protein: the reviewed entry listing it, else the only reviewed entry of its gene
mapping = pd.read_csv(
    RAW / dataset / "uniprot/MOUSE_10090_idmapping_selected.tab.gz",
    sep="\t",
    header=None,
    usecols=[0, 18, 20],
    names=["accession", "gene", "protein"],
    dtype=str,
    keep_default_na=False,
)
mapping["reviewed"] = mapping.accession.isin(reviewed.index)


def listing(column):
    # accession per Ensembl id of a column, reviewed entries first
    rows = mapping.assign(id=mapping[column].str.split("; ")).explode("id")
    rows["id"] = rows.id.str.replace(r"\.\d+$", "", regex=True)
    return rows.sort_values(["reviewed", "accession"], ascending=[False, True])


by_protein = listing("protein").groupby("id").agg(accession=("accession", "first"), reviewed=("reviewed", "first"))
by_gene = listing("gene").query("reviewed").groupby("id").accession.agg(list)
by_gene = by_gene[by_gene.str.len() == 1].str[0]
ids = pd.DataFrame({"protein": annotation.index, "gene": annotation["gene.id"].to_numpy()}, index=annotation.index)
direct = ids.protein.map(by_protein.accession)
lead = direct.where(ids.protein.map(by_protein.reviewed).fillna(False).astype(bool), ids.gene.map(by_gene))
lead = lead.fillna(direct).fillna(ids.protein)
TISSUES = {"Striatum": "Str", "Cerebellum": "Cer", "Hippocampus": "Hip"}
PEPTIDES = RAW / dataset / "figshare/data_zip/data/PeptideQuantitation"
parts, frames = [], []
for tissue, short in TISSUES.items():
    table = pd.read_excel(book, sheet_name=f"{tissue} proteins", dtype={"mouse.id": str}).set_index("mouse.id")
    # summed reporter signal of each protein in each channel; zero signal is not a measurement
    channels = [f"rq_{tag}_sn" for tag in [*mice.tag.unique(), "131c"]]
    psm = pd.read_csv(next(PEPTIDES.glob(f"*_{tissue}.tsv")), sep="\t", usecols=["ProteinId", "Class", *channels])
    psm = psm[psm.ProteinId.str.contains("ENSMUSP") & ~psm.ProteinId.str.contains("##|contaminant")]
    psm["ProteinId"] = psm.ProteinId.str.replace(r"\.\d+$", "", regex=True)
    signal = psm.groupby(["Class", "ProteinId"])[channels].sum()
    frame = mice.loc[table.index].reset_index()
    empty = [
        (signal.loc[plex, f"rq_{tag}_sn"].reindex(table.columns) == 0).to_numpy()
        | (signal.loc[plex, "rq_131c_sn"].reindex(table.columns) == 0).to_numpy()
        for plex, tag in zip(frame.Batch, frame.tag, strict=True)
    ]
    parts.append(table.mask(np.vstack(empty)).T.set_axis(short + "_" + table.index, axis=1))
    frame = frame.assign(tissue=tissue, sample_id=short + "_" + frame["mouse.id"], donor_id=frame["mouse.id"])
    frame["age_months"] = frame.Age.map({"Y": 8, "O": 18})
    frame["sex"] = frame.Sex.replace(utils.SEX)
    frame["batch"] = tissue.lower() + "_" + frame.Batch
    frames.append(frame)
samples = pd.concat(frames, ignore_index=True).rename(columns={"Batch": "plex", "tag": "tmt_tag"})
samples["column"] = samples.sample_id
samples = samples[["column", "sample_id", "donor_id", "tissue", "sex", "age_months", "batch", "plex", "tmt_tag"]]
values = pd.concat(parts, axis=1)[samples.sample_id]
assert samples.groupby("tissue").size().to_dict() == {"Cerebellum": 20, "Hippocampus": 20, "Striatum": 18}
save_groups(dataset, values, samples, lead.reindex(values.index), annotation.symbol.reindex(values.index), missing=None)

# Scifo 2025 whole brain (Proteome Discoverer 3.0, label-free), mzTab of the jPOST deposit
# quantities: PD protein abundances per run, not normalised; PD contaminant entries leave
dataset = "scifo2025_pxd058684"
lines = gzip.open(RAW / dataset / "jpost/Mar032021_AgeingMouse_Brain_LFQ_120min_60K_TopS.MTD-PRT-PEP.mzTab.gz", "rt")
lines = lines.read().replace("\r", "").split("\n")
table = pd.read_csv(
    io.StringIO("\n".join(line for line in lines if line.startswith(("PRH", "PRT")))), sep="\t", dtype=str
)
found = (
    re.match(r"MTD\tms_run\[(\d+)\]-location\t.*\\(\w{9})_(?:TL_)?Br_(\d+)_(\d+)M_(c\d+)_", line) for line in lines
)
samples = pd.DataFrame([m.groups() for m in found if m], columns=["assay", "batch", "run_no", "age_months", "donor_id"])
samples["column"] = "protein_abundance_assay[" + samples.assay + "]"
samples = samples.assign(tissue="Whole brain", sex="male", sample_id="Br_" + samples.donor_id)
samples["age_months"] = samples.age_months.astype(int)
samples = samples[["column", "sample_id", "donor_id", "tissue", "sex", "age_months", "batch", "run_no"]]
assert len(samples) == 45 and samples.groupby("age_months").size().to_dict() == {3: 7, 5: 9, 8: 8, 14: 9, 20: 7, 26: 5}
# contaminant entries are not Swiss-Prot mouse entries
table = table[table.description.str.endswith("[OS=Mus musculus]")]
values = np.log2(table[samples.column].replace("null", np.nan).astype(float).replace(0, np.nan))
keep = values.notna().any(axis=1).to_numpy()
save_groups(dataset, values[keep], samples, table.accession[keep], gene_symbols(table.accession[keep]), missing=None)

# Yu 2020 whole brain (TMT10, one plex), male mice at 16 and 80 weeks
# quantities: Dataset S5 relative abundances (rows sum to 100), zeros missing
dataset = "yu2020_pxd017385"
table = pd.read_excel(RAW / dataset / "supp/pnas.1919410117.sd05.xlsx", sheet_name="Brn", header=1)
samples = match(table.columns, r"^(?P<group>Young|Old) (?P<replicate>\d)$")
samples = samples.assign(tissue="Brain", sex="male", batch="tmt10_plex1")
samples["sample_id"] = "Brn_" + samples.group + samples.replicate
samples["donor_id"] = samples.group.str[0] + samples.replicate
samples["age_months"] = samples.group.map({"Young": 16, "Old": 80}).mul(7 / 30.4375).round(2)
samples = samples[["column", "sample_id", "donor_id", "tissue", "sex", "age_months", "batch", "group"]]
values = np.log2(table[samples.column].astype(float).replace(0, np.nan))
accessions = table["Protein.Id"].str.split("|").str[1]
save_groups(dataset, values, samples, accessions, table["Gene.Symbol"], missing=None)

# Huang 2024 hippocampus (TMTpro, one plex), male mice at 3 and 20 months
# quantities: MaxQuant corrected reporter intensities (reporter ion MS2), groups without signal leave
dataset = "huang2024_pxd043352"
table = maxquant_groups(pd.read_csv(RAW / dataset / "pride/proteinGroups.txt", sep="\t", low_memory=False))
# reporters 1-10 are the TMTpro channels in mass order; the sample of each channel is in the PRIDE protocol
protocol = json.loads((RAW / dataset / "pride/PXD043352.project.json").read_text())["sampleProcessingProtocol"]
design = pd.DataFrame(
    re.findall(r"(Young|Aged), replicate (\d) (\d{3}[NC]?)", protocol), columns=["group", "replicate", "channel"]
)
channels = pd.DataFrame({"channel": "129C 130N 130C 131N 131C 132N 132C 133N 133C 134N".split()})
channels["column"] = "Reporter intensity corrected " + (channels.index + 1).astype(str)
samples = channels.merge(design).assign(tissue="Hippocampus", sex="male", batch="tmtpro_plex1")
samples["sample_id"] = "Hip_" + samples.group + samples.replicate
samples["donor_id"] = samples.group.str[0] + samples.replicate
samples["age_months"] = samples.group.map({"Young": 3, "Aged": 20})
assert len(samples) == 10
samples = samples[["column", "sample_id", "donor_id", "tissue", "sex", "age_months", "batch", "channel"]]
values = table[samples.column].astype(float)
values = np.log2(values.where(values > 0))
seen = values.notna().any(axis=1).to_numpy()
save_groups(
    dataset,
    values[seen],
    samples,
    table["Majority protein IDs"][seen],
    table["Gene names"][seen],
    lead="reviewed_first_gene",
    missing=None,
)

# Jha 2026 cortex (TMTpro 16-plex), male mice: age controls at 3-4, 12-14 and 24 months (3 each)
# the sleep-restricted and recovery groups of young mice are experimental arms and leave
# quantities: MaxQuant corrected reporter intensities, zeros missing
dataset = "jha2026_pxd060735"
table = maxquant_groups(pd.read_csv(RAW / dataset / "pride/proteinGroups.txt", sep="\t", low_memory=False))
channels = [c for c in table.columns if c.startswith("Reporter intensity corrected ")][:15]
# Supplementary Table 2 lists the 15 samples by group; its columns are matched to the reporter channels by value
labels = pd.read_excel(RAW / dataset / "supp/MOESM3_supplementary_table_2.xlsx").set_index("Protein")
labels = labels[~labels.index.duplicated(keep=False)]
first = table["Gene names"].fillna("").str.split(";").str[0]
mine = table[channels].astype(float).where(lambda x: x > 0).set_axis(first)
mine = mine[~mine.index.duplicated(keep=False) & mine.index.isin(labels.index)]
gap = pd.DataFrame(
    {
        label: np.log2(mine).sub(np.log2(labels.loc[mine.index, label].where(lambda x: x > 0)), axis=0).abs().median()
        for label in labels
    }
)
assert (gap.min() < 1e-3).all() and list(gap.idxmin()) == channels
samples = pd.DataFrame(
    {"column": channels, "group": [c[:-1] for c in labels.columns], "replicate": [c[-1] for c in labels.columns]}
)
samples = samples[samples.group.isin(["NS", "MID", "OLD"])]
samples = samples.assign(
    tissue="Cortex", sex="male", batch="tmtpro_plex1", sample_id="Ctx_" + samples.group + samples.replicate
)
samples["donor_id"] = samples.group + samples.replicate
samples["age_months"] = samples.group.map({"NS": 3.5, "MID": 13, "OLD": 24})
samples = samples[["column", "sample_id", "donor_id", "tissue", "sex", "age_months", "batch", "group"]]
values = table[samples.column].astype(float)
values = np.log2(values.where(values > 0))
seen = values.notna().any(axis=1).to_numpy()
save_groups(
    dataset,
    values[seen],
    samples,
    table["Majority protein IDs"][seen],
    table["Gene names"][seen],
    lead="reviewed_first_gene",
    missing=None,
)

# Oliviero 2022 whole-brain total lysate (label-free, Progenesis), male mice at 3, 5, 10 and 15 months (3 each)
# the chromatin fraction is a partial proteome and is not used
# quantities: Supplementary file 2 run abundances, zeros missing; the two injections of a mouse are averaged
dataset = "oliviero2022_elife"
book = RAW / dataset / "elife/elife-73524-supp2-v2.xlsx"
header = pd.read_excel(book, sheet_name="Sheet2", header=None, nrows=3)
table = pd.read_excel(book, sheet_name="Sheet2", header=2)
runs = pd.DataFrame({"column": header.iloc[2], "mouse": header.iloc[1].ffill()})
runs = runs[runs.column.astype(str).str.startswith("QEHF2_")]
runs["mouse"] = runs.mouse.map("{:g}".format)
values = table[runs.column].astype(float).replace(0, np.nan)
# one injection column of the file is a copy of another; the copy (labelled mouse 3.2) leaves
copy = values.T.duplicated().to_numpy()
assert runs.column[copy].tolist() == ["QEHF2_03823_GO"]
runs, values = runs[~copy], values.loc[:, ~copy]
samples = runs.drop_duplicates("mouse").assign(tissue="Whole brain", sex="male", batch="progenesis_lysate")
samples["age_months"] = samples.mouse.str.split(".").str[0].astype(int)
samples["donor_id"] = samples.mouse
samples["sample_id"] = "Lys_" + samples.age_months.astype(str) + "mo_" + samples.mouse.str.split(".").str[1]
samples = samples[["mouse", "sample_id", "donor_id", "tissue", "sex", "age_months", "batch"]].rename(
    columns={"mouse": "column"}
)
mean = pd.DataFrame({m: values[g.column].mean(axis=1) for m, g in runs.groupby("mouse")})
values = np.log2(mean[samples.column])
genes = table.Description.str.extract(r"\bGN=(\S+)")[0]
seen = values.notna().any(axis=1).to_numpy()
save_groups(dataset, values[seen], samples, table.Accession[seen], genes[seen], missing=None)

# Mohallem 2024 whole brain (label-free MaxQuant, four protease digests combined), male mice by age group
# exact ages are not given: adults 3-4, mid-age 10 and old 19-21 months, group midpoints used
# quantities: MaxQuant LFQ intensities of Supplemental Table S1, zeros missing
dataset = "mohallem2024_msv000094441"
table = maxquant_groups(pd.read_excel(RAW / dataset / "supp/mmc2.xlsx", sheet_name="STable1", header=1))
samples = match(table.columns, r"^LFQ intensity (?P<group>Adult|MidAge|Old)-(?P<replicate>\d)$")
samples = samples.assign(tissue="Whole brain", sex="male", batch="maxquant_lfq")
samples["sample_id"] = samples.group + "-" + samples.replicate
samples["donor_id"] = samples.sample_id
samples["age_months"] = samples.group.map({"Adult": 3.5, "MidAge": 10, "Old": 20})
samples = samples[["column", "sample_id", "donor_id", "tissue", "sex", "age_months", "batch", "group"]]
assert samples.group.value_counts().to_dict() == {"Adult": 5, "Old": 5, "MidAge": 3}
values = np.log2(table[samples.column].astype(float).replace(0, np.nan))
seen = values.notna().any(axis=1).to_numpy()
save_groups(
    dataset,
    values[seen],
    samples,
    table["Majority protein IDs"][seen],
    table["Gene names"][seen],
    lead="reviewed_first_gene",
    missing=None,
)

# Marchionini 2025 striatum of wild-type LacQ140-background mice (DIA, DIA-NN 1.7.12), 2, 6, 9 and 12 months, 10 each
# all mice carry the Lac repressor transgene and lifelong IPTG; the Htt LacQ140 and untreated groups leave
# age is an acquisition block (2 and 9 months series X2, 6 and 12 months series X1); sex is not deposited
# quantities: MaxLFQ (algorithm of the diann R package) of the precursors in the DIA-NN report, which is
# streamed from PRIDE and deposited with the data release; Precursor.Normalised at 1% precursor and group FDR
dataset = "marchionini2025_pxd054495"
MARGIN = -10.0  # log quantities below this are missing (diann_maxlfq default)


def maxlfq(x):
    # MaxLFQ of one protein group, x = precursors x runs in natural log, NaN missing
    # median log ratios of all run pairs solved by least squares, with a ridge towards the best precursor
    x = np.where(x > MARGIN, x, np.nan)
    best = np.nanmax(np.where(np.isnan(x), -np.inf, x), axis=0)
    use = best > MARGIN
    result = np.where(use, best, np.nan)
    if len(x) == 1 or use.sum() < 2:
        return result
    x = x[:, use]
    gap = np.sort(x[:, :, None] - x[:, None, :], axis=0)  # precursors x run i x run j, NaN last
    count = (~np.isnan(gap)).sum(axis=0)
    middle = [np.take_along_axis(gap, np.maximum(k, 0)[None], axis=0)[0] for k in ((count - 1) // 2, count // 2)]
    pair = count > 0
    np.fill_diagonal(pair, False)
    partners = pair.sum(axis=1)
    ridge = 1e-4 * np.maximum(partners, 1)
    a = -pair.astype(float) + np.diag(partners + ridge)
    b = np.where(pair, (middle[0] + middle[1]) / 2, 0).sum(axis=1) + best[use] * ridge
    result[use] = np.linalg.solve(a, b)
    return np.where(result > MARGIN, result, np.nan)


runs = pd.read_excel(RAW / dataset / "pride/sample_overview.xlsx", dtype=str)
runs = runs[runs.Condition.str.match(r"^M\d+_WT_IPTG$")].assign(Run=lambda t: t["Raw file name"] + ".raw")
precursors = pd.read_parquet(
    RAW / dataset / "pride_stream/precursors.parquet",
    columns=["Run", "Protein.Group", "Precursor.Id", "Precursor.Normalised", "Q.Value", "PG.Q.Value"],
    filters=[("Run", "in", list(runs.Run))],
)
precursors = precursors[
    (precursors["Q.Value"] <= 0.01) & (precursors["PG.Q.Value"] <= 0.01) & (precursors["Precursor.Normalised"] > 0)
]
precursors = precursors.assign(
    run=pd.Index(runs.Run).get_indexer(precursors.Run),
    precursor=pd.factorize(precursors["Precursor.Id"])[0],
    log=np.log(precursors["Precursor.Normalised"].astype(float)),
).sort_values("Protein.Group", kind="stable")
group = precursors["Protein.Group"].to_numpy()
bounds = np.flatnonzero(np.r_[True, group[1:] != group[:-1], True])
quantities = {}
for start, stop in zip(bounds[:-1], bounds[1:], strict=True):
    rows = precursors.iloc[start:stop]
    _, row = np.unique(rows.precursor, return_inverse=True)
    x = np.full((row.max() + 1, len(runs)), np.nan)
    x[row, rows.run.to_numpy()] = rows["log"].to_numpy()
    quantities[group[start]] = maxlfq(x)
values = pd.DataFrame(quantities, index=runs.Run).T / np.log(2)
samples = runs.rename(columns={"Run": "column"}).assign(tissue="Striatum", genotype="WT")
samples["donor_id"] = "CHDI_" + samples.Sample_No_CHDI
samples["age_months"] = samples["Age group"].str.removesuffix("M").astype(int)
samples["batch"] = samples.column.str.extract(r"_(X\d)_")[0]
samples["run_number"] = samples.column.str.extract(r"^I(\d+)_")[0].astype(int)
samples["sample_id"] = (
    "Str_" + samples.age_months.astype(str) + "m_" + samples.Sample_No_CHDI.str.split("_").str[0].str.zfill(2)
)
# sex from the Y-chromosome protein Ddx3y (Q62095), detected in the runs of males
samples["sex"] = np.where(values.loc["Q62095", samples.column].notna().to_numpy(), "male", "female")
samples = samples[["column", "sample_id", "donor_id", "tissue", "sex", "age_months", "batch", "genotype", "run_number"]]
assert len(samples) == 40 and samples.groupby("age_months").size().eq(10).all()
members = values.index.to_series().str.split(";")
genes = members.map(
    lambda accessions: ";".join(dict.fromkeys(g for g in map(gene_of.get, (a.split("-")[0] for a in accessions)) if g))
)
seen = values.notna().any(axis=1).to_numpy()
save_groups(
    dataset,
    values[seen],
    samples,
    values.index.to_series()[seen],
    genes[seen],
    lead="reviewed_first_gene",
    missing=None,
)

# Vignane 2026 frontal cortex (label-free DDA, PEAKS), wild-type male mice at 10 weeks, 10 and 18 months (4 each)
# quantities: log2 TIC-normalised intensities published in PRIDE, the 18-month samples from the 10-week comparison sheet
dataset = "filipovic2026_pxd045785"
book = RAW / dataset / "pride/Results_Totalproteome.xlsx"
young = pd.read_excel(book, sheet_name="WT 10 weeks vs 10 months", na_values=["NaN"]).set_index("Accession", drop=False)
old = pd.read_excel(book, sheet_name="WT 10 weeks vs 18 months", na_values=["NaN"]).set_index("Accession")
table = young.join(old[[c for c in old.columns if c.startswith("WT 18 months")]])
samples = match(table.columns, r"^WT (?P<age>10 weeks|10 months|18 months) (?P<replicate>\d)$")
samples = samples.assign(tissue="Frontal cortex", sex="male", batch="timstof_series1")
samples["sample_id"] = "FC_" + samples.age.str.replace(" ", "_") + "_" + samples.replicate
samples["donor_id"] = samples.sample_id.str[3:]
samples["age_months"] = samples.age.map({"10 weeks": round(10 * 7 / 30.4375, 1), "10 months": 10, "18 months": 18})
samples = samples[["column", "sample_id", "donor_id", "tissue", "sex", "age_months", "batch"]]
assert len(samples) == 12
genes = table.Description.str.extract(r"\bGN=(\S+)")[0]
save_groups(dataset, table[samples.column].astype(float), samples, table.Accession, genes, missing=None)

# Gebru 2024 hippocampus and amygdala (diaPASEF), wild-type mice of both sexes at 3 and 18 months
# the Fkbp5 knockouts leave; the 18-month mice went through a 30-day circadian protocol
# quantities: DIA-NN protein-group MaxLFQ of one search per region (pg_matrix filter), derived from the PRIDE reports
dataset = "gebru2024_pxd049439"
pattern = (
    r"^LB_NG_\w+?_(?P<sex>[MF])(?P<age_months>\d+)(?P<genotype>WT|KO)_(?P<replicate>\d+)_DIA_.*_(?P<run_number>\d+)$"
)
parts, frames, annotation = [], [], []
for tissue, short in (("Hippocampus", "Hip"), ("Amygdala", "Amy")):
    table = pd.read_csv(
        RAW / dataset / f"derived/{tissue}_pg_maxlfq.tsv.gz", sep="\t", keep_default_na=False, na_values=[""]
    )
    frame = match(table.columns, pattern).query("genotype == 'WT'")
    frame = frame.assign(tissue=tissue, batch=f"diann_{tissue.lower()}_search")
    frame["sex"] = frame.sex.replace({"M": "male", "F": "female"})
    frame["age_months"] = frame.age_months.astype(int)
    frame["donor_id"] = frame.column.str.extract(r"_([MF]\d+WT_\d+)_DIA")[0]
    frame["sample_id"] = short + "_" + frame.donor_id
    parts.append(np.log2(table.set_index("Protein.Group")[frame.column].astype(float)))
    annotation.append(table.set_index("Protein.Group")[["Genes"]])
    frames.append(frame)
samples = pd.concat(frames, ignore_index=True)
samples = samples[["column", "sample_id", "donor_id", "tissue", "sex", "age_months", "batch", "genotype", "run_number"]]
values = pd.concat(parts, axis=1)
annotation = pd.concat(annotation).groupby(level=0).first().loc[values.index]
assert samples.groupby("tissue").size().to_dict() == {"Amygdala": 22, "Hippocampus": 21}
save_groups(
    dataset,
    values,
    samples,
    values.index.to_series(),
    annotation.Genes.fillna(""),
    lead="reviewed_first_gene",
    missing=None,
)

# Yarbro 2025 cortex (TMT, nine batches), wild-type littermates of the amyloidosis models, both sexes
# the 5xFAD and App knock-in mice leave; batch is largely confounded with age
# quantities: authors' normalised log2-scale reporter intensities of Supplementary Data sheet 3 (JUMP), NA missing
dataset = "yarbro2025_nat_commun"
table = pd.read_excel(
    RAW / dataset / "supp/MOESM3_supplementary_data_1-24.xlsx", sheet_name="3", header=None, na_filter=False
)
names = table.iloc[4, 8:]
rows = table.iloc[5:, 1:].reset_index(drop=True)
rows.columns = ["accession", "gene", "description", "psms", "peptides", "unique_peptides", "coverage", *names]
samples = match(rows.columns, r"^(?P<batch>B\d)_WT_(?P<sex>[MF])_(?P<age_months>\d+)m\d*$")
samples = samples.assign(tissue="Cortex", sample_id=samples.column, donor_id=samples.column)
samples["sex"] = samples.sex.replace(utils.SEX)
samples["age_months"] = samples.age_months.astype(int)
samples = samples.assign(genotype="WT")[
    ["column", "sample_id", "donor_id", "tissue", "sex", "age_months", "batch", "genotype"]
]
assert samples.groupby("age_months").size().to_dict() == {3: 6, 6: 13, 12: 14, 18: 6}
# the custom human Abeta42 entry of the transgenes is no mouse protein
rows = rows[~rows.accession.str.startswith("cu|")]
values = np.log2(rows[samples.column].replace("NA", np.nan).astype(float))
accessions = rows.accession.str.split("|").str[1]
genes = rows.gene.astype(str).str.strip().replace("NA", "")
save_groups(dataset, values, samples, accessions, genes, missing=None)

# Chen 2025 prefrontal cortex of AD-BXD F1 mice (DIA, 35 batches), non-transgenic littermates at 6 and 14 months
# the 5xFAD carriers leave; a label is dropped when the human APP signal of the run contradicts it
# quantities: Level 3B batch-adjusted log2 peptide-area sums per protein group, as published
dataset = "chen2025_pxd045403"
meta = pd.read_excel(RAW / dataset / "panorama/metadata/MouAD-PFC-meta-SAGE-2.0.xlsx", sheet_name="Metadata", dtype=str)
meta = meta[meta.Genotype.notna()].set_index("Replicate", drop=False)
table = pd.read_csv(RAW / dataset / "panorama/level3b/MouAD_PFC_AXBDX_Protein-BatchAdj-v40125.csv")
table = table.drop(columns="Peptide").set_index("Protein")
# human APP (P05067, transgene) against mouse App: carriers lie above zero in the log2 difference
transgene = table.loc["sp|P05067|A4", meta.index] - table.loc["sp|P12023|A4", meta.index]
samples = meta[(meta.Genotype == "WT") & (transgene <= 0)].rename(columns={"Replicate": "column"})
samples = samples.assign(tissue="Prefrontal cortex", sample_id=samples.column, donor_id=samples.column)
samples["age_months"] = samples.Age.map({"Young": 6, "Old": 14})
samples["sex"] = samples.Sex.str.lower()
samples["batch"] = samples.Batch
samples["sire_strain"] = samples.Strain
samples = samples.assign(genotype="Ntg")[
    ["column", "sample_id", "donor_id", "tissue", "sex", "age_months", "batch", "sire_strain", "genotype"]
]
assert samples.groupby("age_months").size().to_dict() == {6: 83, 14: 105}
# library peptides without a protein, the yeast enolase spike-in and the human transgene products leave
table = table.drop(index=[g for g in table.index if re.search(r"Library Peptides|\|P00924\||\|P05067\||\|P49768\|", g)])
members = table.index.to_series().str.split(", ").map(lambda groups: [g.split("|")[1].split("-")[0] for g in groups])
accessions = members.str.join(";")
save_groups(
    dataset,
    table[samples.column],
    samples,
    accessions,
    members.str[0].map(gene_of).fillna(""),
    lead="reviewed_first_gene",
    missing=None,
)
