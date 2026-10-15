# Bulk transcriptome matrices: RNA-seq gene counts, log2 expression where only that is released
# Sample tables of the paired cohorts come from 01
# reads data/raw (00_download.sh); writes output/prepared/<dataset>/rna.parquet, rna_features.tsv, rna_samples.tsv
import glob
import gzip
import re
import subprocess
import tempfile

import numpy as np
import pandas as pd

import utils

RAW = utils.DATA / "raw"
hgnc = pd.read_csv(
    RAW / "hgnc/hgnc_complete_set.txt", sep="\t", dtype=str, usecols=["status", "entrez_id", "ensembl_gene_id"]
)
entrez = (
    hgnc[(hgnc.status == "Approved") & hgnc.entrez_id.notna()]
    .drop_duplicates("entrez_id")
    .set_index("entrez_id")
    .ensembl_gene_id
)


def save_rna(dataset, data, samples, id_type, symbols=None):
    # Ensembl versions dropped; duplicates keep the most expressed
    data = data[samples.column].set_axis(samples.sample_id.to_numpy(), axis=1)
    if id_type == "ensembl":
        data.index = data.index.str.split(".").str[0]
    symbols = np.full(len(data), np.nan, dtype=object) if symbols is None else symbols.to_numpy()
    if data.index.has_duplicates:
        keep = np.sort(pd.Series(data.mean(axis=1).to_numpy()).groupby(data.index.to_numpy()).idxmax().to_numpy())
        data, symbols = data.iloc[keep], symbols[keep]
    features = pd.DataFrame(index=data.index.rename("gene_id"))
    features["ensembl"] = {"ensembl": data.index.to_numpy(), "entrez": data.index.map(entrez)}.get(id_type, np.nan)
    features["symbol"] = data.index.to_numpy() if id_type == "symbol" else symbols
    assert data.columns.is_unique
    utils.write(data, utils.PREPARED / dataset / "rna.parquet")
    utils.write(features, utils.PREPARED / dataset / "rna_features.tsv")
    utils.write(samples, utils.PREPARED / dataset / "rna_samples.tsv")


def series_matrix(path):
    # GEO series matrix: one row per sample
    rows = {}
    for line in gzip.open(path, "rt", errors="replace"):
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!Sample_"):
            rows.setdefault(line.split("\t")[0][8:], []).append(
                [v.strip().strip('"') for v in line.rstrip("\n").split("\t")[1:]]
            )
    table = pd.DataFrame(index=pd.Index(rows.pop("geo_accession")[0], name="gsm"))
    for key, lines in rows.items():
        if key.startswith("characteristics"):
            for values in lines:
                for gsm, value in zip(table.index, values, strict=True):
                    name, _, text = value.partition(":")
                    if value:
                        table.loc[gsm, re.sub(r"\W+", "_", name.strip().lower()).strip("_")] = text.strip()
        else:
            table[key] = [" | ".join(v for v in column if v) for column in zip(*lines, strict=True)]
    return table.reset_index()


def read_rse(file):
    # counts, row ranges and sample table of a SummarizedExperiment
    with tempfile.TemporaryDirectory() as folder:
        subprocess.run(["Rscript", utils.ROOT / "pipeline/export_rse.R", file, folder], check=True)
        options = {"sep": "\t", "keep_default_na": False, "na_values": ["NA"]}
        rows = pd.read_csv(f"{folder}/rows.tsv.gz", index_col="row_id", **options)
        columns = pd.read_csv(f"{folder}/coldata.tsv.gz", **options)
        counts = np.fromfile(f"{folder}/counts.int32", dtype="<i4").reshape((len(columns), len(rows))).T
    return pd.DataFrame(counts, index=rows.index, columns=columns["sample"]), rows, columns


# Paired cohorts: Langfelder counts, Takasugi normalised counts and Luo log2 expression
luo = pd.read_excel(
    RAW / "luo2024_syn32136022/supp/41380_2024_2576_MOESM2_ESM.xlsx",
    sheet_name="Table S3",
    header=2,
    index_col="Ensembl ID",
)
paired = {
    "takasugi2024_pxd040722": (
        pd.read_csv(RAW / "takasugi2024_pxd040722/geo/GSE225576_Brain.txt.gz", sep="\t", index_col="Gene_ID"),
        "symbol",
        None,
    ),
    "langfelder2016_pxd006302": (
        pd.read_csv(
            RAW / "langfelder2016_pxd006302/geo/GSE65774_Striatum_mRNA_counts_processedData.txt.gz",
            sep="\t",
            index_col=0,
        )
        .round()
        .astype("int64"),
        "ensembl",
        None,
    ),
    "luo2024_syn32136022": (luo, "ensembl", luo["Gene symbol"]),
}
for dataset, (table, id_type, symbols) in paired.items():
    save_rna(
        dataset,
        table,
        utils.read(utils.PREPARED / dataset / "rna_samples.tsv", dtype={"sample_id": str}),
        id_type,
        symbols,
    )

# Markunas 2021 nucleus accumbens (Salmon counts, TPM gene ids)
dataset = "markunas2020_gse171936"
stage = f"{RAW / dataset}/stage2/GSE171936_nac_rnaseq_salmon_dge_gene_data_"
table = pd.read_csv(stage + "counts_filters.csv.gz").iloc[:, 1:]
table.index = pd.read_csv(stage + "tpm_filters.csv.gz", index_col=0).index
meta = series_matrix(RAW / dataset / "stage2/GSE171936_series_matrix.txt.gz")
meta["_meta"] = meta.title.str.extract(r"^(\d+_R\d+C\d+)_")[0]
samples = pd.DataFrame({"column": [c for c in table.columns if re.match(r"^\d+_R\d+C\d+$", c)]}).merge(
    meta, left_on="column", right_on="_meta", how="left"
)
samples = samples.assign(tissue="Nucleus accumbens", sample_id=samples.column, donor_id=samples.column)
samples["age_years"] = pd.to_numeric(samples.age_at_death, errors="coerce")
samples["smoking"] = samples.current_smoking_status
samples["arm"] = np.where(samples.smoking == "control", "aging", "disease")
prefix = {
    "sex_reported": "sex_m_male",
    "sex_imputed": "nac_assigned_sex",
    "race": "ancestry_cauc",
    "rin": "rin_rna_integrity",
}
for name, start in prefix.items():
    samples[name] = samples.filter(regex=f"^{start}").iloc[:, 0]
for name in ("sex_reported", "sex_imputed", "race"):
    samples[name] = samples[name].replace("NA", None)
samples["rin"] = pd.to_numeric(samples.rin, errors="coerce")
samples["sex_letter"] = samples.sex_reported.fillna(samples.sex_imputed).replace(utils.SEX)
save_rna(dataset, table.round().astype("int64"), samples, "ensembl")

# Matlik 2024 sorted nuclei (counts per sample, Entrez ids)
dataset = "matlik2024_gse227729"
files = sorted(glob.glob(f"{RAW / dataset}/stage2/GSM*_Gene.Counts.txt.gz"))
parts = {f.split("/")[-1].split("_")[0]: pd.read_csv(f, sep="\t", index_col="IDs", dtype={"IDs": str}) for f in files}
table = pd.DataFrame({gsm: x.iloc[:, -1] for gsm, x in parts.items()})
samples = pd.DataFrame({"column": table.columns}).merge(
    series_matrix(RAW / dataset / "stage2/GSE227729_series_matrix.txt.gz"), left_on="column", right_on="gsm", how="left"
)
samples = samples.assign(sample_id=samples.column, donor_id=samples.donor)
samples["age_years"] = pd.to_numeric(samples.age, errors="coerce")
samples["arm"] = np.where(samples.disease == "CTRL", "aging", "disease")
samples["pmi_h"] = pd.to_numeric(samples.pmi_hrs, errors="coerce")
samples["sex"] = samples.sex.replace(utils.SEX)
samples["cell_type"] = samples.cell_type.replace(
    {
        "Striatonigral medium spiny neuron": "D1 MSN",
        "Striatopallidal medium spiny neuron": "D2 MSN",
        "Interneuron, Other": "Interneuron TAC3+",
        "Interneuron, Fast Spiking": "Interneuron FS",
        "Interneuron, Non-Fast Spiking": "Interneuron SST+",
        "Cholinergic interneuron": "Interneuron cholinergic",
    }
)
save_rna(dataset, table, samples, "entrez", parts[list(parts)[-1]]["Symbols"])

# GTEx v10 brain tissues (RNA-SeQC counts)
dataset = "gtexv10_phs000424"
parts, descriptions = [], []
for path in sorted(glob.glob(f"{RAW / dataset}/gene_reads/*.gct.gz")):
    part = pd.read_csv(path, sep="\t", skiprows=2, index_col="Name")
    descriptions.append(part.pop("Description"))
    parts.append(part)
table = pd.concat(parts, axis=1)
annotations = [
    pd.read_csv(
        RAW / dataset / f"annotations/GTEx_Analysis_v10_Annotations_{n}DS.txt",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    for n in ("SampleAttributes", "SubjectPhenotypes")
]
samples = pd.DataFrame({"column": table.columns})
samples["subject"] = samples.column.str.extract(r"^(GTEX-[^-]+)-.+$")[0]
samples = samples.merge(annotations[0], left_on="column", right_on="SAMPID", how="left").merge(
    annotations[1], left_on="subject", right_on="SUBJID", how="left"
)
samples = samples.assign(sample_id=samples.column, donor_id=samples.subject, tissue=samples.SMTSD)
samples["sex"] = samples.SEX.replace({"1": "male", "2": "female"})
samples["age_bin"] = samples.AGE
samples["age_midpoint"] = pd.to_numeric(
    samples.AGE.replace({f"{a}-{a + 9}": a + 4.5 for a in range(20, 80, 10)}), errors="coerce"
)
assert samples.SEX.notna().all() and samples.AGE.notna().all()
samples["dthhrdy"] = samples.DTHHRDY.replace({"": None})
samples["rin"] = pd.to_numeric(samples.SMRIN, errors="coerce")
samples["ischemic_min"] = pd.to_numeric(samples.SMTSISCH, errors="coerce")
samples["center"] = samples.SMCENTER
samples["batch"] = samples.SMNABTCH
samples["extraction_date"] = samples.SMNABTCHD
save_rna(dataset, table, samples, "ensembl", pd.concat(descriptions).groupby(level=0, sort=False).first())

# BrainSeq II DLPFC and hippocampus (Collado-Torres 2019)
dataset = "colladotorres2019_phs000979"
counts, rows, coldata = read_rse(RAW / dataset / "libd/rse_gene_unfiltered.Rdata")
samples = pd.DataFrame({"column": counts.columns}).merge(coldata, left_on="column", right_on="sample", how="left")
assert counts.shape[1] == 900 and samples.BrNum.notna().all()
samples = samples.assign(
    sample_id=samples["sample"],
    donor_id=samples.BrNum,
    age_years=samples.Age,
    sex=samples.Sex.replace(utils.SEX),
)
samples["tissue"] = samples.Region.replace({"HIPPO": "Hippocampus"})
samples["arm"] = samples.Dx.replace({"Control": "aging", "Schizo": "disease"})
samples["race"] = samples.Race
samples["rin"] = samples.RIN
samples["mitorate"] = samples.mitoRate
samples["totalassignedgene"] = samples.totalAssignedGene
save_rna(dataset, counts, samples, "ensembl", rows.Symbol.replace("", np.nan))

# BrainSeq III caudate (Benjamin 2022)
dataset = "benjamin2022_brainseq3_caudate"
counts, rows, coldata = read_rse(RAW / dataset / "stage2/caudate_brainseq_phase3_hg38_rseGene_merged_n464.rda")
samples = pd.DataFrame({"column": counts.columns}).merge(coldata, left_on="column", right_on="sample", how="left")
assert counts.shape[1] == 464 and samples.Age.notna().all()
samples["tissue"] = "Caudate"
samples = samples.assign(
    sample_id=samples["sample"],
    donor_id=samples.BrNum,
    age_years=samples.Age,
    sex=samples.Sex.replace(utils.SEX),
)
samples["race"] = samples.Race
samples["diagnosis"] = samples.Dx
samples["arm"] = np.where(samples.Dx == "Control", "aging", "disease")
samples["batch"] = samples.FlowCell
samples["rin"] = samples.RIN
samples["pmi_h"] = samples.PMI
samples["ph"] = samples.pH
samples["mitorate"] = samples.mitoRate
samples["rrna_rate"] = samples.rRNA_rate
samples["totalassignedgene"] = samples.totalAssignedGene
samples["overallmaprate"] = samples.overallMapRate
save_rna(dataset, counts, samples, "ensembl", rows.Symbol.replace("", np.nan))

# Hahn 2023 bulk RNA-seq of 15 brain regions (Smart-seq2 counts)
dataset = "hahn2023_gse212336"
table = pd.read_csv(RAW / dataset / "geo/GSE212336_BulkSeq_Counttable.txt.gz", sep="\t", index_col=0)
# GEO samples named by the count columns
meta = series_matrix(RAW / dataset / "geo/GSE212336_series_matrix.txt.gz")
meta["column"] = meta.description.str.extract(r"(BulkSeq_\w+)$")[0]
samples = pd.DataFrame({"column": table.columns}).merge(meta, on="column", how="left")
# region code and mouse id from the column name
samples["region"] = samples.column.str.split("_").str[1]
samples = samples.assign(sample_id=samples.column, mouse_id=samples.mouse_id_number, arm="aging", batch="single")
samples["donor_id"] = "mouse" + samples.mouse_id
samples["age_months"] = samples.age_months.astype(int)
samples["sex"] = samples.sex.replace(utils.SEX)
assert table.shape[1] == len(meta) == 847 and samples.gsm.notna().all()
assert samples.sample_id.is_unique and not samples.duplicated(["donor_id", "tissue"]).any()
assert samples.groupby("region").tissue.nunique().eq(1).all() and samples.tissue.nunique() == 15
assert (samples.mouse_id == samples.column.str.split("_").str[2]).all() and samples.mouse_id.nunique() == 59
assert samples[["sex", "age_months"]].notna().all().all() and (table.dtypes == "int64").all()
columns = "column sample_id donor_id mouse_id gsm tissue sex age_months arm batch".split()
save_rna(dataset, table, samples[columns], "symbol")

# Shavlakadze 2026 bulk RNA-seq of four brain regions (Regeneron total RNA counts)
dataset = "shavlakadze2026_gse314817"
BRAIN = ["cortex", "hippocampus", "hypothalamus", "cerebellum"]
# library number links each count column to its GEO sample
meta = series_matrix(RAW / dataset / "geo/GSE314817_series_matrix.txt.gz")
meta["column"] = meta.description.str.removeprefix("Library name: ")
meta = meta[meta.tissue.isin(BRAIN)]
# brain libraries only, fractional multi-mapper counts rounded
table = pd.read_csv(
    RAW / dataset / "geo/GSE314817_mouse_raw_counts_wo_muscle.txt.gz",
    sep="\t",
    index_col="GeneName",
    usecols=["GeneName", *meta.column],
).round()
samples = pd.DataFrame({"column": table.columns}).merge(meta, on="column", how="left")
samples = samples.assign(sample_id=samples.column, donor_id=samples.column, arm="aging", batch="single")
samples["age_months"] = samples.age.str.removesuffix(" months").astype(int)
assert table.shape[1] == len(meta) == 414 and samples.gsm.notna().all()
assert samples.sample_id.is_unique and samples.tissue.nunique() == 4 and set(samples.sex) == {"female", "male"}
assert (samples.title == samples.sex + "_" + samples.tissue + "_" + samples.age + "_" + samples.column).all()
assert samples.age_months.isin([6, 12, 18, 21, 24, 27]).all() and (table.to_numpy() >= 0).all()
columns = "column sample_id donor_id gsm tissue sex age_months arm batch".split()
save_rna(dataset, table.astype("int64"), samples[columns], "symbol")

# Ramasamy 2014 UKBEC exon arrays of ten brain regions (RMA log2 per transcript cluster)
dataset = "ramasamy2014_gse46706"
SERIES = RAW / dataset / "geo/GSE46706-GPL5175_series_matrix.txt.gz"
# expression table below the header lines of the series matrix
table = pd.read_csv(SERIES, sep="\t", comment="!", index_col="ID_REF", dtype={"ID_REF": str}).astype("float32")
meta = series_matrix(SERIES)
# first listed gene of each core cluster, platform table below 14 header lines
platform = pd.read_csv(
    RAW / dataset / "platform/GPL5175_data.txt.gz",
    sep="\t",
    skiprows=14,
    index_col="ID",
    usecols=["ID", "gene_assignment"],
    dtype=str,
)
first = platform.gene_assignment.reindex(table.index).str.split(" /// ").str[0].str.split(" // ")
genes = pd.DataFrame({"entrez": first.str[4], "symbol": first.str[1]}).dropna()
# sample fields from the GEO characteristics
samples = pd.DataFrame({"column": table.columns}).merge(meta, left_on="column", right_on="gsm", how="left")
samples = samples.assign(sample_id=samples.column, donor_id=samples.individual_id, tissue=samples.brain_region)
samples["sex"] = samples.sex.replace(utils.SEX)
samples["age_years"] = pd.to_numeric(samples.age_at_death_in_years)
samples["diagnosis"] = samples.disease_status
samples["arm"] = np.where(samples.diagnosis == "neurologically and neuropathologically normal", "aging", "disease")
samples["pmi_h"] = pd.to_numeric(samples.post_mortem_interval_in_hours)
samples["ph"] = pd.to_numeric(samples.brain_ph, errors="coerce")
samples["rin"] = pd.to_numeric(samples.rin, errors="coerce")
assert table.shape == (26493, 1231) and samples.gsm.notna().all() and table.notna().all().all()
assert samples.sample_id.is_unique and not samples.duplicated(["donor_id", "tissue"]).any()
assert samples.donor_id.nunique() == 134 and samples.tissue.nunique() == 10 and set(samples.arm) == {"aging"}
assert samples.groupby("donor_id")[["sex", "age_years", "brain_bank"]].nunique().eq(1).all().all()
assert len(genes) == 24885 and genes.entrez.str.isdigit().all()
columns = "column sample_id donor_id tissue sex age_years arm diagnosis brain_bank pmi_h ph rin".split()
save_rna(dataset, table.loc[genes.index].set_axis(genes.entrez.to_numpy()), samples[columns], "entrez", genes.symbol)
