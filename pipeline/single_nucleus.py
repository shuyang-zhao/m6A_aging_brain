# Single-nucleus units of the coupling model for script 38 (IGF2BP2 in striatal populations): every population of Green
# 2025 and Handsaker 2025 (caudate) with at least 18 donors that have at least 20 nuclei is a unit of donor pseudobulks
# (log2 CPM, genes with 1 CPM in half of the samples; aging arm only), in the form of script 33's unit data
# reads the donor pseudobulks of script 03; writes nothing (imported)
import numpy as np
import pandas as pd

import model
import network as net
import utils

MIN_DONORS, MIN_NUCLEI = 18, 20
# dataset: region column, population column, numeric and categorical covariates of the unit
DATASETS = {
    "green2025_gse277313": dict(region="region", population="cell_type", numeric=["pmi_h"],
                                categorical=["sex", "batch"]),
    "handsaker2025_nemo_broadhdcaudate": dict(region="tissue", population="cell_type", numeric=["pmi_h"],
                                              categorical=["sex"]),
}


def find_units():
    """Populations that meet the inclusion rule: (dataset, region, population, sample ids)."""
    units = []
    for dataset, spec in DATASETS.items():
        samples = utils.read(utils.PREPARED / dataset / "pseudobulk_samples.tsv", low_memory=False,
                             dtype={"sample_id": str})
        if "arm" in samples:
            samples = samples[samples.arm == "aging"]
        samples = samples[samples.n_cells >= MIN_NUCLEI]
        region = samples[spec["region"]].astype(str).str.split(" | ", regex=False).str[0]
        for (reg, population), rows in samples.groupby([region, samples[spec["population"]].astype(str)]):
            if rows.donor_id.nunique() >= MIN_DONORS and rows.donor_id.is_unique:
                units.append(dict(dataset=dataset, region=reg, population=population, species="human",
                                  unit=f"{dataset}|{reg}|{population}", sample_ids=rows.sample_id.tolist()))
    return units


def to_ensembl(logcpm, folder, stem="pseudobulk"):
    """Rows of a pseudobulk matrix by Ensembl id (the id of the features table, else the current symbol's); duplicated
    ids dropped."""
    features = utils.read(folder / f"{stem}_features.tsv", index_col=0, dtype={0: str}).reindex(logcpm.index)
    names = utils.ensg_to_symbol("human")
    ensembl = pd.Series(names.index, index=names.to_numpy())
    ensembl = ensembl[~ensembl.index.duplicated(keep=False)]
    symbol = pd.Series(utils.current_symbol(features.symbol, "human", features.ensembl), index=features.index)
    gene = features.ensembl.where(features.ensembl.notna() & (features.ensembl != ""), symbol.map(ensembl))
    logcpm = logcpm[gene.notna().to_numpy()].set_axis(gene.dropna().to_numpy(), axis=0)
    return logcpm[~logcpm.index.duplicated(keep=False)]


def logcpm_of(dataset, sample_ids):
    """log2 CPM of the unit (Ensembl identifiers) for the genes expressed at 1 CPM in half of the samples."""
    counts = utils.read(utils.PREPARED / dataset / "pseudobulk.parquet", columns=sample_ids)
    cpm = counts.div(counts.sum(), axis=1) * 1e6
    logcpm = np.log2((counts + 0.5).div(counts.sum() + 1, axis=1) * 1e6)[(cpm >= 1).mean(axis=1) >= 0.5]
    return to_ensembl(logcpm, utils.PREPARED / dataset)


def load_unit(unit, controls):
    """Data of a unit in the form of script 33: mRNA matrix of the genes with annotation, design, regulators, strata."""
    spec = DATASETS[unit["dataset"]]
    samples = utils.read(utils.PREPARED / unit["dataset"] / "pseudobulk_samples.tsv", low_memory=False,
                         dtype={"sample_id": str})
    samples = samples.set_index("sample_id").loc[unit["sample_ids"]]
    samples["age"] = samples.age_years / 10
    logcpm = logcpm_of(unit["dataset"], unit["sample_ids"])
    numeric = [c for c in spec["numeric"] if samples[c].notna().all() and samples[c].nunique() > 1]
    categorical = [c for c in spec["categorical"]
                   if samples[c].notna().all() and 1 < samples[c].nunique() <= len(samples) / 5]
    X = model.design(samples, dict(numeric=numeric, categorical=categorical))
    level = logcpm.mean(axis=1)
    genes = net.features["human"].drop(columns="level", errors="ignore").join(level.rename("level"), how="inner")
    genes = genes[genes[net.COVARIATES].notna().all(axis=1)]
    regulators, own = {}, {}
    for name in net.READERS + net.WRITERS + net.ERASERS + controls:
        gene = net.symbols["human"].get(name)
        if gene in logcpm.index and logcpm.loc[gene].std() > 0:
            regulators[name], own[name] = net.z(logcpm.loc[gene].to_numpy(float)), gene
    return dict(unit=unit["unit"], layer="rna_only", species="human", samples=samples, X=X.to_numpy(float), genes=genes,
                Y=logcpm.loc[genes.index].to_numpy(float), M=None, load=unit["load"],
                covariates=np.column_stack([net.z(genes[c]) for c in net.COVARIATES]), regulators=regulators, own=own,
                strata=net.strata_of(samples))
