# Pseudobulk counts per donor and cell population
# reads data/raw (00_download.sh); writes output/prepared/<dataset>/pseudobulk*
import glob

import anndata
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

import utils

RAW = utils.DATA / "raw"


def pseudobulk(blocks, populations, donor, min_cells):
    # sum counts of nuclei per donor and population; samples need min_cells nuclei
    sums, n_cells, annotation = {}, {}, []
    for counts, cells, genes in blocks:
        annotation.append(genes)
        for name, mask in populations(cells).items():
            for who, members in cells[mask].groupby(donor).groups.items():
                key = "|".join([*(who if isinstance(who, tuple) else [who]), name])
                total = pd.Series(
                    np.asarray(counts[cells.index.get_indexer(members)].sum(axis=0)).ravel(), index=genes.index
                )
                sums[key] = total if key not in sums else total.add(sums[key], fill_value=0)
                n_cells[key] = n_cells.get(key, 0) + len(members)
    keep = [key for key in sums if n_cells[key] >= min_cells]
    matrix = pd.DataFrame({key: sums[key] for key in keep}).fillna(0).astype("int64")
    assert matrix.min().min() >= 0
    genes = pd.concat(annotation).groupby(level=0, sort=False).first().loc[matrix.index]
    return matrix, [n_cells[key] for key in keep], genes


def labels(column, neurons, names=None):
    # population per label, plus the neuron union
    def populations(cells):
        found = names or cells[column].dropna().unique()
        return {**{label: cells[column] == label for label in found}, "neuron": cells[column].isin(neurons)}

    return populations


def save_pseudobulk(dataset, matrix, samples, genes):
    # counts, features and samples of a dataset
    features = genes[["ensembl", "symbol"]].rename_axis("gene_id")
    assert matrix.index.equals(features.index) and samples.sample_id.is_unique
    folder = utils.PREPARED / dataset
    utils.write(matrix.set_axis(samples.sample_id.to_numpy(), axis=1), folder / "pseudobulk.parquet")
    utils.write(features, folder / "pseudobulk_features.tsv")
    utils.write(samples, folder / "pseudobulk_samples.tsv")


def read_h5ad(path):
    # raw counts, nuclei and genes of an h5ad file
    adata = anndata.read_h5ad(path)
    raw = adata.raw if adata.raw is not None else adata
    names = raw.var_names.astype(str)
    counts = raw.X.tocsr() if sp.issparse(raw.X) else sp.csr_matrix(raw.X)
    genes = pd.DataFrame({"ensembl": np.nan, "symbol": names}, index=names)
    return counts, adata.obs.astype(str).reset_index(drop=True), genes


def read_h5(path):
    # 10x h5 as nuclei x genes, expression only
    with h5py.File(path, "r") as h5:
        group = h5["matrix"]
        shape = tuple(group["shape"][:])
        counts = sp.csc_matrix((group["data"][:], group["indices"][:], group["indptr"][:]), shape=shape)
        barcodes = pd.Index(group["barcodes"][:].astype(str))
        features = {key: group["features"][key][:].astype(str) for key in group["features"]}
    ids = features.get("id", features["name"])
    genes = pd.DataFrame({"ensembl": ids if "id" in features else np.nan, "symbol": features["name"]}, index=ids)
    expression = features["feature_type"] == "Gene Expression"
    return counts[np.where(expression)[0]].T.tocsr(), barcodes, genes[expression]


def read_10x(cells, files, key, barcode):
    # one block per h5 file with its labelled nuclei
    for value, rows in cells.groupby(key, sort=True):
        counts, barcodes, genes = read_h5(glob.glob(files.format(key=value))[0])
        rows = rows.assign(_position=barcodes.get_indexer(rows[barcode]))
        rows = rows[rows._position >= 0].reset_index(drop=True)
        yield counts[rows._position.to_numpy()], rows, genes


def fields(matrix, pattern=r"^(?P<donor>[^|]+)\|(?P<cell_type>.+)$"):
    # sample table from the <donor>|<population> column names
    samples = pd.DataFrame({"column": matrix.columns})
    return samples.join(samples.column.str.extract(pattern))


# Sloan 2026 prefrontal cortex (h5ad)
dataset = "sloan2026_snrna"
block = read_h5ad(RAW / "sloan2026_snt456vhrv853/zenodo/adata_snrna_rawcounts_0807_filtered.h5ad")
matrix, n_cells, genes = pseudobulk([block], labels("cc_broadclass", ["Exc", "Inh"]), "donor", 20)
meta = pd.read_csv(RAW / "sloan2026_snt456vhrv853/supp/zenodo_donor_metadata_sheet.csv")
samples = fields(matrix).merge(meta, left_on="donor", right_on="h5ad ID", how="left")
samples["tissue"] = "DLPFC"
samples["age_text"] = samples["*Age at Death"].astype(str).str.extract(r"^(\d+)")[0]
samples = samples.assign(
    sample_id=samples.column, donor_id=samples.donor, age_years=pd.to_numeric(samples.age_text, errors="coerce")
)
samples["sex"] = samples["*Sex"].replace(utils.SEX)
assert samples.age_years.notna().all()
samples["site"] = samples.Site
samples["n_cells"] = n_cells
save_pseudobulk(dataset, matrix, samples, genes)

# Zemke 2026 hippocampus (10x h5 per sample)
dataset = "zemke2026_gse278576"
geo = f"{RAW / dataset}/geo"
cells = pd.read_csv(
    f"{geo}/GSE278576_hippocampus_RNA_seurat_object_filtered_cells_metadata.tsv.gz", sep="\t", low_memory=False
)
cells["_barcode"] = cells.bacrode.str.extract(r"([ACGT]{16}-\d+)$")[0]  # column name as released
neurons = ["CA1", "CA2-CA3", "DG", "SUB", "LAMP5", "NR2F2", "PVALB", "SST", "VIP", "Chandelier"]
blocks = read_10x(cells, geo + "/GSE278576_{key}_raw_feature_bc_matrix.h5", "orig.ident", "_barcode")
matrix, n_cells, genes = pseudobulk(blocks, labels("subclass", neurons), "orig.ident", 20)
donors = cells.drop_duplicates("orig.ident")[["orig.ident", "Age", "Gender"]]
samples = fields(matrix).merge(donors, left_on="donor", right_on="orig.ident", how="left")
samples = samples.assign(
    tissue="Hippocampus",
    sample_id=samples.column,
    donor_id=samples.donor,
    age_years=pd.to_numeric(samples.Age, errors="coerce"),
)
assert samples.age_years.notna().all()
samples["sex"] = samples.Gender.replace(utils.SEX)
samples["n_cells"] = n_cells
save_pseudobulk(dataset, matrix, samples, genes)

# Handsaker 2025 caudate (NeMO h5 per reaction)
dataset = "handsaker2025_nemo_broadhdcaudate"
stage = f"{RAW / dataset}/stage2"
cells = pd.read_csv(f"{stage}/10x_village_cell_metadata_with_SPN_subtypes.txt", sep="\t", low_memory=False)
CLASSES = {
    "MSN": "SPN",
    "Astrocyte": "astrocyte",
    "Oligodendrocyte": "oligodendrocyte",
    "OPC": "polydendrocyte",
    "Microglia": "microglia",
    "Interneuron": "interneuron",
    "Endothelial": "endothelia",
}


def populations(cells):
    # broad classes, neurons and overlapping SPN subtypes
    out = {name: cells.predClass == label for name, label in CLASSES.items()}
    out["neuron"] = cells.predClass.isin(["SPN", "interneuron"])
    spn = cells.predClass == "SPN"
    regular = spn & (cells.predClass_SPN_EN == "noneccentric")
    direct, indirect = cells.predClass_SPN_DI == "direct", cells.predClass_SPN_DI == "indirect"
    out |= {
        "D1_MSN": regular & direct,
        "D2_MSN": regular & indirect,
        "eSPN": spn & (cells.predClass_SPN_EN == "eccentric"),
    }
    for pathway, selected in (("D1", direct), ("D2", indirect)):
        for zone in ("patch", "matrix"):
            out[f"{pathway}_MSN_{zone}"] = regular & selected & (cells.predClass_SPN_MP == zone)
    return out


blocks = read_10x(cells, stage + "/CHDI_{key}.*.donors.h5", "PREFIX", "CBC")
matrix, n_cells, genes = pseudobulk(blocks, populations, "DONOR", 10)
meta = pd.read_csv(f"{stage}/donor_metadata.txt", sep="\t", dtype=str)
samples = fields(matrix).merge(meta, left_on="donor", right_on="SID", how="left")
samples = samples.assign(tissue="Head of caudate", sample_id=samples.column, donor_id=samples.donor)
samples["age_years"] = pd.to_numeric(samples.Age.replace(">89", 90), errors="coerce")
samples["sex"] = samples.Sex.replace(utils.SEX)
samples["arm"] = samples.Status.replace({"Control": "aging", "Case": "disease"})
samples["pmi_h"] = pd.to_numeric(samples.PMI, errors="coerce")
samples["n_cells"] = n_cells
save_pseudobulk(dataset, matrix, samples, genes)

# Mesecar 2026 putamen, cortex, SVZ: donor x region x type
dataset = "mesecar2026_zenodo7803697"
counts, cells, genes = read_h5ad(RAW / dataset / "stage2/aging.pegasus.leiden_085.subclustered.h5ad")
broad = ["Astrocyte", "Endothelial", "Ependymal", "ExN", "InN", "Microglia", "Mural", "OPC", "Oligodendrocyte", "SPN"]
by_class = labels("broad_celltype", ["SPN", "InN", "ExN"], broad)
matrix, n_cells, genes = pseudobulk([(counts, cells, genes)], by_class, ["donor_id", "Brain_region"], 10)
samples = fields(matrix, r"^(?P<donor>[^|]+)\|(?P<region>[^|]+)\|(?P<cell_type>.+)$")
donors = cells[["donor_id", "Tissue_source", "Clinical_diagnosis", "Age", "Sex", "Age_group"]]
samples = samples.merge(donors.drop_duplicates("donor_id"), left_on="donor", right_on="donor_id", how="left")
samples = samples.assign(
    sample_id=samples.column,
    donor_id=samples.donor,
    tissue=samples.region,
    age_years=pd.to_numeric(samples.Age, errors="coerce"),
)
samples["sex"] = samples.Sex.replace(utils.SEX)
samples["n_cells"] = n_cells
save_pseudobulk(dataset, matrix, samples, genes)


def deposited(dataset, unions):
    # released pseudobulk; listed populations summed per donor
    folder = RAW / dataset / "pseudobulk"
    matrix = pd.read_parquet(folder / "matrix.parquet")
    meta = pd.read_csv(folder / "samples.tsv", sep="\t", dtype={"sample_id": str}).set_index("sample_id")
    genes = pd.read_csv(folder / "genes.tsv", sep="\t", keep_default_na=False, na_values=[""], dtype={"gene_id": str})
    for name, members in unions.items():
        part = meta[meta.cell_type.isin(members)]
        sums = matrix[part.index].T.groupby(part.animal_id.to_numpy()).sum().T
        added = part.groupby("animal_id", as_index=False).agg("first")
        added = added.assign(cell_type=name, n_cells=part.groupby("animal_id").n_cells.sum().to_numpy())
        matrix = pd.concat([matrix, sums.add_suffix("|" + name)], axis=1)
        meta = pd.concat([meta, added.set_index(added.animal_id + "|" + name)])
    samples = pd.DataFrame({"column": matrix.columns}).merge(
        meta.reset_index(names="key"), left_on="column", right_on="key", how="left"
    )
    genes = pd.DataFrame(
        {"ensembl": genes.ensembl.to_numpy(), "symbol": genes.gene_symbol_src.to_numpy()},
        index=pd.Index(genes.gene_id, name="gene_id"),
    )
    assert samples.key.notna().all()
    return matrix, samples.assign(sample_id=samples.column, donor_id=samples.animal_id), genes.loc[matrix.index]


# Hahn 2023 mouse caudate putamen; MSN and neuron sums
dataset = "hahn2023_gse212576"
unions = {"MSN": ["MSN_D1", "MSN_D2"], "neuron": ["MSN_D1", "MSN_D2", "eSPN", "interneuron"]}
matrix, samples, genes = deposited(dataset, unions)
samples["batch"] = samples.library
samples["tissue"] = samples.region
save_pseudobulk(dataset, matrix, samples, genes)

# Green 2025 human caudate in alcohol use disorder
dataset = "green2025_gse277313"
matrix, samples, genes = deposited(dataset, {"neuron": ["D1_MSN", "D2_MSN", "eSPN", "MSN_unassigned", "Interneuron"]})
samples["tissue"] = samples.region
save_pseudobulk(dataset, matrix, samples, genes)

# Allen aging mouse atlas (Jin 2025), released pseudobulk
dataset = "jin2025_zengagingmouse10xv3"
save_pseudobulk(dataset, *deposited(dataset, {}))
