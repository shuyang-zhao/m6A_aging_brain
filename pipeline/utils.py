# Shared helpers: paths, tables, prepared units, gene identifiers and a few statistics
from functools import cache
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "output"
PREPARED = OUT / "prepared"
RESOURCES = PREPARED / "resources"

SEED = 20261001  # seed of the permutations, the matched gene sampling and the models
N_PERMUTATIONS = 1000
Q = 0.05  # BH q-value cutoff
Z95 = 1.96  # rounded normal quantile of the 95% intervals computed here (limma and statsmodels give their own)

MODIFICATIONS = ["m6A", "Psi", "m5C", "m7G"]
PEAK_WINDOW = 50  # nt from a CLIP binding region to a single-base site of the reader's modification
# m6A-seq peaks of five human brain regions give the m6A load of a human transcript (load models and coupling model),
# not sites; the single-base m6A maps give the sites
BRAIN_M6A = "m6a_brain_shafik2026"
SEX = {"M": "male", "F": "female", "MALE": "male", "FEMALE": "female", "Male": "male", "Female": "female"}
HUMAN_PROTEIN_UNITS = [  # human units with a protein layer and no mRNA layer
    "koopmans2026_dpfc_l1_3",
    "koopmans2026_dpfc_l4_6",
    "anderton2026_hippocampus",
    "eggers2025_substantia_nigra",
    "lee2023_hippocampus_ca1",
    "ping2020_dlpfc",
    "teng2023_amygdala",
    "teng2023_hippocampus",
    "teng2023_hypothalamus",
    "teng2023_nucleus_accumbens",
    "teng2023_prefrontal_cortex",
    "teng2023_ventral_tegmental_area",
    "merrihew2023_smtg",
    "merrihew2023_hippocampus",
    "merrihew2023_caudate",
    "merrihew2023_ipl",
    "blumenreich2024_striatum",
    "blumenreich2024_substantia_nigra",
    "blumenreich2024_cingulate_gyrus",
    "blumenreich2024_middle_temporal_gyrus",
    "blumenreich2024_occipital_cortex",
]
# context, layer, units; each context and layer is one family of tests
# eight units fit no context and belong to none: GTEx amygdala, hypothalamus, spinal cord C1 and substantia nigra,
# UKBEC medulla, substantia nigra, thalamus and white matter
CONTEXTS = [
    (
        "human_striatum",
        "mrna_only",
        "gtex_caudate gtex_putamen gtex_nucleus_accumbens brainseq3_caudate markunas2020_nac ramasamy2014_putamen",
    ),
    (
        "human_cortex",
        "mrna_only",
        "gtex_cortex gtex_frontal_cortex_ba9 gtex_anterior_cingulate_ba24 brainseq2_dlpfc luo2024_bankadj "
        "ramasamy2014_frontal_cortex ramasamy2014_temporal_cortex ramasamy2014_occipital_cortex",
    ),
    ("human_hippocampus", "mrna_only", "gtex_hippocampus brainseq2_hippocampus ramasamy2014_hippocampus"),
    ("human_cerebellum", "mrna_only", "gtex_cerebellum gtex_cerebellar_hemisphere ramasamy2014_cerebellar_cortex"),
    (
        "human_neurons",
        "mrna_only",
        "sloan2026_dlpfc_neurons zemke2026_hippocampus_neurons handsaker2025_caudate_neurons "
        "matlik2024_neurons green2025_caudate_neurons mesecar2026_putamen_neurons",
    ),
    ("human_protein", "protein_only", " ".join(["luo2024_bankadj", *HUMAN_PROTEIN_UNITS])),
    ("human_protein", "mrna_adjusted", "luo2024_bankadj"),
    (
        "mouse_mrna",
        "mrna_only",
        "takasugi2024 chdi2016_rna hahn2023_bulk_motor_cortex hahn2023_bulk_visual_cortex "
        "hahn2023_bulk_entorhinal_cortex hahn2023_bulk_hippocampus_anterior hahn2023_bulk_hippocampus_posterior "
        "hahn2023_bulk_caudate_putamen hahn2023_bulk_olfactory_bulb hahn2023_bulk_hypothalamus hahn2023_bulk_thalamus "
        "hahn2023_bulk_pons hahn2023_bulk_medulla hahn2023_bulk_cerebellum hahn2023_bulk_corpus_callosum "
        "hahn2023_bulk_choroid_plexus hahn2023_bulk_subventricular_zone shavlakadze2026_cortex "
        "shavlakadze2026_hippocampus shavlakadze2026_hypothalamus shavlakadze2026_cerebellum",
    ),
    ("mouse_neurons", "mrna_only", "hahn2023_msn hahn2023_neurons"),
    (
        "mouse_protein",
        "protein_only",
        "chdi2016_expanded takasugi2024 lu2025_brainstem lu2025_cerebellum lu2025_hippocampus lu2025_olfactory_bulb "
        "tsumagari2023_cortex tsumagari2023_hippocampus wang2025_brain keele2025_cortex keele2025_hippocampus "
        "keele2025_striatum marino2025_whole_brain",
    ),
    ("mouse_protein", "mrna_adjusted", "chdi2016_expanded takasugi2024"),
    # sensitivity: the mouse proteomes of cohorts that reach 20 months or more (Langfelder 2016, 2-10 months, left out)
    (
        "mouse_protein_aged",
        "protein_only",
        "takasugi2024 lu2025_brainstem lu2025_cerebellum lu2025_hippocampus lu2025_olfactory_bulb tsumagari2023_cortex "
        "tsumagari2023_hippocampus wang2025_brain keele2025_cortex keele2025_hippocampus keele2025_striatum "
        "marino2025_whole_brain",
    ),
]
# covariates of the load models (scripts 21, 25 and 28) and of stage 2 of the coupling model: transcript, genomic,
# intronic and 3' UTR length, mean level in the layer (mRNA in the mRNA-adjusted layer), neuron share, programmes
LOAD_COVARIATES = [
    "log_tx_len",
    "log_gene_len",
    "log1p_intron_len",
    "log1p_utr3_len",
    "level",
    "ct_frac_neuron",
    "myelin",
    "lysosome",
    "mitochondrion",
    "synapse",
]
# repeats units of mouse_protein: a unit's context is its first one, and tests pooled over contexts leave it out
SENSITIVITY = "mouse_protein_aged"
# single-nucleus cohorts: sample filter of the neuron unit, neuron and non-neuronal cell-type labels
SINGLE_NUCLEUS = {
    "sloan2026_snrna": dict(keep="tissue == 'DLPFC'", neuron="neuron", other="Ast,Oli,Mic,OPC,End,Peri"),
    "zemke2026_gse278576": dict(
        keep="tissue == 'Hippocampus'", neuron="neuron", other="Astro,Oligo,Microglia,OPC,Endo,Macro,VLMC,T-Cell"
    ),
    "handsaker2025_nemo_broadhdcaudate": dict(
        keep="arm == 'aging'", neuron="neuron", other="Astrocyte,Oligodendrocyte,OPC,Microglia,Endothelial"
    ),
    "matlik2024_gse227729": dict(
        keep="arm == 'aging'",
        neuron="D1 MSN,D2 MSN,Interneuron TAC3+,Interneuron FS,Interneuron SST+,Interneuron cholinergic",
        other="Microglia,Astrocyte,Oligodendrocyte",
    ),
    "green2025_gse277313": dict(
        keep="arm == 'aging'",
        neuron="neuron",
        other="Oligodendrocyte,OPC,Astrocyte,Endothelial,Microglia,Ependymal,Pericyte,Macrophage,Vascular",
    ),
    "mesecar2026_zenodo7803697": dict(
        keep="tissue == 'Putamen'",
        neuron="neuron",
        other="Astrocyte,Endothelial,Ependymal,Microglia,Mural,OPC,Oligodendrocyte",
    ),
}
BATCH_CORRECTED = ["luo2024_bankadj"]  # TMT batch removed from protein as a random effect (scripts 10 and 12)
# human units code age per decade, mouse units the late against the early age group of the unit
AGE_SCALE = {"human": "per decade", "mouse": "late vs early"}


def read(path, **kwargs):
    """Read a tsv, tsv.gz or parquet table."""
    if str(path).endswith(".parquet"):
        return pd.read_parquet(path, **kwargs)
    return pd.read_csv(path, sep="\t", **kwargs)


def write(df, path):
    """Write a table by its suffix; text files keep the index only when it is named."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".parquet":
        df.to_parquet(path)
    else:
        df.to_csv(path, sep="\t", index=any(name is not None for name in df.index.names))
    return path


def load_unit(unit):
    """Protein and mRNA matrices (None if absent), sample table and model design of a unit."""
    folder = PREPARED / unit
    samples = read(folder / "samples.tsv", index_col="sample_id", dtype={"sample_id": str})
    layers = [folder / f"{layer}.parquet" for layer in ("protein", "rna")]
    protein, rna = [read(path) if path.exists() else None for path in layers]
    row = read(PREPARED / "units.tsv", index_col="unit", keep_default_na=False).loc[unit]
    design = {
        "numeric": [c for c in row.numeric.split(",") if c],
        "categorical": [c for c in row.categorical.split(",") if c],
        "hidden_factors": int(row.hidden_factors),
    }
    return protein, rna, samples, design


def ensg_to_symbol(species):
    """Series mapping unversioned Ensembl gene ids to symbols."""
    return read(PREPARED / "ref" / f"gene_symbol_{species}.tsv").set_index("gene_id").symbol


def orthologs():
    """One-to-one mouse-human ortholog pairs (columns mouse, human; upper-case symbols)."""
    return read(PREPARED / "ref" / "orthologs.tsv")


@cache
def symbol_steps(species):
    """Name lookups of output/prepared/ref/symbols_<species>.tsv (script 06): one Series per step, indexed by the
    queried name."""
    table = read(PREPARED / "ref" / f"symbols_{species}.tsv", keep_default_na=False)
    return {step: rows.set_index("query").value for step, rows in table.groupby("step")}


def current_symbol(symbol, species, ensembl):
    """Current symbols by Ensembl id, official name or old symbol; an old symbol never takes a listed gene's name."""
    steps = symbol_steps(species)
    names = pd.Series(symbol, dtype=object).fillna("").astype(str)
    direct = names.map(steps["official"]).fillna(names.str.upper().map(steps.get("official_upper", {})))
    direct = pd.Series(ensembl, dtype=object).map(steps["ensembl"]).fillna(direct)
    synonym = names.map(steps["synonym"]).fillna(names.str.upper().map(steps["synonym_upper"]))
    resolved = direct.fillna(synonym).fillna(names.replace("", np.nan))
    shadowed = direct.isna() & synonym.notna() & resolved.isin(set(direct.dropna()))
    return resolved.where(~shadowed, names).to_numpy()


def read_gtf(path, features):
    """GENCODE GTF rows of the given features with gene and transcript ids, names, types and tags."""
    names = ["chrom", "source", "feat", "start", "end", "score", "strand", "frame", "attr"]
    columns = [0, 2, 3, 4, 6, 8]
    dtype = {"chrom": str, "strand": str}
    gtf = pd.read_csv(path, sep="\t", comment="#", header=None, names=names, usecols=columns, dtype=dtype)
    gtf = gtf[gtf.feat.isin(features)].reset_index(drop=True)
    for name, key in (
        ("gene_id", "gene_id"),
        ("gene_name", "gene_name"),
        ("gene_type", "gene_type"),
        ("tx_id", "transcript_id"),
        ("tx_type", "transcript_type"),
    ):
        gtf[name] = gtf.attr.str.extract(rf'{key} "([^"]*)"')[0]
    gtf["tags"] = gtf.attr.str.findall(r'tag "([^"]*)"').str.join(";")
    return gtf.drop(columns="attr")


def lift(chain, chrom, pos, strand):
    """Lift 1-based positions with a pyliftover chain; unmapped or cross-chromosome sites are NaN."""
    strand = np.broadcast_to(np.asarray(strand, dtype=object), len(pos))
    lifted = []
    for contig, position, sign in zip(chrom, pos, strand, strict=True):
        hit = chain.convert_coordinate(contig, int(position) - 1, sign)
        lifted.append(hit[0][:3] if hit and hit[0][0] == contig else (np.nan,) * 3)
    chroms, positions, strands = zip(*lifted, strict=True)
    return pd.DataFrame({"chrom": chroms, "pos": [position + 1 for position in positions], "strand": strands})


def bh(p):
    """Benjamini-Hochberg q-values; NaN p-values stay NaN."""
    values = np.asarray(p, float)
    q = np.full(values.shape, np.nan)
    ok = np.isfinite(values)
    if ok.any():
        q[ok] = sps.false_discovery_control(values[ok])
    return pd.Series(q, index=p.index, name=p.name) if isinstance(p, pd.Series) else q


def pct(log2fc):
    """Percentage change for a log2 fold change."""
    return 100 * np.expm1(np.log(2) * log2fc)


def random_effects(est, se):
    """DerSimonian-Laird random-effects pooling: estimate, se, 95% CI, z, p, tau2, I2 (percent) and cohort count."""
    est, var = np.asarray(est, float), np.asarray(se, float) ** 2
    weight, k = 1 / var, len(est)
    q = (weight * (est - (weight * est).sum() / weight.sum()) ** 2).sum()
    tau2 = i2 = np.nan  # not estimable from one cohort
    if k > 1:
        tau2 = max(0.0, (q - (k - 1)) / (weight.sum() - (weight**2).sum() / weight.sum()))
        i2 = 100 * max(0.0, 1 - (k - 1) / q) if q > 0 else 0.0
    weight = 1 / (var + np.nan_to_num(tau2))
    pooled = (weight * est).sum() / weight.sum()
    error = 1 / np.sqrt(weight.sum())
    z = pooled / error
    p = 2 * sps.norm.sf(abs(z))
    ci = dict(ci_low=pooled - Z95 * error, ci_high=pooled + Z95 * error)
    return dict(est=pooled, se=error, **ci, z=z, p=p, tau2=tau2, i2=i2, k=k)


def stouffer(z, w):
    """Weighted Stouffer combination of z-scores: combined z and two-sided p."""
    z = np.asarray(z, float)
    w = np.asarray(w, float)
    combined = (w * z).sum() / np.sqrt((w**2).sum())
    return dict(z=combined, p=2 * sps.norm.sf(abs(combined)))


def by_symbol(table):
    """Table indexed by upper-case symbol, ambiguous symbols dropped."""
    table = table.dropna(subset=["symbol"]).assign(key=lambda d: d.symbol.str.upper())
    return table[~table.key.duplicated(keep=False)].set_index("key")


def near_peaks(sites, peaks, window):
    """Sites (1-based atlas positions) within window nt of a peak on the same chromosome and strand."""
    near = pd.Series(False, index=sites.index)
    by_strand = dict(list(peaks.sort_values("start").groupby(["chrom", "strand"])))
    for key, group in sites.groupby(["chrom", "strand"]):
        if key not in by_strand:
            continue
        peak, base = by_strand[key], group.pos.to_numpy() - 1
        reach = np.maximum.accumulate(peak.end.to_numpy() - 1 + window)  # last base within reach of the peaks so far
        last = np.searchsorted(peak.start.to_numpy() - window, base, side="right") - 1
        near[group.index] = (last >= 0) & (reach[np.maximum(last, 0)] >= base)
    return near


def perm_pvalue(obs, null, alternative="two-sided"):
    """Permutation p-value, (1 + null draws at least as extreme) / (1 + draws); two-sided or greater."""
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    extreme = np.abs(null) >= abs(obs) if alternative == "two-sided" else null >= obs
    return (1 + extreme.sum()) / (1 + len(null))
