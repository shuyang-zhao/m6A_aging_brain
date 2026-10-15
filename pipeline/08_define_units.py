# Analysis units: samples, age coding, covariates, blocks, hidden factors
# and gene-level matrices (protein median-normalised, RNA TMM log-CPM)
# reads 01-03 and ref/ of 06; writes output/prepared/<unit>/ and units.tsv
import numpy as np
import pandas as pd
from scipy import stats
from scipy.linalg import null_space

import utils

MIN_DONORS = 8  # units with fewer donors are not written
MARKERS = {  # marker genes of the cell-type indices
    "rna": {
        "neuron": "RBFOX3 SNAP25 SYT1 SLC17A7",
        "astro": "ALDH1L1 AQP4 SLC1A3 GLUL",
        "oligo": "MBP MOG PLP1 MAG",
        "micro": "P2RY12 TMEM119 CX3CR1 AIF1",
        "vascular": "CLDN5 PECAM1 KDR EMCN",
    },
    "protein": {
        "Astro": "GFAP AQP4 ALDH1L1 SLC1A3 SLC1A2 GJA1",
        "Micro": "CX3CR1 P2RY12 HEXB CSF1R C1QA TMEM119 AIF1",
        "Neuron": "SNAP25 SYT1 STMN2 NEFL NEFM SYN1",
        "Oligo": "MBP PLP1 MOBP MOG MAG CNP CLDN11",
    },
}
MOUSE = ("takasugi", "langfelder", "lu2025", "tsumagari", "wang2025", "keele", "marino", "hahn", "shavlakadze")
CELLS = "cell_neuron cell_astro cell_oligo cell_micro"
PROTEIN_CELLS = ("protein", "Neuron Oligo Astro Micro", 3)
ADULT = "age_years >= 18"
AGE10 = dict(column="age_years", per=10)  # age per ten years

# unit -> dataset, sample queries, donor pairing, columns,
# age coding, covariates, permutation blocks, hidden factors
# cols: name=[log1p:][layer.]source
# pair: column=equal or a tolerance
# cells: layer, wanted indices, minimum markers
# fmt: off
UNITS = {
    "takasugi2024": dict(
        data="takasugi2024_pxd040722", protein="tissue == 'Cerebrum'", rna="tissue == 'Cerebrum'",
        value="norm_counts", pair="sex=equal age_months=12", min_side=2, cells=("rna", None, 2),
        cols="sex=protein.sex age_months=protein.age_months protein_batch=protein.batch rna_batch=rna.batch",
        age=dict(column="age_months", early=[6], late=[30], middle=[15, 24]),
        numeric=f"middle_15 middle_24 {CELLS} cell_vascular", categorical="sex", blocks="sex"),
    "chdi2016_expanded": dict(
        data="langfelder2016_pxd006302", hidden=1, min_side=2, cells=("rna", None, 2),
        protein="tissue == 'Striatum' and age_months in [2, 10]", rna="tissue == 'Striatum' and age_months in [2, 10]",
        pair="sex=equal age_months=equal genotype=equal", groups=("genotype", "age_months", 2),
        cols="sex=protein.sex age_months=protein.age_months protein_batch=protein.batch rna_batch=rna.batch "
             "genotype=protein.genotype arm=protein.arm",
        age=dict(column="age_months", early=[2], late=[10]),
        numeric=CELLS, categorical="sex genotype", blocks="sex genotype"),
    "chdi2016_rna": dict(
        data="langfelder2016_pxd006302", hidden=1, min_side=2, cells=("rna", None, 2), donors_of="chdi2016_expanded",
        rna="tissue == 'Striatum' and age_months in [2, 10]",
        cols="sex age_months rna_batch=batch genotype arm", age=dict(column="age_months", early=[2], late=[10]),
        numeric=CELLS, categorical="sex genotype", blocks="sex genotype"),
    "luo2024_bankadj": dict(
        data="luo2024_syn32136022", hidden=2, value="log2", cells=("rna", None, 2),
        protein="tissue == 'Frontal cortex' and arm == 'aging'",
        rna=f"tissue == 'Frontal cortex' and arm == 'aging' and {ADULT}",
        pair="sex=equal age_years=0.5 brain_bank=equal",
        cols="sex=protein.sex age_years=protein.age_years protein_batch=protein.batch rna_batch=rna.batch "
             "brain_bank=protein.brain_bank rin=rna.rin log_pmi=log1p:rna.pmi_h",
        age=AGE10, numeric=f"rin log_pmi {CELLS}", categorical="brain_bank sex rna_batch", blocks="sex protein_batch"),
    "brainseq3_caudate": dict(
        data="benjamin2022_brainseq3_caudate", rna=f"arm == 'aging' and {ADULT}", complete="rin log_pmi ph",
        cells=("rna", None, 2), cols="sex age_years batch rin log_pmi=log1p:pmi_h ph", age=AGE10,
        numeric=f"rin log_pmi ph {CELLS}", categorical="sex batch", blocks="sex batch"),
}
for layers, name in (("L1-3", "l1_3"), ("L4-6", "l4_6")):
    UNITS[f"koopmans2026_dpfc_{name}"] = dict(
        data="koopmans2026_pxd058441", protein=f"tissue == 'DPFC {layers}' and cohort == 'control' and {ADULT}",
        complete="log_pmd", cells=PROTEIN_CELLS,
        cols="sex age_years batch brain_bank gel_lane log_pmd=log1p:pmd_h ms_run_order", age=AGE10,
        numeric="cell_Neuron cell_Oligo cell_Astro log_pmd ms_run_order", categorical="sex batch brain_bank gel_lane",
        blocks="sex batch")
UNITS["eggers2025_substantia_nigra"] = dict(
    data="eggers2025_pxd051145", protein=f"tissue == 'Substantia nigra' and {ADULT} and pmi_h.notna()",
    cells=PROTEIN_CELLS, cols="sex age_years batch log_pmi=log1p:pmi_h ms_run_order", age=AGE10,
    numeric="log_pmi cell_Neuron cell_Oligo cell_Astro cell_Micro ms_run_order", categorical="sex", blocks="sex")
UNITS["anderton2026_hippocampus"] = dict(
    data="anderton2026_biorxiv714902", protein=f"tissue == 'Hippocampus' and {ADULT}", cells=PROTEIN_CELLS,
    cols="sex age_years batch biobank log_pmi=log1p:pmi_h", age=AGE10,
    numeric="log_pmi cell_Oligo cell_Astro", categorical="batch biobank", blocks="sex batch")
for tissue in ("Brainstem", "Cerebellum", "Hippocampus", "Olfactory bulb"):
    UNITS[f"lu2025_{tissue.lower().replace(' ', '_')}"] = dict(
        data="lu2025_pxd053154", protein=f"tissue == '{tissue}' and arm == 'aging' and age_months >= 4",
        cols="sex age_months", age=dict(column="age_months", early=[4, 7], late=[21, 24], middle=[10, 13, 17]),
        numeric="middle_10 middle_13 middle_17 early_1 late_1")
for tissue in ("Cortex", "Hippocampus"):
    UNITS[f"tsumagari2023_{tissue.lower()}"] = dict(
        data="tsumagari2023_pxd041485", protein=f"tissue == '{tissue}'", cols="sex age_months batch",
        age=dict(column="age_months", early=[3], late=[24], middle=[15]),
        numeric="middle_15", categorical="batch", blocks="batch")
UNITS["wang2025_brain"] = dict(
    data="wang2025_pxd047296", protein="tissue == 'Brain'", cols="sex age_months",
    age=dict(column="age_months", early=[4], late=[20], middle=[8, 12]),
    numeric="middle_8 middle_12", categorical="sex", blocks="sex")
for tissue in ("Cortex", "Hippocampus", "Striatum"):
    UNITS[f"keele2025_{tissue.lower()}"] = dict(
        data="keele2025_pxd058160", protein=f"tissue == '{tissue}'", cols="sex age_months",
        age=dict(column="age_months", early=[3], late=[20], middle=[12]),
        numeric="middle_12", categorical="sex", blocks="sex")
UNITS["marino2025_whole_brain"] = dict(
    data="marino2025_pxd047951", protein="tissue == 'Whole brain'", cols="sex age_months",
    age=dict(column="age_months", early=[3.5], late=[33.0]))
for tissue in ("DLPFC", "Hippocampus"):
    UNITS[f"brainseq2_{tissue.lower()}"] = dict(
        data="colladotorres2019_phs000979", rna=f"tissue == '{tissue}' and arm == 'aging' and {ADULT}",
        complete="rin mitorate totalassignedgene sex race", cols="sex age_years race rin mitorate totalassignedgene",
        age=AGE10, numeric="rin mitorate totalassignedgene", categorical="sex race", blocks="sex")
GTEX = {
    "amygdala": "Amygdala", "anterior_cingulate_ba24": "Anterior cingulate cortex (BA24)",
    "caudate": "Caudate (basal ganglia)", "cerebellar_hemisphere": "Cerebellar Hemisphere", "cerebellum": "Cerebellum",
    "cortex": "Cortex", "frontal_cortex_ba9": "Frontal Cortex (BA9)", "hippocampus": "Hippocampus",
    "hypothalamus": "Hypothalamus", "nucleus_accumbens": "Nucleus accumbens (basal ganglia)",
    "putamen": "Putamen (basal ganglia)", "spinal_cord_c1": "Spinal cord (cervical c-1)",
    "substantia_nigra": "Substantia nigra",
}
for region, tissue in GTEX.items():
    UNITS[f"gtex_{region}"] = dict(
        data="gtexv10_phs000424", rna=f"tissue == 'Brain - {tissue}'", complete="rin log_ischemic",
        cells=("rna", None, 2), age=dict(AGE10, column="age_midpoint"),
        cols="sex age_midpoint rin log_ischemic=log1p:ischemic_min dthhrdy center batch",
        numeric=f"rin log_ischemic {CELLS}", categorical="sex dthhrdy center", blocks="sex center")
UNITS["markunas2020_nac"] = dict(
    data="markunas2020_gse171936", rna="arm == 'aging'", complete="rin race sex",
    cols="sex=sex_letter age_years race rin", age=AGE10, numeric="rin", categorical="sex race", blocks="sex")


def neurons(data, query, cols, age=AGE10, **extra):
    # neuron pseudobulk unit per donor; sex is the only covariate
    age = dict(age, column=next(name for name in cols.split() if name.startswith("age")))
    return dict(data=data, rna=query, cols=cols, age=age, numeric="", categorical="sex", blocks="sex", **extra)


early_late = dict(early=[3], late=[21])
matlik = ["D1 MSN", "D2 MSN", "Interneuron TAC3+", "Interneuron FS", "Interneuron SST+", "Interneuron cholinergic"]
UNITS["sloan2026_dlpfc_neurons"] = neurons(
    "sloan2026_snrna", "tissue == 'DLPFC' and cell_type == 'neuron'", "sex age_years site")
UNITS["zemke2026_hippocampus_neurons"] = neurons(
    "zemke2026_gse278576", "tissue == 'Hippocampus' and cell_type == 'neuron'", "sex age_years")
UNITS["handsaker2025_caudate_neurons"] = neurons(
    "handsaker2025_nemo_broadhdcaudate", "tissue == 'Head of caudate' and cell_type == 'neuron' and arm == 'aging'",
    "sex age_years")
UNITS["matlik2024_neurons"] = neurons(  # the neuron populations and regions of a donor summed into one sample
    "matlik2024_gse227729", f"cell_type in {matlik} and arm == 'aging'", "sex age_years", pool=True)
UNITS["hahn2023_msn"] = neurons("hahn2023_gse212576", "cell_type == 'MSN'", "sex age_months batch", early_late)
UNITS["hahn2023_neurons"] = neurons("hahn2023_gse212576", "cell_type == 'neuron'", "sex age_months batch", early_late)
UNITS["green2025_caudate_neurons"] = neurons(
    "green2025_gse277313", "cell_type == 'neuron' and arm == 'aging'", "sex age_years batch")
UNITS["mesecar2026_putamen_neurons"] = neurons(
    "mesecar2026_zenodo7803697", "tissue == 'Putamen' and cell_type == 'neuron'", "sex age_years")
# human protein cohorts
# the donors the authors excluded are left out
UNITS["lee2023_hippocampus_ca1"] = dict(
    data="lee2023_pxd043880",
    protein=f"tissue == 'Hippocampus CA1' and arm == 'aging' and {ADULT} and author_excluded.isna()",
    complete="log_pmi", cells=PROTEIN_CELLS, cols="sex age_years brain_bank log_pmi=log1p:pmi_h ms_run_order",
    age=AGE10, numeric="log_pmi cell_Neuron cell_Oligo cell_Astro ms_run_order", categorical="sex brain_bank",
    blocks="sex")
# ten controls, batch and sex leave five residual df
UNITS["ping2020_dlpfc"] = dict(
    data="ping2020_pxd020296", protein=f"tissue == 'DLPFC' and arm == 'aging' and {ADULT}",
    cols="sex age_years batch", age=AGE10, categorical="sex batch", blocks="batch")
for tissue in (
    "Amygdala", "Hippocampus", "Hypothalamus", "Nucleus accumbens", "Prefrontal cortex", "Ventral tegmental area",
):
    micro = " cell_Micro" if tissue in ("Hippocampus", "Hypothalamus", "Prefrontal cortex") else ""
    UNITS[f"teng2023_{tissue.lower().replace(' ', '_')}"] = dict(
        data="teng2023_pxd040884", protein=f"tissue == '{tissue}' and arm == 'aging' and {ADULT}", complete="log_pmi",
        cells=PROTEIN_CELLS, cols="age_years batch log_pmi=log1p:pmi_h", age=AGE10,
        numeric="log_pmi cell_Neuron cell_Oligo cell_Astro" + micro, categorical="batch", blocks="batch")
for name, tissue in (
    ("smtg", "Superior and middle temporal gyri (SMTG)"),
    ("hippocampus", "Hippocampus"),
    ("caudate", "Caudate"),
):
    UNITS[f"merrihew2023_{name}"] = dict(
        data="merrihew2023_pxd034525", protein=f"tissue == '{tissue}' and arm == 'aging' and {ADULT}",
        complete="log_pmi", cols="sex age_years log_pmi=log1p:pmi_h", age=AGE10,
        numeric="log_pmi", categorical="sex", blocks="sex")
UNITS["merrihew2023_ipl"] = dict(
    data="merrihew2023_pxd034525", protein=f"tissue == 'Inferior parietal lobule (IPL)' and arm == 'aging' and {ADULT}",
    cols="sex age_years", age=AGE10, categorical="sex", blocks="sex")
for region in ("Striatum", "Substantia nigra", "Cingulate gyrus", "Middle temporal gyrus", "Occipital cortex"):
    UNITS[f"blumenreich2024_{region.lower().replace(' ', '_')}"] = dict(
        data="blumenreich2024_pxd047134", protein=f"tissue == '{region}' and arm == 'aging' and {ADULT}",
        cells=PROTEIN_CELLS, cols="sex age_years batch", age=AGE10,
        numeric="cell_Neuron cell_Oligo cell_Astro" + ("" if region == "Substantia nigra" else " cell_Micro"),
        categorical="sex batch", blocks="sex batch")
# bulk RNA cohorts
HAHN = {
    "motor_cortex": "Motor cortex", "visual_cortex": "Visual Cortex", "entorhinal_cortex": "Entorhinal cortex",
    "hippocampus_anterior": "Hippocampus (anterior)", "hippocampus_posterior": "Hippocampus (posterior)",
    "caudate_putamen": "Caudate putamen", "olfactory_bulb": "Olfactory bulb", "hypothalamus": "Hypothalamus",
    "thalamus": "Thalamus", "pons": "Pons", "medulla": "Medulla", "cerebellum": "Cerebellum",
    "corpus_callosum": "Corpus callosum", "choroid_plexus": "Choroid Plexus",
    "subventricular_zone": "Subventricular zone",
}
# the 26- and 28-month mice of the bulk series fall outside the age groups of the units (3 to 21 months)
for region, tissue in HAHN.items():
    UNITS[f"hahn2023_bulk_{region}"] = dict(
        data="hahn2023_gse212336", rna=f"tissue == '{tissue}' and age_months <= 21", cells=("rna", None, 2),
        cols="sex age_months", age=dict(column="age_months", early=[3], late=[21], middle=[12, 15, 18]),
        numeric=f"middle_12 middle_15 middle_18 {CELLS} cell_vascular", categorical="sex", blocks="sex")
for tissue in ("cortex", "hippocampus", "hypothalamus", "cerebellum"):
    UNITS[f"shavlakadze2026_{tissue}"] = dict(
        data="shavlakadze2026_gse314817", rna=f"tissue == '{tissue}'", cells=("rna", None, 2),
        cols="sex age_months", age=dict(column="age_months", early=[6], late=[27], middle=[12, 18, 21, 24]),
        numeric=f"middle_12 middle_18 middle_21 middle_24 {CELLS} cell_vascular", categorical="sex", blocks="sex")
UKBEC = (
    "frontal_cortex", "temporal_cortex", "occipital_cortex", "hippocampus", "medulla", "thalamus", "putamen",
    "white_matter", "cerebellar_cortex", "substantia_nigra",
)
for region in UKBEC:
    tissue = region.replace("_", " ")
    UNITS[f"ramasamy2014_{region}"] = dict(
        data="ramasamy2014_gse46706", rna=f"tissue == '{tissue}' and arm == 'aging' and {ADULT}", value="log2",
        complete="rin log_pmi", cells=("rna", None, 2), age=AGE10,
        cols="sex age_years brain_bank rin log_pmi=log1p:pmi_h",
        numeric=f"rin log_pmi {CELLS} cell_vascular", categorical="sex brain_bank", blocks="sex brain_bank")
# fmt: on


def tmm_factor(values, column, ref, lib):
    # TMM factor of a sample against the reference (edgeR)
    sample, reference = values[:, column], values[:, ref]
    with np.errstate(all="ignore"):
        ratio = np.log2((sample / lib[column]) / (reference / lib[ref]))
        level = (np.log2(sample / lib[column]) + np.log2(reference / lib[ref])) / 2
        variance = (lib[column] - sample) / lib[column] / sample + (lib[ref] - reference) / lib[ref] / reference
    finite = np.isfinite(ratio) & np.isfinite(level)
    ratio, level, variance = ratio[finite], level[finite], variance[finite]
    if np.abs(ratio).max() < 1e-6:
        return 1.0
    size = len(ratio)
    # edgeR defaults: 30% of the log ratios and 5% of the mean log expression trimmed at each end
    trim, trim_level = int(np.floor(size * 0.3)), int(np.floor(size * 0.05))
    keep = (stats.rankdata(ratio) >= trim + 1) & (stats.rankdata(ratio) <= size - trim)
    keep &= (stats.rankdata(level) >= trim_level + 1) & (stats.rankdata(level) <= size - trim_level)
    return 2 ** (np.sum(ratio[keep] / variance[keep]) / np.sum(1 / variance[keep]))


def log_expression(counts, value):
    # log2 expression of genes with 1 CPM in half the samples
    # value log2: already logged; norm_counts: library-size log-CPM
    # otherwise TMM log-CPM with a prior count of 0.5 (as limma voom; edgeR cpm uses 2)
    if value == "log2":
        return counts.astype(float)
    library = counts.sum(axis=0)
    expressed = (counts.div(library, axis=1) * 1e6 >= 1).mean(axis=1) >= 0.5
    if value == "norm_counts":
        return np.log2((counts.loc[expressed] + 0.5).div(library + 1, axis=1) * 1e6)
    values, lib = counts.to_numpy(float), library.to_numpy(float)
    values = values[(values > 0).any(axis=1)]
    # reference: the sample whose upper quartile is closest to the mean upper quartile
    quantile = np.quantile(values, 0.75, axis=0) / lib
    ref = np.argmin(np.abs(quantile - quantile.mean()))
    factors = np.array([tmm_factor(values, column, ref, lib) for column in range(len(lib))])
    effective = library * (factors / np.exp(np.mean(np.log(factors))))
    prior = 0.5 * effective / effective.mean()
    scaled = (counts.to_numpy(float) + prior.to_numpy()) / (effective + 2 * prior).to_numpy() * 1e6
    return pd.DataFrame(np.log2(scaled), index=counts.index, columns=counts.columns).loc[expressed]


def protein_genes(matrix, features, chosen, kept, samples, species):
    # median-normalised on proteins detected in 70% of the samples, then one protein group per Ensembl gene
    values = matrix[chosen]
    detected = values.notna().mean(axis=1) >= 0.7
    shift = values[detected].sub(values[detected].median(axis=1), axis=0).median()
    values = values.sub(shift, axis=1)[kept].set_axis(samples.index, axis=1)
    isoform = pd.Series(values.index).str.match(r"^[A-Z0-9]+-\d+(_|$)").fillna(False).to_numpy()
    rank = pd.DataFrame(
        {
            "symbol": features.gene_symbol.reindex(values.index),
            "detected": values.notna().mean(axis=1),
            "specific": ~features.genes.reindex(values.index).fillna("").str.contains(";"),
            "canonical": ~isoform,
            "mean": values.mean(axis=1),
        }
    )
    rank = rank[rank.symbol.notna() & ~rank.symbol.str.contains(";")]
    rank = rank.sort_values(["detected", "specific", "canonical", "mean"], ascending=False, kind="stable")
    rank = rank.drop_duplicates("symbol")
    rank["gene_id"] = rank.symbol.map(utils.symbol_steps(species)["to_ensembl"])
    rank = rank.dropna(subset=["gene_id"]).drop_duplicates("gene_id")
    return values.loc[rank.index].set_axis(rank.gene_id), rank.set_index("gene_id").symbol


def rna_genes(counts, features, value, samples, species):
    # log expression, then the most expressed row per Ensembl gene
    values = log_expression(counts, value).set_axis(samples.index, axis=1)
    features = features.reindex(values.index)
    symbol = utils.current_symbol(features.symbol, species, features.ensembl)
    prefix = "ENSMUSG" if species == "mouse" else "ENSG"
    ids = pd.Series(values.index, index=values.index)
    gene = features.ensembl.where(features.ensembl.astype(str).str.startswith(prefix))
    gene = gene.fillna(ids.where(ids.astype(str).str.startswith(prefix)))
    gene = gene.fillna(pd.Series(symbol, index=values.index).map(utils.symbol_steps(species)["to_ensembl"]))
    rows = pd.DataFrame(
        {"id": values.index, "gene_id": gene.to_numpy(), "mean": values.mean(axis=1).to_numpy(), "symbol": symbol}
    )
    rows = rows.dropna(subset=["gene_id"]).sort_values("mean", ascending=False).drop_duplicates("gene_id")
    return values.loc[rows.id].set_axis(rows.gene_id.to_numpy()), rows.set_index("gene_id").symbol


def observed(matrices, samples, spec):
    # genes seen in half of the donors of each layer
    # and, with min_side, in enough donors per age group
    keep = pd.Series(True, index=matrices[next(iter(matrices))].index)
    for values in matrices.values():
        keep &= values.notna().mean(axis=1) >= 0.5
        for side in ("early", "late") if spec.get("min_side") else ():
            donors = samples[spec["age"]["column"]].isin(spec["age"][side]).to_numpy()
            keep &= values.loc[:, donors].notna().sum(axis=1) >= spec["min_side"]
    return {name: values[keep.to_numpy()] for name, values in matrices.items()}


def pair_donors(selected, rules):
    # donors in both layers whose listed fields agree
    # (equal text, or within a tolerance)
    both = {name: frame.set_index("donor_id", drop=False) for name, frame in selected.items()}
    donors = both["protein"].index.intersection(both["rna"].index)
    keep = pd.Series(True, index=donors)
    for rule in rules.split():
        column, _, limit = rule.partition("=")
        protein, rna = both["protein"].loc[donors, column], both["rna"].loc[donors, column]
        if limit == "equal":
            keep &= (protein.astype(str).str.lower() == rna.astype(str).str.lower()) | protein.isna() | rna.isna()
        else:
            keep &= (protein - rna).abs() <= float(limit)
    return {name: frame.loc[sorted(donors[keep.to_numpy()])] for name, frame in both.items()}


def code_age(samples, age):
    # age per ten years, or late versus early groups
    # plus middle groups and contrasts within a side
    values = samples[age["column"]]
    if "per" in age:
        samples["age"] = values / age["per"]
        return
    samples["age"] = values.isin(age["late"]).astype(float)
    for group in age.get("middle", []):
        samples[f"middle_{group:g}"] = (values == group).astype(float)
    for side in ("early", "late"):
        bases = null_space(np.ones((1, len(age[side])))).T if len(age[side]) > 1 else []
        for number, basis in enumerate(bases, start=1):
            samples[f"{side}_{number}"] = values.map(dict(zip(age[side], basis, strict=True))).fillna(0.0)


def add_cell_indices(samples, matrices, labels, source, wanted, minimum):
    # cell index: mean z-score of its marker genes
    values, upper = matrices[source], labels[source].str.upper()
    for cell, markers in MARKERS[source].items():
        if wanted and cell not in wanted.split():
            continue
        rows = values.loc[upper.index[upper.isin(markers.split())].intersection(values.index)]
        if source == "protein":  # markers quantified in 80% of the samples
            rows = rows[rows.notna().mean(axis=1) >= 0.8]
            rows = rows[rows.std(axis=1) > 1e-8]
        else:
            rows = rows.dropna(how="any")
        if len(rows) >= minimum:
            z_scores = rows.sub(rows.mean(axis=1), axis=0).div(rows.std(axis=1), axis=0)
            samples[f"cell_{cell}"] = z_scores.mean(axis=0).reindex(samples.index).to_numpy()


def read_layers(spec):
    # matrices, features and selected samples per layer
    folder = utils.PREPARED / spec["data"]
    layers, selected = {}, {}
    for name in [name for name in ("protein", "rna") if name in spec]:
        stem = name if name == "protein" else "rna" if (folder / "rna.parquet").exists() else "pseudobulk"
        table = utils.read(folder / f"{stem}_samples.tsv", index_col="sample_id", dtype={"sample_id": str})
        table["sample_id"], table["donor_id"] = table.index.astype(str), table.donor_id.astype(str)
        features = utils.read(folder / f"{stem}_features.tsv", index_col=0, dtype={0: str})
        layers[name] = utils.read(folder / f"{stem}.parquet"), features
        selected[name] = table.query(spec[name], engine="python")
    assert all(len(table) > 0 for table in selected.values())
    if spec.get("pool"):  # counts of a donor's selected samples summed into one sample
        chosen, (matrix, features) = selected["rna"], layers["rna"]
        counts = matrix[chosen.sample_id].T.groupby(chosen.donor_id.to_numpy()).sum().T
        selected["rna"] = chosen.drop_duplicates("donor_id").assign(sample_id=lambda frame: frame.donor_id)
        layers["rna"] = counts, features
    if "donors_of" in spec:  # the donors of another unit
        donors = utils.read(utils.PREPARED / spec["donors_of"] / "samples.tsv", dtype=str).donor_id
        selected["rna"] = selected["rna"][selected["rna"].donor_id.isin(donors)]
    return layers, selected


def donor_table(spec, paired):
    # one row per donor: covariates, age coding, blocks
    first = next(iter(paired))
    samples = pd.DataFrame(index=pd.Index(paired[first].donor_id, name="sample_id"))
    samples["donor_id"] = paired[first].donor_id.to_numpy()
    for token in spec.get("cols", "").split():
        name, _, source = token.partition("=")
        transform, _, source = (source or name).rpartition(":")
        layer, _, column = source.rpartition(".")
        values = paired[layer or first][column]
        samples[name] = np.log1p(values.astype(float)).to_numpy() if transform else values.to_numpy()
    code_age(samples, spec["age"])
    blocks = spec.get("blocks", "").split()
    samples["permutation_block"] = "single"
    if blocks:
        samples["permutation_block"] = samples[blocks].fillna("single").astype(str).agg("|".join, axis=1)
    return samples


units, excluded = [], []
for unit, spec in UNITS.items():
    species = "mouse" if spec["data"].startswith(MOUSE) else "human"
    layers, selected = read_layers(spec)
    # donors of both layers; groups need enough donors
    paired = pair_donors(selected, spec["pair"]) if "pair" in spec else dict(selected)
    if "groups" in spec:  # genotype groups need enough donors at every age
        by, column, minimum = spec["groups"]
        counts = paired["protein"].groupby([by, column]).size().unstack(fill_value=0)
        enough = counts.index[(counts >= minimum).all(axis=1)]
        paired = {name: frame[frame[by].isin(enough)] for name, frame in paired.items()}
    # one row per donor
    samples = donor_table(spec, paired)
    # gene-level matrices of the paired donors
    matrices, labels = {}, {}
    for name, (matrix, features) in layers.items():
        if name == "protein":
            chosen, kept = selected[name].sample_id.to_numpy(), paired[name].sample_id.to_numpy()
            matrices[name], labels[name] = protein_genes(matrix, features, chosen, kept, samples, species)
        else:
            counts = matrix[paired[name].sample_id.to_numpy()]
            matrices[name], labels[name] = rna_genes(counts, features, spec.get("value"), samples, species)
    # genes in every layer and observed in enough donors
    common = matrices[next(iter(paired))].index
    for values in matrices.values():
        common = common.intersection(values.index)
    matrices = observed({name: values.loc[common] for name, values in matrices.items()}, samples, spec)
    # cell-type indices, then complete cases
    if "cells" in spec:
        add_cell_indices(samples, matrices, labels, *spec["cells"])
    complete = samples[spec.get("complete", "").split()].notna().all(axis=1).to_numpy()
    if not complete.all():  # donors without a covariate leave and the gene filter is repeated
        samples = samples[complete]
        matrices = observed({name: values.loc[:, complete] for name, values in matrices.items()}, samples, spec)
    if len(samples) < MIN_DONORS:
        excluded.append(unit)
        continue
    assert samples.index.is_unique and all(values.columns.equals(samples.index) for values in matrices.values())
    assert all(len(values) > 2000 for values in matrices.values())
    for layer, values in matrices.items():
        utils.write(values, utils.PREPARED / unit / f"{layer}.parquet")
    utils.write(samples, utils.PREPARED / unit / "samples.tsv")
    # units.tsv lists the covariates that exist in samples.tsv
    wanted = [spec.get("numeric", "").split(), spec.get("categorical", "").split()]
    present = [[name for name in group if name in samples] for group in wanted]
    for name in set(sum(wanted, [])) - set(sum(present, [])):
        print(unit, "covariate not in the samples:", name)
    row = dict(unit=unit, dataset=spec["data"], numeric=",".join(present[0]), categorical=",".join(present[1]))
    units.append(row | dict(blocks=",".join(spec.get("blocks", "").split()), hidden_factors=spec.get("hidden", 0)))
    print(
        f"{unit}: {len(samples)} samples, "
        + ", ".join(f"{len(values)} {name} genes" for name, values in matrices.items())
    )
assert len(units) + len(excluded) == len(UNITS)
print(f"fewer than {MIN_DONORS} donors, not written: {', '.join(excluded) or 'none'}")
utils.write(pd.DataFrame(units).set_index("unit"), utils.PREPARED / "units.tsv")
