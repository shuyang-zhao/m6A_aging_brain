# Source data of Figures 1-6 and Additional file 1 Figs S1-S9: one small table per panel with the plotted values and
# their statistics, and Tables S1-S7 (Additional file 2); per-gene values only for the units shown
# reads the outputs of scripts 10, 13, 20, 21, 28, 30, 33, 34, 37, 38, 50, 52, 53 and 54, the prepared units, the site
# atlas and the CLIP provenance
# writes source_data/fig<N><panel>.tsv, sfig<N><panel>.tsv and supplementary_table<N>.tsv, replacing the tables there
import re

import numpy as np
import pandas as pd
from scipy import stats as sps

import utils

OUT = utils.ROOT / "source_data"
MODS = {"m6A": "m6A", "Psi": "Ψ", "m5C": "m5C", "m7G": "m7G"}
LAYERS = {"protein_only": "Protein", "mrna_adjusted": "mRNA-adjusted protein", "mrna_only": "mRNA"}
STUDIES = {"gtex": "GTEx", "brainseq2": "BrainSeq II", "brainseq3": "BrainSeq III"}
STUDIES |= {"matlik2024": "Mätlik 2024", "colladotorres2019": "Collado-Torres 2019", "filipovic2026": "Vignane 2026"}
PERCENT = {"effect_pct": "effect", "effect_pct_ci_low": "ci_low", "effect_pct_ci_high": "ci_high"}
ACRONYMS = dict(
    dlpfc="DLPFC", dpfc="DPFC", smtg="SMTG", ipl="IPL", ca1="CA1", ba9="BA9", ba24="BA24", c1="C1", nac="NAc", rna="RNA"
)
TARGETS = {"mouse_protein": "Mouse protein", "human_mrna": "Human mRNA", "human_neurons": "Human neuronal mRNA"}
TARGETS |= {"mouse_mrna_adjusted": "Mouse mRNA-adjusted protein", "human_protein": "Human protein"}
PROTEIN = ["mouse_protein", "human_protein"]  # the targets of the main-figure panels
WEIGHTS = {  # weight sets of the clock (script 52)
    "model": "Model",
    "without_m6a": "Without m6A features",
    "length_only": "Length only",
    "m6a_load": "Linear m6A load",
    "DEMAGALHAES_AGING": "MSigDB DEMAGALHAES_AGING",
    "LU_AGING_BRAIN": "MSigDB LU_AGING_BRAIN",
    "observed": "Observed z of other cohorts",
}
EXTERNAL = utils.EXTERNAL.split()
CONTEXTS = [*utils.CONTEXTS, ("mouse_protein_external", "protein_only", utils.EXTERNAL)]
CONTEXT_UNITS = {(context, layer): listed.split() for context, layer, listed in CONTEXTS}
PANELS = {  # context and layer of the protein panels of Figure 3
    ("mouse_protein", "protein_only"): "Mouse discovery",
    ("mouse_protein_external", "protein_only"): "Mouse external",
    ("human_protein", "protein_only"): "Human protein",
}


def result(script, name):
    return utils.read(utils.OUT / script / name)


def label(unit):
    # koopmans2026_dpfc_l1_3 -> Koopmans 2026 DPFC L1-3
    study, *words = unit.split("_")
    name = STUDIES.get(study) or re.sub(r"^(.)(\D*)(\d{4})$", lambda m: m[1].upper() + m[2] + " " + m[3], study)
    text = " ".join(ACRONYMS.get(word, word) for word in words if word != "bankadj")
    return f"{name} {re.sub(r'\bl(\d) (\d)', r'L\1-\2', text)}".strip()


def cohort_name(dataset):
    # benjamin2022_brainseq3_caudate -> Benjamin 2022
    return label(dataset.split("_")[0].replace("gtexv10", "gtex"))


def named(table, **names):
    # human-readable names for the values of the given columns
    return table.assign(**{column: table[column].map(mapping) for column, mapping in names.items()})


def rho_interval(table):
    # 95% interval of a Spearman rho on the Fisher scale, the scale of the pooling in scripts 52 and 54
    z, half = np.arctanh(table.rho), utils.Z95 / np.sqrt(table.n - 3)
    return table.assign(ci_low=np.tanh(z - half), ci_high=np.tanh(z + half))


LAYER_OF_ARM = {"Output": "protein", "Translation": "mRNA-adjusted protein", "Stability": "mRNA"}


def save(name, table):
    # the network scripts name the layers of the coupling model by their reading (arm); the figures use the layers
    if "arm" in table:
        table = table.assign(arm=table.arm.map(LAYER_OF_ARM)).rename(columns={"arm": "layer"})
    table.to_csv(OUT / f"{name}.tsv", sep="\t", index=False, float_format="%.4g")
    print(name, table.shape)


OUT.mkdir(exist_ok=True)
for stale in OUT.glob("*.tsv"):
    stale.unlink()
cohort_of = utils.read(utils.PREPARED / "units.tsv", index_col="unit").dataset  # dataset of each unit

# Supplementary Table 1: every dataset considered (ref/datasets.tsv) with its set (discovery: mouse; validation: human;
# external: the mouse proteomes of utils.EXTERNAL; disease) or the reason it was not used; supplementary_table1_units
# and Figure 1b: the analysis units with dataset, layers, samples, ages, model terms and contexts
datasets = utils.read(utils.ROOT / "ref" / "datasets.tsv", keep_default_na=False)
save("supplementary_table1", datasets)
units = utils.read(utils.PREPARED / "units.tsv", keep_default_na=False)
sources = utils.read(utils.ROOT / "data_sources.tsv", keep_default_na=False)
doi = sources.assign(dataset=sources.id.str.split("/").str[0]).drop_duplicates("dataset").set_index("dataset").doi
contexts = {}
for context, layer, listed in CONTEXTS:
    for unit in listed.split():
        contexts.setdefault(unit, []).append(f"{context} ({LAYERS[layer]})")
rows = []
for unit in units.itertuples():
    folder = utils.PREPARED / unit.unit
    samples = utils.read(folder / "samples.tsv")
    species = "mouse" if "age_months" in samples else "human"
    ages = samples.age_months if species == "mouse" else samples.filter(["age_years", "age_midpoint"]).iloc[:, 0]
    layers = [name for name, file in (("protein", "protein"), ("RNA", "rna")) if (folder / f"{file}.parquet").exists()]
    # discovery: mouse; validation: human; external: mouse proteomes scored by the clock only
    kind = "external" if unit.unit in EXTERNAL else "discovery" if species == "mouse" else "validation"
    rows.append(
        dict(unit=label(unit.unit), unit_id=unit.unit, dataset=cohort_name(unit.dataset), dataset_id=unit.dataset)
        | dict(species=species, set=kind)
        | dict(layers=" and ".join(layers), samples=len(samples), donors=samples.donor_id.nunique())
        | dict(age_min=ages.min(), age_max=ages.max(), age_unit="months" if species == "mouse" else "years")
        | dict(age_term=utils.AGE_SCALE[species], covariates=",".join(filter(None, [unit.numeric, unit.categorical])))
        | dict(hidden_factors=unit.hidden_factors, permutation_blocks=unit.blocks)
        | dict(contexts="; ".join(contexts.get(unit.unit, [])), doi=doi[unit.dataset])
    )
table1 = pd.DataFrame(rows)
assert set(table1.dataset_id) <= set(datasets.dataset_id[datasets.set != "not used"])
save("supplementary_table1_units", table1)

# 1b: one row per set, dataset and molecular layer (protein, bulk RNA, single-nucleus RNA): units, donors of its largest
# unit and age range; mouse bulk and single-nucleus RNA make one group
NUCLEI = CONTEXT_UNITS["human_neurons", "mrna_only"] + CONTEXT_UNITS["mouse_neurons", "mrna_only"]
GROUPS = {
    ("human", "Protein"): "Human protein",
    ("human", "Bulk RNA"): "Human bulk RNA",
    ("human", "Single-nucleus RNA"): "Human single-nucleus RNA",
    ("mouse", "Protein"): "Mouse protein",
    ("mouse", "Bulk RNA"): "Mouse RNA",
    ("mouse", "Single-nucleus RNA"): "Mouse RNA",
}
long = table1.assign(layer=table1.layers.str.split(" and ")).explode("layer")
long["layer"] = np.where(
    long.layer == "protein", "Protein", np.where(long.unit_id.isin(NUCLEI), "Single-nucleus RNA", "Bulk RNA")
)
long["group"] = [GROUPS[key] for key in zip(long.species, long.layer, strict=True)]
long["set"] = np.where(
    long.unit_id.isin(EXTERNAL), "external", np.where(long.species == "mouse", "discovery", "validation")
)
long["dataset_id"] = long.unit_id.map(cohort_of)
fig1b = long.groupby(["set", "group", "layer", "dataset_id"], sort=False).agg(
    species=("species", "first"),
    n_units=("unit_id", "size"),
    donors=("donors", "max"),
    age_min=("age_min", "min"),
    age_max=("age_max", "max"),
    age_unit=("age_unit", "first"),
)
fig1b = fig1b.reset_index().assign(dataset=lambda d: d.dataset_id.map(cohort_name))
twice = fig1b.duplicated(["group", "dataset"], keep=False)  # Hahn 2023: bulk atlas and nuclei
fig1b["dataset"] += np.where(twice, fig1b.layer.map({"Bulk RNA": " bulk", "Single-nucleus RNA": " nuclei"}), "")
assert fig1b.dataset_id.nunique() == units.dataset.nunique() and not fig1b.duplicated(["group", "dataset"]).any()
save("fig1b", fig1b[["set", "group", "layer", "dataset", "n_units", "donors", "age_min", "age_max", "age_unit"]])


# Figure 2: the models of gene-level ageing (script 50)
def model(name):
    return result("50_transcript_ageing_model", name)


MODELS = {
    "tabpfn": "TabPFN",
    "lightgbm": "LightGBM",
    "elastic_net": "Elastic net",
    "shuffled_target": "Shuffled target",
}
REGIONS = {"all": "All sites", "5UTR": "5' UTR", "CDS": "CDS", "3UTR": "3' UTR"}
SHARED = {"shared": "Mouse protein model", "shared_shuffled": "Shuffled target"}
# a, held-out Spearman rho of every model on the five targets
performance = model("performance.tsv")
columns = ["target", "model", "n", "rho", "rho_ci_low", "rho_ci_high", "p"]
save("fig2a", named(performance, target=TARGETS, model=MODELS)[columns])

# b, held-out R2 lost when a feature group is removed (the m6A block against its permuted copies); the other targets
# in Fig. S3c
ablation = model("block_ablation.tsv")
columns = ["target", "block", "delta_r2", "null_low", "null_high", "p"]
save("fig2b", named(ablation[ablation.target.isin(PROTEIN)], target=TARGETS)[columns])

# c, in silico deletion of the m6A sites: mean change of the prediction of the marked genes, all sites and by region,
# against the permuted copies of the m6A block
deletion = model("deletion.tsv")
columns = ["target", "region", "n_genes", "effect", "null_low", "null_high", "p"]
save("fig2c", named(deletion[deletion.target.isin(PROTEIN)], target=TARGETS, region=REGIONS)[columns])
permuted = model("deletion_permuted.tsv")
save("fig2c_permuted", named(permuted[permuted.target.isin(PROTEIN)], target=TARGETS, region=REGIONS))

# d, the mouse protein model applied to human genes: rho against the human protein and mRNA targets
transfer = model("transfer.tsv")
columns = ["target", "model", "n", "rho", "rho_ci_low", "rho_ci_high", "p"]
save("fig2d", named(transfer, target=TARGETS, model=SHARED)[columns])

# Figure 3: m6A load marks proteins that decline with age (script 21; mouse discovery and external, human protein)
sites = result("21_modification_load", "site_count_effects.tsv").query(
    "modification == 'm6A' and layer == 'protein_only'"
)
sites = sites.assign(cohort=sites.unit.map(cohort_of))
rows = []
for key, panel in PANELS.items():
    # mean age effect of the genes by number of m6A sites: cohort means of its units, 95% interval across cohorts
    by_cohort = sites[sites.unit.isin(CONTEXT_UNITS[key])].groupby(["cohort", "sites"]).effect_pct.mean()
    for count, group in by_cohort.groupby("sites"):
        mean, width = group.mean(), sps.t.ppf(0.975, len(group) - 1) * group.sem()
        rows.append(
            dict(panel=panel, sites=count, n_cohorts=len(group), effect=mean, ci_low=mean - width, ci_high=mean + width)
        )
save("fig3a", pd.DataFrame(rows).replace({"sites": {">=4": "≥4"}}))

# b, m6A load per cohort (coefficient of the units of a cohort averaged, as in the pooling) and its random-effects
# estimate over cohorts, model without protein lifetime
loads = result("21_modification_load", "modification_load_effects.tsv")
loads = loads.query("modification == 'm6A' and layer == 'protein_only' and model == 'without protein lifetime'")
combined = result("21_modification_load", "combined_load_effects.tsv")
combined = combined.query("modification == 'm6A' and layer == 'protein_only' and model == 'without protein lifetime'")
rows = []
for (context, layer), panel in PANELS.items():
    per_cohort = loads[loads.unit.isin(CONTEXT_UNITS[context, layer])].groupby("cohort")[["estimate", "se"]].mean()
    for cohort, fit in per_cohort.iterrows():
        low, high = (utils.pct(fit.estimate + sign * utils.Z95 * fit.se) for sign in (-1, 1))
        rows.append(
            dict(panel=panel, cohort=cohort_name(cohort), kind="cohort", effect=utils.pct(fit.estimate))
            | dict(ci_low=low, ci_high=high)
        )
    total = combined[(combined.context == context) & (combined.layer == layer)].rename(columns=PERCENT).iloc[0]
    rows.append(
        dict(panel=panel, cohort="Random effects", kind="pooled")
        | total[["effect", "ci_low", "ci_high", "p", "i2"]].to_dict()
    )
save("fig3b", pd.DataFrame(rows))


def panel_of(table):
    # the protein panels of Figure 3 as the label of each row; the other contexts and layers are None
    return [PANELS.get(key) for key in zip(table.context, table.layer, strict=True)]


def core(table):
    # the rows of the discovery mouse and the human proteins
    table["panel"] = panel_of(table)
    return table[table.panel.isin(["Mouse discovery", "Human protein"])].rename(columns=PERCENT)


# c, core robustness in mouse and human protein: the covariates added in turn (the fifth step is the primary model),
# the load beside DRACH motifs as covariate and in ubiquitous genes, and its age-permutation null
ladder = core(result("21_modification_load", "covariate_ladder.tsv"))
ladder = ladder[ladder.step.isin(["lengths", "expression", "neuron share", "programmes"])]
save("fig3c_ladder", ladder[["panel", "step", "n_cohorts", "effect", "ci_low", "ci_high", "p", "retained_pct"]])
main = combined.assign(check="m6A load")
checks = core(pd.concat([main, result("21_modification_load", "load_checks.tsv")]))
checks = checks[checks.check.isin(["m6A load", "DRACH as covariate", "ubiquitous genes"])]
save("fig3c_checks", checks[["panel", "check", "n_cohorts", "effect", "ci_low", "ci_high", "p", "i2"]])
null = core(result("21_modification_load", "permutation_null_draws.tsv"))
save("fig3c_null", null[["panel", "draw", "effect"]])
observed = core(result("21_modification_load", "permutation_null.tsv"))
save(
    "fig3c_null_observed",
    observed[["panel", "n_draws", "effect", "null_pct_low", "null_pct_high", "p", "p_permutation"]],
)

# Figure 4: zero-shot clock (script 52): Spearman rho of the sample score with age; sets: held-out (cohorts scored with
# weights fitted without them) and external (mouse, scored with the weights of all cohorts)
cohorts, pooled = result("52_proteome_clock", "cohorts.tsv"), result("52_proteome_clock", "pooled.tsv")
scored = rho_interval(cohorts.query("weights == 'model'")).assign(kind="cohort")
scored["cohort"] = scored.cohort.map(cohort_name)
total = pooled.query("weights == 'model'").assign(cohort="Pooled", kind="pooled")
columns = ["set", "species", "cohort", "kind", "n", "rho", "ci_low", "ci_high", "p", "i2", "k"]
save("fig4a", pd.concat([scored, total])[columns].astype({"n": "Int64", "k": "Int64"}))

# b, score against age in four cohorts, with the rho of panel a
SHOWN = ["scifo2025_pxd058684", "chen2025_pxd045403", "wang2025_pxd047296", "koopmans2026_pxd058441"]
scores = result("52_proteome_clock", "scores.tsv.gz")
scores = scores[scores.cohort.isin(SHOWN)].merge(cohorts.query("weights == 'model'")[["cohort", "n", "rho", "p"]])
scores = scores.assign(age_unit=np.where(scores.species == "mouse", "months", "years"))
scores = scores.assign(order=scores.cohort.map(SHOWN.index)).sort_values("order", kind="stable")
scores = scores.assign(cohort=scores.cohort.map(cohort_name)).rename(columns={"model": "score"})
save("fig4b", scores[["cohort", "set", "species", "donor", "age", "age_unit", "score", "n", "rho", "p"]])

# c, pooled rho of the model and the control weights per set, the difference of each control from the model, and the
# pooled rho of random weights (1,000 draws) as the null; Fig. S4c repeats the draws
compared = result("52_proteome_clock", "compared.tsv").rename(
    columns={
        "control": "weights",
        "difference": "diff",
        "ci_low": "diff_ci_low",
        "ci_high": "diff_ci_high",
        "p": "diff_p",
    }
)
fig4c = pooled[pooled.weights != "without_m6a"].merge(
    compared[["set", "species", "weights", "diff", "diff_ci_low", "diff_ci_high", "diff_p"]], how="left"
)
columns = ["set", "species", "weights", "rho", "ci_low", "ci_high", "p", "i2", "k", "p_random"]
save("fig4c", named(fig4c, weights=WEIGHTS)[[*columns, "diff", "diff_ci_low", "diff_ci_high", "diff_p"]])
random_null = result("52_proteome_clock", "random_null.tsv")
save("fig4c_null", random_null)

# Figure 5: the clock in neurodegeneration (script 53): donor scores (the mean over the regions of a donor), and the
# coefficients of the disease on the score (SD of the score; OLS with age, sex and batch; AD studies combined by random
# effects); HD: per 100 CAG repeats
disease = result("53_disease_clock", "scores.tsv.gz")
disease["study"] = disease.study.map(label)
save(
    "fig5a",
    disease.query("disease_name == 'AD'")
    .query("study in ['Ping 2020', 'Merrihew 2023']")
    .assign(group=lambda d: np.where(d.disease == 1, "AD", "Control"))
    .rename(columns={"model": "score"})[["study", "donor_id", "group", "score"]],
)
hd = disease.query("disease_name == 'HD'")
save("fig5b", hd.rename(columns={"disease": "cag", "model": "score"})[["donor_id", "cag", "age", "score"]])


def disease_effects(outcome):
    tests = result("53_disease_clock", "tests.tsv").query("outcome == @outcome")
    tests = tests.assign(study=tests.study.map(lambda study: study if study == "AD combined" else label(study)))
    per_cag = tests.disease == "HD"
    scale = np.where(per_cag, 100, 1)
    tests = tests.assign(contrast=np.where(per_cag, "per 100 CAG repeats", "disease vs control"))
    return tests.assign(**{c: tests[c] * scale for c in ["est", "se", "ci_low", "ci_high"]})[
        ["study", "disease", "contrast", "n", "est", "se", "ci_low", "ci_high", "p"]
    ]


save("fig5c", disease_effects("model"))

# Figure 6: the reading layer (script 54); a, individual m6A load index against age (Spearman rho per cohort and pooled)
index = rho_interval(result("54_reader_titration", "index_cohorts.tsv")).assign(kind="cohort")
index["cohort"] = index.cohort.map(cohort_name)
index_total = result("54_reader_titration", "index_pooled.tsv").assign(cohort="Pooled", kind="pooled")
columns = ["species", "cohort", "kind", "n", "rho", "ci_low", "ci_high", "p", "i2", "k", "n_negative"]
save("fig6a", pd.concat([index, index_total])[columns].astype({"n": "Int64", "k": "Int64", "n_negative": "Int64"}))

# b, coefficient of the index on the reader level pooled over the units of a species, against the control RBPs
titration = result("54_reader_titration", "titration.tsv")
save(
    "fig6b",
    titration[["species", "regulator", "n_units", "est", "ci_low", "ci_high", "p", "q", "control_pct", "n_controls"]],
)
controls = result("54_reader_titration", "titration_controls.tsv")
save("fig6b_controls", controls[["species", "regulator", "n_units", "est", "ci_low", "ci_high", "p"]])

screen = result("30_regulator_screen", "regulator_effects.tsv")
net = result("33_regulator_network", "edges.tsv")


def bridge(name):
    return result("38_reader_bridges", name)


# c, FMR1 edge to its m6A-marked (s) and unmarked (u) targets of the two brain CLIP lists (mouse; human through
# orthologs) in every unit and arm
FMRP_LISTS = {"Darnell2011_brain": "Darnell 2011", "Maurin2018_brain": "Maurin 2018"}
FMRP_LISTS |= {f"{source}_mapped": name for source, name in FMRP_LISTS.items()}
fmr1 = net[(net.regulator == "FMR1") & net.source.isin(FMRP_LISTS)].assign(mice="all")
fmr1 = fmr1.assign(list=lambda d: d.source.map(FMRP_LISTS), unit=lambda d: d.unit.map(label))
save("fig6c_fmr1", fmr1[["list", "species", "arm", "unit", "mice", "b_s", "ci_low", "ci_high", "b_u", "control_p95"]])

# percentile of the YTHDF module and YTHDF3 load edge among the load edges of the control RBPs (pooled over units)
hits = bridge("yth_coupling_pooled.tsv").query("target_type == 'm6A_load' and regulator in ['YTHDF', 'YTHDF3']")
save("fig6c_ythdf", hits[["species", "arm", "regulator", "pctile_signed", "p_control_two_sided"]])

# Additional file 1 Fig. S1: resource and concordance. a, age t of the orthologs of one pair: Koopmans 2026 DPFC L1-3
# protein, Takasugi 2024 mRNA-adjusted protein; b, human protein age effects against mouse age effects on each layer
EFFECT_COLUMNS = ["unit", "gene_id", "symbol", "estimate", "t", "q", "effect_pct", "status"]
effects = {}
for layer in LAYERS:  # script 10: age effects of the genes each unit estimated
    table = utils.read(utils.OUT / "10_protein_age_effects" / f"{layer}.tsv.gz", usecols=EFFECT_COLUMNS)
    effects[layer] = table[table.status == "ok"]
concordance = result("13_species_concordance", "species_concordance.tsv").query("human_layer == 'protein_only'")
PAIR = ("koopmans2026_dpfc_l1_3", "takasugi2024")
protein = effects["protein_only"]
human_t = utils.by_symbol(protein[protein.unit == PAIR[0]]).t
orthologs = utils.orthologs().set_index("human").mouse
human_t = human_t[human_t.index.isin(orthologs.index)].rename(index=orthologs)
paired = {layer: utils.by_symbol(table[table.unit == PAIR[1]]) for layer, table in effects.items()}
paired = {layer: table[["estimate", "t", "q"]] for layer, table in paired.items()}
paired = pd.concat(paired, axis=1, join="inner")  # genes of all three mouse layers, as in script 13
pairs = pd.concat({"human_t": human_t, "mouse_t": paired["mrna_adjusted", "t"]}, axis=1, join="inner")
row = concordance[(concordance.human_unit == PAIR[0]) & (concordance.mouse_unit == PAIR[1])]
row = row[row.mouse_layer == "mrna_adjusted"].iloc[0]
assert len(pairs) == row.n and np.isclose(sps.spearmanr(pairs.human_t, pairs.mouse_t).statistic, row.rho)
save("sfig1a", pairs.rename_axis("gene").reset_index())
mouse = {"takasugi2024": "Mouse cerebrum (Takasugi 2024)"}
sfig1b = concordance.assign(
    human_unit=concordance.human_unit.map(label),
    mouse_cohort=concordance.mouse_unit.map(mouse),
    mouse_layer=concordance.mouse_layer.map(LAYERS),
)
assert not sfig1b.duplicated(["human_unit", "mouse_cohort", "mouse_layer"]).any()
save("sfig1b", sfig1b[["human_unit", "mouse_cohort", "mouse_layer", "rho", "ci_low", "ci_high", "p"]])

# Fig. S2: m6A maps; a, genes with an m6A site in each of the two mouse maps; b, c, the mouse m6A load and signature by
# map (script 28); d, quality of the maps with the human brain peaks: sites (peaks), share in a DRACH motif (peaks that
# contain one), share placed on a canonical transcript, sites per marked gene
CONTEXT_ROWS = {  # context and layer of the Additional file 1 panels: row label
    ("mouse_protein", "protein_only"): "Mouse protein",
    ("mouse_protein", "mrna_adjusted"): "Mouse mRNA-adjusted protein",
    ("mouse_mrna", "mrna_only"): "Mouse mRNA",
    ("mouse_neurons", "mrna_only"): "Mouse neuronal mRNA",
    ("human_protein", "protein_only"): "Human protein",
    ("human_striatum", "mrna_only"): "Human striatum mRNA",
    ("human_cortex", "mrna_only"): "Human cortex mRNA",
    ("human_hippocampus", "mrna_only"): "Human hippocampus mRNA",
    ("human_cerebellum", "mrna_only"): "Human cerebellum mRNA",
    ("human_neurons", "mrna_only"): "Human neuronal mRNA",
}


def rows_of(table):
    # the contexts and layers shown in the Additional file 1 panels, with their row label
    row = pd.Series(list(zip(table.context, table.layer, strict=True)), index=table.index).map(CONTEXT_ROWS)
    return table.assign(row=row)[row.notna()]


MAPS = {"m6a_hippocampus_cyto": "hippocampus cyto", "m6a_hippocampus_syn": "hippocampus syn"}
MAPS |= {utils.BRAIN_M6A: "human brain peaks"}
columns = ["map", "gene_id", "tx_offset", "ref_ok"]
atlas = utils.read(utils.RESOURCES / "site_atlas" / "sites.parquet", columns=columns)
atlas = atlas[atlas["map"].isin(MAPS)].replace({"map": MAPS})
placed = atlas[atlas.gene_id.notna() & atlas.tx_offset.notna() & atlas.ref_ok]
genes = placed.groupby("map").gene_id.apply(set).drop("human brain peaks")
membership = pd.DataFrame({name: pd.Series(True, index=sorted(ids)) for name, ids in genes.items()}).fillna(False)
save("sfig2a", membership.value_counts().rename("genes").reset_index())
for name, table, columns in (
    ("b", result("28_map_sensitivity", "combined_load_effects.tsv"), ["effect_pct", "ci_low_pct", "ci_high_pct"]),
    ("c", result("28_map_sensitivity", "combined_signatures.tsv"), ["z", "q"]),
):
    shown = rows_of(table)
    save(f"sfig2{name}", shown[["context", "layer", "definition", *columns]])
drach = result("28_map_sensitivity", "map_drach.tsv").set_index("map").drop(columns=["species", "n"])
quality = pd.DataFrame(dict(sites=atlas.groupby("map").size(), placed=placed.groupby("map").size()))
quality["sites_per_gene"] = quality.placed / placed.groupby("map").gene_id.nunique()
quality["placed"] /= quality.sites
quality = quality.join(drach.rename(index=MAPS)).rename(columns={"unit": "feature", "drach_fraction": "drach"})
save("sfig2d", quality.reset_index()[["map", "feature", "sites", "drach", "placed", "sites_per_gene"]])

# Fig. S3: details of the models of gene-level ageing; a, held-out rho of every model, and of the mouse protein model
# applied to human genes (transfer); b, per-gene effect of deleting all sites by the number of sites of the gene (mean,
# 95% CI); c, held-out R2 lost when a feature group is removed, on all five targets
transfer = named(
    transfer, target={"human_protein": "Human protein (mouse model)", "human_mrna": "Human mRNA (mouse model)"}
)
transfer = transfer.assign(model=transfer.model.map({"shared": "TabPFN", "shared_shuffled": "Shuffled target"}))
columns = ["target", "model", "n", "rho", "rho_ci_low", "rho_ci_high"]
save("sfig3a", pd.concat([named(performance, target=TARGETS, model=MODELS), transfer])[columns])
genes, rows = model("deletion_genes.tsv.gz"), []
for target in PROTEIN:
    inputs = utils.OUT / "50_transcript_ageing_model" / "inputs" / f"{target}.parquet"
    count = utils.read(inputs, columns=["m6A_all"]).m6A_all
    effect = genes[genes.target == target].set_index("gene_id").effect_all
    count = np.expm1(count.reindex(effect.index)).round()
    classes = pd.cut(count, [0, 1, 3, 7, np.inf], labels=["1", "2-3", "4-7", "≥8"])
    for bin_, values in effect.groupby(classes, observed=True):
        mean, half = values.mean(), utils.Z95 * values.std() / np.sqrt(len(values))
        rows.append(
            dict(target=TARGETS[target], sites=bin_, n=len(values), mean=mean, ci_low=mean - half, ci_high=mean + half)
        )
save("sfig3b", pd.DataFrame(rows))
save("sfig3c", named(ablation, target=TARGETS)[["target", "block", "delta_r2", "null_low", "null_high", "p"]])

# Fig. S4: details of the clock; a, rho of every weight set in every cohort; b, correlation of the scores of the model
# without the m6A features with the scores of the model; c, pooled rho of random weights per set against the model
cohorts = cohorts.assign(cohort=cohorts.cohort.map(cohort_name))
save("sfig4a", named(cohorts, weights=WEIGHTS)[["set", "species", "cohort", "weights", "n", "rho", "p"]])
save("sfig4b", cohorts.query("weights == 'without_m6a'")[["set", "species", "cohort", "n", "r_model"]])
observed = pooled.query("weights == 'model'")[["set", "species", "rho", "p_random"]].rename(
    columns={"rho": "model_rho"}
)
save("sfig4c", random_null.merge(observed))

# Fig. S5: the clock in disease; a, score of every sample by region, cases and controls of the AD, PD and ALS cohorts
# (ALS: one region, frontal cortex); b, coefficients of the m6A part of the score (model minus the model without the
# m6A features), SD of the score
regional = result("53_disease_clock", "sample_scores.tsv.gz").query("disease_name != 'HD'")
regional = regional.assign(
    study=regional.study.map(label),
    region=regional.region.fillna("Frontal cortex"),
    group=np.where(regional.disease == 1, "Case", "Control"),
)
regional = regional.drop(columns="disease").rename(columns={"disease_name": "disease", "model": "score"})
save("sfig5a", regional[["study", "disease", "region", "donor_id", "group", "score"]])
save("sfig5b", disease_effects("m6a_part"))

# Fig. S6: m6A against the other marks and the other marks beside m6A
signatures = rows_of(result("20_modification_signatures", "combined_signatures.tsv"))
save("sfig6a", signatures.assign(modification=signatures.modification.map(MODS))[["row", "modification", "z", "q"]])
loads = result("21_modification_load", "combined_load_effects.tsv").query("model == 'without protein lifetime'")
loads = rows_of(loads)
save("sfig6b", loads.assign(modification=loads.modification.map(MODS))[["row", "modification", "z", "q"]])

# Fig. S7: YTHDF3 and the YTHDF readers; a, YTHDF3 age effect, mouse protein per cohort and combined, human mRNA per
# context
HUMAN = ["human_protein", "human_striatum", "human_cortex", "human_hippocampus", "human_cerebellum", "human_neurons"]
regulator_cohorts = result("30_regulator_screen", "regulator_cohorts.tsv")
regulator_cohorts = regulator_cohorts.query("regulator == 'YTHDF3' and context == 'mouse_protein'")
regulator_cohorts = regulator_cohorts.query("layer == 'protein_only'")
cohorts_ythdf3 = pd.DataFrame(
    dict(
        panel="Mouse protein",
        label=regulator_cohorts.cohort.map(cohort_name),
        effect=utils.pct(regulator_cohorts.estimate),
        ci_low=utils.pct(regulator_cohorts.estimate - utils.Z95 * regulator_cohorts.se),
        ci_high=utils.pct(regulator_cohorts.estimate + utils.Z95 * regulator_cohorts.se),
    )
)
ythdf3 = screen.query("regulator == 'YTHDF3'").drop(columns=["ci_low", "ci_high"])
mouse_protein = (ythdf3.context == "mouse_protein") & (ythdf3.layer == "protein_only")
ythdf3 = ythdf3[mouse_protein | (ythdf3.context.isin(HUMAN[1:]) & (ythdf3.layer == "mrna_only"))]
mouse_row = ythdf3.context == "mouse_protein"
ythdf3 = ythdf3.assign(
    panel=np.where(mouse_row, "Mouse protein", "Human mRNA"),
    label=np.where(mouse_row, "Combined", ythdf3.context.str.removeprefix("human_")),
).rename(columns=PERCENT)
save("sfig7a", pd.concat([cohorts_ythdf3, ythdf3[["panel", "label", "effect", "ci_low", "ci_high", "q"]]]))

# b: pooled load edge of the YTHDF module and YTHDF3 against the control RBPs (random effects over units)
edges_ythdf = bridge("yth_coupling_edges.tsv").query("target_type == 'm6A_load'")
controls = (
    edges_ythdf[edges_ythdf.is_control]
    .groupby(["species", "arm", "regulator"])
    .apply(lambda d: utils.random_effects(d.b_s, d.se_s)["est"] if len(d) > 1 else d.b_s.iloc[0], include_groups=False)
)
save("sfig7b_controls", controls.rename("edge").reset_index())
save("sfig7b", hits[["species", "arm", "regulator", "estimate", "pctile_signed"]])

# c: change of the target score after loss of a YTHDF reader (log2 intervals converted to %)
loss = bridge("yth_perturbation.tsv").query("primary == True and predicted != 0 and effect_pct.notna()")
loss = loss.assign(ci_low_pct=utils.pct(loss.ci_low), ci_high_pct=utils.pct(loss.ci_high))
columns = ["dataset", "reader", "target_reader", "readout", "comparison", "effect_pct", "ci_low_pct", "ci_high_pct"]
save("sfig7c", loss[[*columns, "sign_matches"]])

# Fig. S8: FMRP family and PRRC2A; a, sign check against knockouts
check = result("34_perturbation_signs", "sign_check.tsv").query(
    "regulator in ['FMR1', 'PRRC2A'] and layer == 'protein_only'"
)
columns = ["perturbation", "regulator", "tissue", "readout", "effect", "ci_low", "ci_high", "q", "agrees"]
save("sfig8a", check[columns])

# b, median edge of the regulators and modules over the units of both species with their primary source, in units of
# the 95th percentile of the control RBPs, with the YTHDF load edges of script 38
NETWORK = ["YTHDF1", "YTHDF2", "YTHDF3", "YTHDC1", "YTHDC2", "FMR1", "FXR1", "FXR2", "PRRC2A", "HNRNPC", "RBMX"]
NETWORK += ["IGF2BP2"]
MODULES = {"YTHDF": "YTHDF module", "FMRP_family": "FMRP family", "writer_complex": "Writer complex"}
MODULES |= {"erasers": "Erasers"}
ARMS = ["Output", "Translation", "Stability"]


def primary_edges(edges, regulators):
    # the edges of the primary source of each regulator and species; FMR1: the brain CLIP of Darnell 2011 in both
    # species (human units through orthologs)
    darnell = edges.source.isin(["Darnell2011_brain", "Darnell2011_brain_mapped"])
    keep = np.where(edges.regulator == "FMR1", darnell, edges.primary.eq(True))
    return edges[keep & edges.regulator.isin(regulators)]


edges = primary_edges(net, [*NETWORK, *MODULES])
edges = edges.assign(
    value=edges.b_s / edges.control_p95,
    m6A_specific=edges.edge_class == "m6A-specific",
    binding_only=edges.edge_class == "Binding-only",
)
edges["shared"] = edges.m6A_specific & edges.shared_edge.fillna(False).astype(bool)
cells = (
    edges.groupby(["regulator", "family", "arm"])
    .agg(
        value=("value", "median"),
        m6A_specific=("m6A_specific", "any"),
        shared=("shared", "any"),
        binding_only=("binding_only", "any"),
    )
    .reset_index()
)
GROUP_OF = {"reader": "Readers", "module": "Modules", "writer": "Writers, erasers", "eraser": "Writers, erasers"}
cells["group"] = cells.family.map(GROUP_OF)
cells["row"] = cells.regulator.replace(MODULES)
yth_load = bridge("yth_coupling_pooled.tsv").query(
    "target_type == 'm6A_load' and regulator in ['YTHDF', 'YTHDF3'] and arm in @ARMS"
)
yth_load = yth_load.assign(
    row=yth_load.regulator.map({"YTHDF": "YTHDF module", "YTHDF3": "YTHDF3"}) + " load, " + yth_load.species,
    value=yth_load.estimate / yth_load.control_abs_p95,
    m6A_specific=yth_load.pctile_signed >= 95,
    shared=False,
    binding_only=False,
    group="YTHDF load edge",
)
columns = ["group", "row", "arm", "value", "m6A_specific", "shared", "binding_only"]
save("sfig8b", pd.concat([cells, yth_load])[columns])

# c, edges shared across units
SHARED_EDGES = [
    *NETWORK[5:11],
    "YTHDC1",
    "YTHDF1",
    "YTHDF2",
    "YTHDF3",
    "FMRP_family",
    "YTHDF",
    "writer_complex",
    "erasers",
]
shared = primary_edges(net[net.arm == "Output"], SHARED_EDGES)
shared = shared.assign(unit=shared.unit.map(label), regulator=shared.regulator.replace(MODULES))
save("sfig8c", shared[["regulator", "unit", "species", "posterior_mean", "lfsr"]])

# Fig. S9: the targets of IGF2BP2 do not follow its decline; a, IGF2BP2 mRNA and the target scores of its primary
# tiers in the single-nucleus cohorts and, when GTEx is not in utils.EXCLUDED, in GTEx (region x age)
panels = [bridge("igf2bp2_nuclei_combined.tsv").query("not is_control").assign(panel="Single nucleus")]
if cohort_of.index.str.startswith("gtex_").any():
    regions = bridge("igf2bp2_striatum_x_age.tsv")
    regions = regions[~regions.is_control & regions.source.isin(["IGF2BP2:primary", "IGF2BP2 mRNA"])]
    panels.insert(0, regions.assign(panel="GTEx"))
columns = ["panel", "tier", "comparison", "estimate_pct", "ci_low_pct", "ci_high_pct"]
save("sfig9a", pd.concat(panels)[columns])

# b: genes bound with m6A in each combination of the four neural binding sources (T1, at least two)
NEURAL = ["conway2016_H9_eCLIP", "degrauwe2016_MGG8_PARCLIP", "han2022_A172_eCLIP", "preitner2016_E14brain_CLIP"]
bound = result("37_reader_compendium", "bound_genes.tsv").query("reader == 'IGF2BP2' and T1 and source in @NEURAL")
membership = pd.DataFrame(
    {source: pd.Series(True, index=bound[bound.source == source].human_gene_id.dropna().unique()) for source in NEURAL}
).fillna(False)
save("sfig9b", membership.value_counts().rename("genes").reset_index())
coupling = bridge("igf2bp2_coupling.tsv")
columns = ["dataset", "population", "b_s", "ci_low", "ci_high", "pctile_abs"]
save("sfig9c", coupling.assign(dataset=coupling.dataset.map(cohort_name))[columns])
protein_score = bridge("igf2bp2_proteome_pooled.tsv").query("not is_control and source == 'IGF2BP2:primary'")
save("sfig9d", protein_score[["species", "tier", "comparison", "estimate_pct", "ci_low_pct", "ci_high_pct"]])
loss = bridge("igf2bp2_perturbation.tsv")
loss = loss.assign(ci_low_pct=utils.pct(loss.ci_low), ci_high_pct=utils.pct(loss.ci_high))
save("sfig9e", loss[["dataset", "tier", "comparison", "effect_pct", "ci_low_pct", "ci_high_pct", "q"]])
neuron_glia = bridge("igf2bp2_neuron_glia_combined.tsv").query("not is_control")
save("sfig9f", neuron_glia[["set", "comparison", "pctile_z"]])

# Supplementary Table 2: modification site maps with their source, resolution and numbers of sites and genes
columns = ["map", "modification", "species", "resolution", "source", "gene_id", "tx_offset", "ref_ok"]
atlas = utils.read(utils.RESOURCES / "site_atlas" / "sites.parquet", columns=columns)
first = {column: (column, "first") for column in ["modification", "species", "resolution", "source"]}
maps = atlas.groupby("map").agg(**first, sites=("map", "size"))
placed = atlas[atlas.gene_id.notna() & atlas.tx_offset.notna() & atlas.ref_ok]
maps["genes"] = placed.groupby("map").gene_id.nunique()
maps = maps.assign(modification=maps.modification.map(MODS))
save("supplementary_table2", maps.reset_index())

# Supplementary Table 3: CLIP libraries of the readers with the read counts of each processing step
clip = utils.read(utils.RESOURCES / "reader_clip/provenance.tsv")
clip = clip.assign(deduplicated_reads=clip.deduplicated_reads.astype("Int64"))  # libraries without UMIs have none
save("supplementary_table3", clip.drop(columns="tools"))

# Supplementary Table 4: modification load per context, layer and covariate set, cohorts combined (script 21)
load = result("21_modification_load", "combined_load_effects.tsv")
load = load.assign(layer=load.layer.map(LAYERS), modification=load.modification.map(MODS))
columns = ["context", "layer", "model", "modification", "n_units", "n_cohorts", "n_negative", "effect_pct"]
save("supplementary_table4", load[[*columns, "effect_pct_ci_low", "effect_pct_ci_high", "z", "p", "q"]])

# Supplementary Table 5: age effect of every regulator in every context and layer (script 30), with the aged cohorts
columns = ["context", "layer", "regulator", "n_cohorts", "effect_pct", "effect_pct_ci_low", "effect_pct_ci_high"]
save("supplementary_table5", screen.assign(layer=screen.layer.map(LAYERS))[[*columns, "p", "q", "cohorts"]])

# Supplementary Table 6: binding sources and perturbation datasets of the reader compendium (script 37) with the
# publication of each accession (PMID and DOI from PubMed through the GEO record; GSE171790, GSE266096 and GSE306850
# have no PubMed link in GEO and were matched by the paper's text, or its title and the GEO contributors; ENCODE: Van
# Nostrand 2020)
PAPERS = {  # accession: study, PMID, DOI
    "GSE104867": ("Li 2018", "29855337", "10.1186/s13059-018-1436-y"),
    "GSE106607": ("Shi 2018", "30401835", "10.1038/s41586-018-0666-1"),
    "GSE119801": ("Regué 2019", "30692269", "10.1128/MCB.00528-18"),
    "GSE133898": ("Heck 2020", "32169943", "10.1261/rna.073502.119"),
    "GSE146726": ("Okholm 2020", "33287884", "10.1186/s13073-020-00812-8"),
    "GSE148039": ("Lasman 2020", "32943573", "10.1101/gad.340695.120"),
    "GSE148536": ("Lasman 2020", "32943573", "10.1101/gad.340695.120"),
    "GSE153688": ("Yu 2021", "34643063", "10.1002/advs.202101329"),
    "GSE153689": ("Yu 2021", "34643063", "10.1002/advs.202101329"),
    "GSE158742": ("Dixit 2021", "33023892", "10.1158/2159-8290.CD-20-0331"),
    "GSE171790": ("Zhuang 2023", "36670199", "10.1038/s41380-023-01953-z"),
    "GSE188846": ("Ito-Kureha 2022", "35879451", "10.1038/s41590-022-01268-1"),
    "GSE192705": ("Alvarez-Dominguez 2022", "35036870", "10.1016/j.isci.2021.103680"),
    "GSE192792": ("Han 2022", "35992092", "10.1016/j.isci.2022.104823"),
    "GSE194207": ("Flamand 2022", "35438793", "10.1093/nar/gkac251"),
    "GSE199962": ("Niu 2022", "36039295", "10.1016/j.isci.2022.104908"),
    "GSE201070": ("Zou 2023", "36694229", "10.1186/s13059-023-02862-8"),
    "GSE214882": ("Zou 2023", "37949069", "10.1016/j.molcel.2023.10.028"),
    "GSE230846": ("Nabeel-Shah 2024", "39303720", "10.1016/j.molcel.2024.08.037"),
    "GSE264498": ("Zhang 2024", "39543097", "10.1038/s41467-024-54000-y"),
    "GSE266096": ("Zhang 2025", "39476834", "10.1016/j.stem.2024.10.002"),
    "GSE270417": ("Jang 2025", "39587835", "10.1002/advs.202410806"),
    "GSE292558": ("Schymik 2026", "41943844", "10.7150/ijbs.122142"),
    "GSE306850": ("Zhou 2026", "42697900", "10.1038/s41380-026-03840-9"),
    "GSE315406": ("Her 2026", "42025161", "10.1016/j.cels.2026.101588"),
    "GSE326611": ("Shi 2026", "41973927", "10.1073/pnas.2518250123"),
    "GSE49339": ("Wang 2014", "24284625", "10.1038/nature12730"),
    "GSE63591": ("Wang 2015", "26046440", "10.1016/j.cell.2015.05.014"),
    "GSE73845": ("Degrauwe 2016", "27184842", "10.1016/j.celrep.2016.04.086"),
    "GSE78507": ("Conway 2016", "27068461", "10.1016/j.celrep.2016.03.052"),
    "GSE83822": ("Preitner 2016", "27385015", "10.1242/dev.128348"),
    "GSE86214": ("Shi 2017", "28106072", "10.1038/cr.2017.15"),
    "GSE90639": ("Huang 2018", "29476152", "10.1038/s41556-018-0045-z"),
    "GSE90684": ("Huang 2018", "29476152", "10.1038/s41556-018-0045-z"),
    "ENCSR": ("Van Nostrand 2020", "32728246", "10.1038/s41586-020-2077-3"),
}
PAPERS["GSE147849"] = PAPERS["GSE148039"]
PAPERS["Yu 2021 Table S3"] = PAPERS["GSE153688"]
ACCESSION = r"GSE\d+|ENCSR|Yu 2021 Table S3"  # a GEO series, an ENCODE experiment or a table of the study


def publication(accession):
    first = re.search(ACCESSION, accession)[0]
    return pd.Series(PAPERS.get(first, ("", "", "")), index=["study", "PMID", "DOI"])


compendium = utils.OUT / "37_reader_compendium"
columns = ["reader", "source", "species", "cell", "tissue_class", "rank", "method", "accession", "T1", "gene_id"]
bound = utils.read(compendium / "bound_genes.tsv", usecols=[*columns, "note"])
binding = (
    bound.groupby(["reader", "source"], sort=False)
    .agg(
        species=("species", "first"),
        cell_type=("cell", "first"),
        tissue_class=("tissue_class", "first"),
        rank=("rank", "first"),
        method=("method", "first"),
        accession=("accession", "first"),
        n_genes=("gene_id", "nunique"),
        n_genes_with_m6A=("T1", "sum"),
        note=("note", "first"),
    )
    .reset_index()
    .rename(columns={"source": "dataset"})
    .assign(kind="binding")
)
columns = ["dataset", "reader", "accession", "species", "cell_type", "direction", "readout", "gene_id", "padj", "note"]
trial = utils.read(compendium / "perturbation.tsv", usecols=columns)
perturbation = (
    trial.groupby(["dataset", "reader", "readout"], sort=False)
    .agg(
        species=("species", "first"),
        cell_type=("cell_type", "first"),
        method=("direction", "first"),
        accession=("accession", "first"),
        n_genes=("gene_id", "nunique"),
        n_genes_padj_below_0_05=("padj", lambda p: (p < 0.05).sum()),
        note=("note", "first"),
    )
    .reset_index()
    .assign(kind="perturbation")
)
table6 = pd.concat([binding, perturbation], ignore_index=True)
table6 = pd.concat([table6, table6.accession.apply(publication)], axis=1)
columns = ["kind", "reader", "dataset", "species", "cell_type", "method", "readout", "tissue_class", "rank", "n_genes"]
columns += ["n_genes_with_m6A", "n_genes_padj_below_0_05", "accession", "study", "PMID", "DOI", "note"]
save("supplementary_table6", table6[columns])

# Supplementary Table 7: the 90 control RBPs of the network (ENCODE eCLIP proteins that do not read m6A, quantified in
# the proteome units) and the RBPs left out, with the reason
controls = result("33_regulator_network", "control_rbps.tsv").sort_values("selected", ascending=False, kind="stable")
save("supplementary_table7", controls[["rbp", "selected", "reason", "cell_lines", "n_units_quantified", "n_units"]])
