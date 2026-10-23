# Age effect of genes against the load of a modification (log of 1 + sites on the canonical transcript; human m6A: the
# peaks of the human brain m6A-seq map), per SD of the load over the genes of the fit: OLS with HC3 per unit and layer,
# contexts combined by random effects; in the
# mRNA-adjusted layer also on the genes with a protein lifetime (with and without it as covariate) and on the others
# reads scripts 09, 10 and 12 and the unit levels; writes output/21_modification_load/
import numpy as np
import pandas as pd
import statsmodels.api as sm
from joblib import Parallel, delayed

import model
import utils

OUT = utils.OUT / "21_modification_load"
LAYERS = ["protein_only", "mrna_only", "mrna_adjusted"]
COVARIATES = utils.LOAD_COVARIATES
LIFETIME = "log2_protein_lifetime"
# model: covariates, and the columns a gene needs to enter the fit; genes with a lifetime are fitted with and without it
MODELS = {
    "without protein lifetime": (COVARIATES, COVARIATES),
    "lifetime genes without lifetime": (COVARIATES, [LIFETIME, *COVARIATES]),
    "with protein lifetime": ([LIFETIME, *COVARIATES], [LIFETIME, *COVARIATES]),
    "genes without a lifetime": (COVARIATES, COVARIATES),
}
LENGTHS = ["log_tx_len", "log_gene_len", "log1p_intron_len", "log1p_utr3_len"]
LADDER = {
    "none": [],
    "lengths": LENGTHS,
    "expression": [*LENGTHS, "level"],
    "neuron share": [*LENGTHS, "level", "ct_frac_neuron"],
    "programmes": COVARIATES,
    "other loads": [*COVARIATES, "load_Psi", "load_m5C", "load_m7G"],
}
GLIA_SHARE = 0.15  # glia-enriched genes: neuron share below 0.15, about the lower quartile of genes
UBIQUITOUS = 0.1  # ubiquitous genes: at least this share in each of neurons, astrocytes, oligodendrocytes and microglia
CHECKS = {  # name: predictor, covariates added to the primary ones, genes kept; the unmarked motifs are the total minus
    # the marked ones, so as a covariate they fall with the load by construction
    "DRACH as covariate": ("load_m6A", ["drach"], None),
    "unmarked DRACH as covariate": ("load_m6A", ["drach_unmarked"], None),
    "unmarked DRACH alone": ("drach_unmarked", [], None),
    "unmarked DRACH given load": ("drach_unmarked", ["load_m6A"], None),
    "GLORI stoichiometry": ("glori_stoichiometry", [], None),
    "glia-enriched genes": ("load_m6A", [], lambda d: d.ct_frac_neuron < GLIA_SHARE),
    "ubiquitous genes": ("load_m6A", [], lambda d: d.filter(like="ct_frac_").min(axis=1) >= UBIQUITOUS),
}
SITE_BINS, SITE_CLASSES = [-1, 0, 1, 3, np.inf], ["0", "1", "2-3", ">=4"]
WORKERS = 40

features = {
    species: utils.read(utils.OUT / "09_transcript_features" / f"{species}_features.tsv.gz").set_index("gene_id")
    for species in ("mouse", "human")
}
effects = {
    layer: utils.read(
        utils.OUT / "10_protein_age_effects" / f"{layer}.tsv.gz", usecols=["unit", "gene_id", "estimate", "status"]
    ).query("status == 'ok'")
    for layer in LAYERS
}
cohort_of = utils.read(utils.PREPARED / "units.tsv", index_col="unit").dataset


def per_sd(v):
    return v / v.std()


def draws_fit(data, covariates, Y):
    # OLS of every column of Y on the m6A load and the covariates: load coefficient and its HC3 standard error
    X = np.column_stack([np.ones(len(data)), data.load_m6A, data[covariates]]).astype(float)
    inverse = np.linalg.pinv(X)
    beta = inverse @ Y
    residual = (Y - X @ beta) / (1 - np.einsum("ij,ji->i", X, inverse))[:, None]
    return beta[1], np.sqrt((inverse[1][:, None] ** 2 * residual**2).sum(axis=0))


def pool(est, se):
    # DerSimonian-Laird estimate of every column of cohorts x draws, as utils.random_effects (one cohort: tau2 = 0)
    w = 1 / se**2
    q = (w * (est - (w * est).sum(axis=0) / w.sum(axis=0)) ** 2).sum(axis=0)
    tau2 = 0
    if len(est) > 1:
        tau2 = np.maximum(0, (q - (len(est) - 1)) / (w.sum(axis=0) - (w**2).sum(axis=0) / w.sum(axis=0)))
    w = 1 / (se**2 + tau2)
    return (w * est).sum(axis=0) / w.sum(axis=0)


def unit_fits(unit, layer, table):
    # one regression per modification and covariate set (the lifetime check only in the mRNA-adjusted layer), and the
    # mean age effect by site count
    species = "mouse" if table.gene_id.iloc[0].startswith("ENSMUSG") else "human"
    protein, rna, samples, _ = utils.load_unit(unit)
    level = (protein if layer == "protein_only" else rna).mean(axis=1).rename("level")
    genes = table.set_index("gene_id")[["estimate"]].join(features[species]).join(level)
    key = dict(unit=unit, layer=layer, species=species, n_samples=len(samples))
    counts = {mod: f"n_sites_{mod}" for mod in utils.MODIFICATIONS}  # human m6A: brain peaks (script 09)
    rows, classes, adjusted, ladder, checks = [], [], [], [], []
    for name, (covariates, required) in MODELS.items():
        if name != "without protein lifetime" and layer != "mrna_adjusted":
            continue
        data = genes.dropna(subset=["estimate", *required])
        if name == "genes without a lifetime":
            data = data[data[LIFETIME].isna()]
        for mod in utils.MODIFICATIONS:
            count = counts[mod]
            data = data.assign(load=per_sd(np.log1p(data[count])))
            fit = model.property_effect(data, "load", covariates)
            rows.append(key | dict(modification=mod, model=name, n_carriers=(data.load > 0).sum()) | fit)
            if name == "without protein lifetime":
                sites = pd.cut(data[count], SITE_BINS, labels=SITE_CLASSES)
                by_class = data.groupby(sites, observed=True).estimate.agg(n_genes="size", mean="mean", se="sem")
                classes.append(by_class.rename_axis("sites").reset_index().assign(**key, modification=mod))
                dummies = pd.get_dummies(sites, drop_first=True, dtype=float)  # classes against genes without sites
                X = sm.add_constant(pd.concat([dummies, data[covariates]], axis=1).astype(float))
                fit = sm.OLS(data.estimate.astype(float), X).fit(cov_type="HC3")
                effect = pd.DataFrame(dict(sites=dummies.columns, estimate=fit.params[dummies.columns].to_numpy()))
                effect["se"] = fit.bse[dummies.columns].to_numpy()
                adjusted.append(effect.assign(**key, modification=mod))
        if name == "without protein lifetime":  # m6A load with the covariates added in turn, on the same genes
            data = data.assign(**{f"load_{mod}": per_sd(np.log1p(data[counts[mod]])) for mod in utils.MODIFICATIONS})
            rho = data.load_m6A.corr(data.log_tx_len, method="spearman")
            for step, step_covariates in LADDER.items():
                fit = model.property_effect(data, "load_m6A", step_covariates)
                ladder.append(key | dict(step=step, rho_tx_len=rho) | fit)
            drach = dict(drach=np.log1p(data.n_drach), drach_unmarked=np.log1p(data.n_drach_unmarked))
            data = data.assign(**{k: per_sd(v) for k, v in drach.items()})
            if species == "mouse":
                data = data.assign(glori_stoichiometry=per_sd(np.log1p(data.stoich_m6A_glori)))
            X = sm.add_constant(data[COVARIATES].astype(float))
            residual = {c: sm.OLS(data[c].astype(float), X).fit().resid for c in ("load_m6A", "drach_unmarked")}
            partial = np.corrcoef(residual["load_m6A"], residual["drach_unmarked"])[0, 1]  # falls by construction
            for check, (predictor, extra, keep) in CHECKS.items():
                if predictor in data:
                    subset = data if keep is None else data[keep(data)]
                    fit = model.property_effect(subset, predictor, [*COVARIATES, *extra])
                    checks.append(key | dict(check=check, predictor=predictor, partial_r=partial) | fit)
            # gene length per SD in the same model: the gene property known to mark age-related loss (Stoeger 2022)
            others = ["load_m6A", *[c for c in COVARIATES if c != "log_gene_len"]]
            fit = model.property_effect(data.assign(gene_length=per_sd(data.log_gene_len)), "gene_length", others)
            checks.append(key | dict(check="gene length given load", predictor="gene_length", partial_r=partial) | fit)
            # the same coefficient in the age-permuted fits (genes with every draw)
            null = utils.read(utils.OUT / "12_permutation_fdr" / f"null_estimate_{unit}_{layer}.parquet")
            data = data.loc[data.index.intersection(null.index[null.notna().all(axis=1)])]
            observed = draws_fit(data, COVARIATES, data[["estimate"]].to_numpy(float))
            permuted = draws_fit(data, COVARIATES, null.loc[data.index].to_numpy(float))
            draws = key | dict(observed=observed, null=permuted)
    return rows, classes, adjusted, ladder, checks, draws


# every unit and layer in parallel
jobs = [delayed(unit_fits)(unit, layer, table) for layer in LAYERS for unit, table in effects[layer].groupby("unit")]
results = Parallel(n_jobs=WORKERS)(jobs)


FAMILIES = [*utils.CONTEXTS, ("mouse_protein_external", "protein_only", utils.EXTERNAL)]


def stack(i):
    return pd.DataFrame([row for result in results for row in result[i]]).assign(cohort=lambda d: d.unit.map(cohort_of))


def combine(table, by, medians=()):
    # contexts of utils.CONTEXTS and the external validation set: coefficients of the units of a cohort averaged,
    # cohorts combined by random effects (the z of a unit does not carry the variation between cohorts), per group of
    # the by columns
    rows = []
    for context, layer, units in FAMILIES:
        family = table[(table.layer == layer) & table.unit.isin(units.split())]
        for key, group in family.groupby(by, sort=False):
            cohorts = group.groupby("cohort")[["estimate", "se"]].mean()
            effect = utils.random_effects(cohorts.estimate, cohorts.se)
            row = dict(context=context, layer=layer) | dict(zip(by, key, strict=True))
            row |= {column: group[column].median() for column in medians}
            row |= dict(n_units=len(group), n_cohorts=len(cohorts), n_negative=(cohorts.estimate < 0).sum())
            row |= dict(estimate=effect["est"], effect_pct=utils.pct(effect["est"]))
            row |= dict(effect_pct_ci_low=utils.pct(effect["ci_low"]), effect_pct_ci_high=utils.pct(effect["ci_high"]))
            rows.append(row | dict(z=effect["z"], p=effect["p"], i2=effect["i2"]))
    return pd.DataFrame(rows)


coefficients = stack(0)
coefficients["q"] = coefficients.groupby(["unit", "layer", "model"]).p.transform(utils.bh)  # BH over modifications
utils.write(coefficients, OUT / "modification_load_effects.tsv")
# BH over the modifications of a context, layer and model
combined = combine(coefficients, ["model", "modification"])
combined["q"] = combined.groupby(["context", "layer", "model"]).p.transform(utils.bh)
utils.write(combined.drop(columns="estimate"), OUT / "combined_load_effects.tsv")

# ladder: retained_pct, the coefficient as a percentage of the coefficient without covariates; rho_tx_len, median over
# units of the Spearman correlation of load and transcript length
steps = combine(stack(3), ["step"], medians=["rho_tx_len"])
steps["retained_pct"] = 100 * steps.estimate / steps.groupby(["context", "layer"]).estimate.transform("first")
utils.write(steps.drop(columns="estimate"), OUT / "covariate_ladder.tsv")
# is it m6A? the m6A load against DRACH motifs, GLORI stoichiometry and gene subsets (CHECKS); partial_r: of load
# and unmarked motifs
checks = combine(stack(4), ["check", "predictor"], medians=["partial_r"])
utils.write(checks.drop(columns="estimate"), OUT / "load_checks.tsv")

# age-permutation null of the pooled m6A load coefficient (primary model): the permuted fits of script 12 for every
# unit pooled like the observed ones, draw by draw (permutation_null_draws); two-sided empirical P
draws = pd.DataFrame([result[5] for result in results]).assign(cohort=lambda d: d.unit.map(cohort_of))


def cohort_mean(cohorts, column, i):
    # coefficient (i = 0) or standard error (1) averaged over the units of each cohort
    return np.array([np.mean([fit[i] for fit in group[column]], axis=0) for _, group in cohorts])


null_rows, null_draws = [], []
for context, layer, units in FAMILIES:
    cohorts = draws[(draws.layer == layer) & draws.unit.isin(units.split())].groupby("cohort")
    estimate, se = cohort_mean(cohorts, "observed", 0)[:, 0], cohort_mean(cohorts, "observed", 1)[:, 0]
    observed = utils.random_effects(estimate, se)
    null = pool(cohort_mean(cohorts, "null", 0), cohort_mean(cohorts, "null", 1))
    row = dict(context=context, layer=layer, n_cohorts=cohorts.ngroups, n_draws=len(null))
    row |= dict(effect_pct=utils.pct(observed["est"]), null_pct_low=utils.pct(np.percentile(null, 2.5)))
    row |= dict(null_pct_high=utils.pct(np.percentile(null, 97.5)), p=observed["p"])
    null_rows.append(row | dict(p_permutation=utils.perm_pvalue(observed["est"], null)))
    draws_pct = dict(draw=range(1, len(null) + 1), effect_pct=utils.pct(null))
    null_draws.append(pd.DataFrame(dict(context=context, layer=layer) | draws_pct))
utils.write(pd.DataFrame(null_rows), OUT / "permutation_null.tsv")
utils.write(pd.concat(null_draws), OUT / "permutation_null_draws.tsv")

shown = ["context", "layer", "modification", "n_cohorts", "n_negative", "effect_pct", "z", "i2", "q"]
print(combined[combined.model == "without protein lifetime"][shown].round(4).to_string(index=False))
print(checks[["context", "layer", "check", "n_cohorts", "effect_pct", "p", "i2"]].round(4).to_string(index=False))
print(pd.DataFrame(null_rows).round(4).to_string(index=False))

# mean age effect by site class with a normal 95% interval (site_count_effects), and the classes against genes without
# sites adjusted for the covariates (site_class_effects; combined_site_class_effects: pooled as the load), in percent
classes = pd.concat([table for result in results for table in result[1]])
for column, width in [("effect_pct", 0), ("effect_pct_ci_low", -utils.Z95), ("effect_pct_ci_high", utils.Z95)]:
    classes[column] = utils.pct(classes["mean"] + width * classes.se)
utils.write(classes.drop(columns=["mean", "se"]), OUT / "site_count_effects.tsv")
adjusted = pd.concat([table for result in results for table in result[2]]).assign(
    effect_pct=lambda d: utils.pct(d.estimate)
)
utils.write(adjusted, OUT / "site_class_effects.tsv")
combined_classes = combine(adjusted.assign(cohort=lambda d: d.unit.map(cohort_of)), ["modification", "sites"])
utils.write(combined_classes.drop(columns="estimate"), OUT / "combined_site_class_effects.tsv")
print(combined_classes.query("modification == 'm6A'").round(4).to_string(index=False))
