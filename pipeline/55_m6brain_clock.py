# m6Brain clock: elastic net of real age (months in mouse, years in human) on the proteins standardised within their
# unit (missing = 0), the samples of the training units stacked, whose penalty factors come from the m6Brain gene
# predictions. Penalty factor p_j = exp(-gamma * s_j), s = percentile rank (0-1) of the prior over the genes, p
# normalised to mean 1; feature j is divided by p_j and the coefficients back-transformed, which is exact for the L1
# term while the L2 term scales with p_j^2. m6brain_clock: s = |m6Brain weight| of script 52 fitted without the scored
# cohort; elastic_net: gamma = 0 (unguided); m6a_load_prior: s = the m6A load of the gene (m6A alone as the prior).
# l1_ratio {0.1, 0.5, 0.9}, alpha (30-value path down to 3e-3 of the alpha that zeroes all coefficients) and, for the
# priors, gamma {0.5, 1, 2, 4} minimise the mean MSE of an inner CV (leave-one-cohort-out among the training cohorts;
# 5-fold by donor within a cohort), then the model is refitted on all training samples; at_floor flags a chosen alpha
# at the end of its path.
# Scored as script 52: donor mean over the units of a cohort of the predictions (Pearson r: unit-standardised scores,
# MAE: the raw predictions, which keep the offset between the cohort mean age and the training mean that within-unit
# standardisation cannot recover); r pooled on the Fisher scale; the clocks compared with the unguided net per cohort
# (Meng's SE) and pooled. Held-out: each cohort of mouse_protein / human_protein with the model of the other cohorts;
# external: utils.EXTERNAL with the mouse model of all cohorts. Within-cohort nested CV (outer 10-fold by donor x 10
# seeds, inner 5-fold by donor) in Luo 2024 (human) and Lu 2025 (mouse): MAE and r over the seeds, and the paired
# Wilcoxon test of the donor errors (averaged over the seeds) of each variant against the unguided net; the unguided
# net chose the end of the path in 11% of the Luo folds, and on a path down to 1e-3 (36 values) its MAE was the same
# (6.02 against 5.99 years).
# Every published clock (ref/published_clocks.tsv: proteomic, transcriptomic, single-cell, senescence and half-life
# weights trained elsewhere) is scored on the same cohorts by utils.clock_scores: genes outside the model table of the
# held-out cohort are unmeasured, mouse clocks reach human genes and human clocks mouse genes through one-to-one
# orthologs, half-lives are centred, and the input is the one of the clock's own preprocessing (utils.CLOCK_INPUT).
# native_modality.tsv: the transcriptomic clocks on the RNA and the protein models on the protein of the same samples
# (Takasugi 2024, Luo 2024), and the m6Brain clock against each on the shared donors (Meng's SE). The null of the MAE
# predicts the mean age of the training samples (cohorts.tsv: mae_null; within_cohort.tsv: the same folds and seeds);
# failures.tsv lists the cohorts where the m6Brain clock or its label-free mode has an r CI that includes 0, an MAE
# not below the null or a published clock of higher r. The benchmark adds these rows to the rows of script 52 and
# assigns every weight set a tier; benchmark_pooled.tsv ranks the pooled r within each set and species, overall and
# within the tier. theory.tsv (Luo only, last stage): the training donors of each fold subsampled to 25, 50, 100 and
# all, the unguided net, the m6Brain-guided net and the same prior permuted across the genes, with l1_ratio and gamma
# fixed at the values chosen most often in the nested run (unguided l1_ratio 0.9; guided and permuted l1_ratio 0.9,
# gamma 2) and alpha tuned by inner CV; MAE and r over 5 seeds, non-zero proteins, paired Wilcoxon tests.
# Run with OMP_NUM_THREADS=1 PYTHONWARNINGS=ignore: the lowest alphas of a path stop at max_iter, and the coordinate
# descent otherwise spins up every core. The first run takes hours; cache/ keeps the finished inner-CV cells and folds.
# reads output/52_proteome_clock/ (runs, inputs, benchmark.tsv, cohorts.tsv), output/prepared/ref/published_clocks.tsv
# (script 06) and the prepared units; writes output/55_m6brain_clock/ (weights/{target}__{cohort|all}__{clock}.tsv:
# gene_id and weight of the non-zero genes and the "(intercept)" row, score = intercept + standardised proteins @ weight;
# scores.tsv.gz (predicted age of every donor: donor mean of the raw predictions over the units of its cohort), cohorts,
# pooled, compared, benchmark, benchmark_pooled, native_modality, failures, theory, hyperparameters, within_cohort and
# within_cohort_seeds)
import joblib
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import stats as sps
from sklearn.linear_model import ElasticNet, enet_path

import utils

S52 = utils.OUT / "52_proteome_clock"
OUT = utils.OUT / "55_m6brain_clock"
TARGETS = {"mouse": "mouse_protein", "human": "human_protein"}
VARIANTS = ["elastic_net", "m6brain_clock", "m6a_load_prior"]
GROUPS = {"m6brain_clock": "m6Brain clock", "elastic_net": "reference", "m6a_load_prior": "reference"}
L1_RATIO, GAMMA, PATH = [0.1, 0.5, 0.9], [0.5, 1, 2, 4], (3e-3, 30)  # path: floor x alpha_max, number of alphas
WITHIN = {"luo2024_syn32136022": "human_protein", "lu2025_pxd053154": "mouse_protein"}
N_OUTER, N_INNER, N_SEEDS, WORKERS = 10, 5, 10, 48
# benchmark tiers: of the published clocks by source (the others are proteomic clocks), of the other weight sets by name
SOURCE_TIERS = dict.fromkeys(utils.CLOCK_INPUT, "RNA clocks applied to protein") | dict(
    palmer2021="signatures", saul2022="signatures", fornasiero2018="label-free"
)
WEIGHT_TIERS = dict.fromkeys(["DEMAGALHAES_AGING", "LU_AGING_BRAIN"], "signatures") | dict.fromkeys(
    ["model", "without_m6a", "length_only", "m6a_load", "random"], "label-free"
)
# units holding protein and RNA of the same samples
PAIRED = {"takasugi2024": ("mouse", "mouse_protein"), "luo2024_bankadj": ("human", "human_protein")}

cohort_of = utils.read(utils.PREPARED / "units.tsv", index_col="unit").dataset
clock_w = utils.clock_weights()
tiers = {c: SOURCE_TIERS.get(source, "protein clocks") for c, (source, _, _) in clock_w.items()} | WEIGHT_TIERS


def units_of(target):
    return next(units.split() for name, layer, units in utils.CONTEXTS if (name, layer) == (target, "protein_only"))


def standardised(x):
    return (x - x.mean()) / x.std()


def meng_se(r1, r2, r12, n):
    # SE of atanh r1 - atanh r2 for two correlations with a shared variable (Meng, Rosenthal & Rubin 1992)
    mean_sq = (r1**2 + r2**2) / 2
    f = np.minimum((1 - r12) / (2 * (1 - mean_sq)), 1)
    h = (1 - f * mean_sq) / (1 - mean_sq)
    return np.sqrt(2 * (1 - r12) * h / (n - 3))


def cached(name, compute):
    path = OUT / "cache" / f"{name}.pkl"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        joblib.dump(compute(), path)
    return joblib.load(path)


UNITS = {}


def unit_data(unit):
    # samples x all proteins (log abundance; and standardised within the unit, NaN kept), age in real units, donor
    if unit not in UNITS:
        protein, _, samples, _ = utils.load_unit(unit)
        samples = samples.loc[protein.columns]
        age = samples.age_months if "age_months" in samples else samples.age_years
        UNITS[unit] = protein.T, protein.T.pipe(standardised), age.to_numpy(), samples.donor_id.astype(str).to_numpy()
    return UNITS[unit]


def stack(units, genes):
    parts = [unit_data(u) for u in units]
    X = np.vstack([z.reindex(columns=genes).fillna(0).to_numpy() for _, z, _, _ in parts])
    cohort = np.concatenate([[cohort_of[u]] * len(age) for u, (_, _, age, _) in zip(units, parts, strict=True)])
    return X, np.concatenate([a for _, _, a, _ in parts]), np.concatenate([d for _, _, _, d in parts]), cohort


def rank01(x):
    r = pd.Series(np.asarray(x)).rank()
    return ((r - 1) / (len(r) - 1)).to_numpy()


def penalty(s, gamma):
    p = np.exp(-gamma * s)
    return p / p.mean()


def alpha_grid(Xs, y, l1):
    floor, n = PATH
    Xc = Xs - Xs.mean(0)
    return np.abs(Xc.T @ (y - y.mean())).max() / (len(y) * l1) * np.logspace(0, np.log10(floor), n)


def path_mse(X, y, folds, p, l1):
    # validation MSE along the alpha path, mean over the inner folds
    Xs = X / p
    alphas = alpha_grid(Xs, y, l1)
    mse = []
    for tr, va in folds:
        mu, ybar = Xs[tr].mean(0), y[tr].mean()
        coefs = enet_path(Xs[tr] - mu, y[tr] - ybar, l1_ratio=l1, alphas=alphas)[1]
        mse.append((((Xs[va] - mu) @ coefs + ybar - y[va, None]) ** 2).mean(0))
    return np.mean(mse, 0)


def cells(variant):
    return [(g, l1) for g in ([0] if variant == "elastic_net" else GAMMA) for l1 in L1_RATIO]


def refit(X, y, s, mse):
    # best (gamma, l1_ratio, alpha) of the inner CV, refitted on all training samples
    err, (gamma, l1), i = min((m.min(), c, m.argmin()) for c, m in mse.items())
    p = penalty(s, gamma)
    Xs = X / p
    alpha = alpha_grid(Xs, y, l1)[i]
    fit = ElasticNet(alpha=alpha, l1_ratio=l1, max_iter=5000).fit(Xs, y)
    w = fit.coef_ / p
    floor, n = PATH
    hyper = dict(l1_ratio=l1, alpha=alpha, alpha_rel=floor ** (i / (n - 1)), at_floor=i == n - 1, gamma=gamma)
    return w, fit.intercept_, hyper | dict(n_nonzero=(w != 0).sum(), inner_mse=err)


def priors_of(target, k, genes):
    runs = utils.read(S52 / "runs" / f"{target}__{k}.parquet")
    m6a = utils.read(S52 / "inputs" / f"{target}__{k}.parquet").m6A_all
    priors = [np.zeros(len(genes)), rank01(runs.model.abs()), rank01(m6a.reindex(genes))]
    return dict(zip(VARIANTS, priors, strict=True))


def loco_data(target, k):
    # training samples of the held-out cohort k (or of all cohorts), inner folds = the training cohorts
    genes = utils.read(S52 / "runs" / f"{target}__{k}.parquet").index
    X, y, _, cohort = stack([u for u in units_of(target) if k == "all" or cohort_of[u] != k], genes)
    folds = [(np.flatnonzero(cohort != c), np.flatnonzero(cohort == c)) for c in np.unique(cohort)]
    return genes, X, y, folds, priors_of(target, k, genes)


def loco_cell(target, k, variant, gamma, l1):
    def compute():
        _, X, y, folds, priors = loco_data(target, k)
        return path_mse(X, y, folds, penalty(priors[variant], gamma), l1)

    return cached(f"{target}__{k}__{variant}__{gamma}__{l1}", compute)


def loco_fit(target, k, variant, mse):
    genes, X, y, _, priors = loco_data(target, k)
    w, b, hyper = refit(X, y, priors[variant], mse)
    keep = np.flatnonzero(w)
    weights = pd.DataFrame(dict(gene_id=[*genes[keep], "(intercept)"], weight=[*w[keep], b]))
    utils.write(weights, OUT / "weights" / f"{target}__{k}__{variant}.tsv")
    return genes, w, b, dict(stage="loco", target=target, cohort=k, variant=variant) | hyper


def group_folds(groups, k, seed):
    ids = np.random.default_rng(seed).permutation(np.unique(groups))
    fold = pd.Series(np.arange(len(ids)) % k, index=ids)[groups].to_numpy()
    return [(np.flatnonzero(fold != i), np.flatnonzero(fold == i)) for i in range(k)]


def outer_fold(cohort, X, y, donor, tr, te, s, variant, seed, fold):
    def compute():
        inner = group_folds(donor[tr], N_INNER, 1000 * seed + fold)
        mse = {c: path_mse(X[tr], y[tr], inner, penalty(s, c[0]), c[1]) for c in cells(variant)}
        w, b, hyper = refit(X[tr], y[tr], s, mse)
        return X[te] @ w + b, dict(stage="nested", seed=seed, fold=fold, variant=variant) | hyper

    return cached(f"{cohort}__{seed}__{fold}__{variant}", compute)


def theory_fold(name, X, y, donor, tr, te, s, gamma, l1, seed, fold):
    # one training subset: alpha tuned by inner CV at the fixed gamma and l1_ratio; predictions and non-zero proteins
    def compute():
        inner = group_folds(donor[tr], N_INNER, 1000 * seed + fold)
        w, b, hyper = refit(X[tr], y[tr], s, {(gamma, l1): path_mse(X[tr], y[tr], inner, penalty(s, gamma), l1)})
        return X[te] @ w + b, hyper["n_nonzero"]

    return cached(name, compute)


def donor_scores(units, genes, fits, species):
    # donor mean over the units of the raw predictions, and of the unit-standardised ones (script 52's scores), of the
    # variants and of the published clocks (utils.clock_scores)
    raw, std = [], []
    for unit in units:
        x, z, age, donor = unit_data(unit)
        Z = z.reindex(columns=genes).fillna(0).to_numpy()
        pred = pd.DataFrame({v: Z @ w + b for v, (w, b, _) in fits.items()}, index=donor)
        pred = pd.concat([pred, utils.clock_scores(x.T, species, genes=genes).set_axis(donor)], axis=1).assign(age=age)
        raw.append(pred)
        std.append(pred.assign(**{c: standardised(pred[c]) for c in pred if c != "age"}))
    return [pd.concat(x).groupby(level=0).mean() for x in (raw, std)]


# leave-one-cohort-out fits of every cohort and of all cohorts: inner-CV cells in parallel (human first, the slowest),
# then the refits
cohort_keys = {t: ["all", *sorted({cohort_of[u] for u in units_of(t)})] for t in ["human_protein", "mouse_protein"]}
keys = [(t, k, v) for t, ks in cohort_keys.items() for k in ks for v in VARIANTS]
grid = [(*key, *c) for key in keys for c in cells(key[2])]
parallel = Parallel(n_jobs=WORKERS, batch_size=1)
mse = dict(zip(grid, parallel(delayed(loco_cell)(*g) for g in grid), strict=True))
refits = parallel(delayed(loco_fit)(*key, {g[3:]: m for g, m in mse.items() if g[:3] == key}) for key in keys)
FITS, hyper_rows = {}, []
for (target, k, variant), (genes, w, b, hyper) in zip(keys, refits, strict=True):
    FITS.setdefault((target, k), (genes, {}))[1][variant] = (w, b, hyper)
    hyper_rows.append(hyper)

SETS = [("held-out", species, target, units_of(target), False) for species, target in TARGETS.items()]
SETS.append(("external", "mouse", "mouse_protein", utils.EXTERNAL.split(), True))
cohorts, pooled, compared, benchmark, scores = [], [], [], [], []
for scored, species, target, units, external in SETS:
    rows = []
    for k in sorted({cohort_of[u] for u in units}):
        genes, fits = FITS[target, "all" if external else k]
        k_units = [u for u in units if cohort_of[u] == k]
        raw, std = donor_scores(k_units, genes, fits, species)
        donors = raw[["age", *VARIANTS]].rename_axis("donor").reset_index()
        scores.append(donors.assign(set=scored, species=species, cohort=k))
        quantified = [unit_data(u)[0].reindex(columns=genes).notna().any() for u in k_units]
        quantified = pd.concat(quantified, axis=1).any(axis=1).to_numpy()
        names = utils.clock_names(species, genes)
        r_scores = std.drop(columns="age").corr()["elastic_net"]
        weights = {v: w for v, (w, _, _) in fits.items()}
        weights |= {c: names[c].map(w).fillna(0).to_numpy() for c, (_, _, w) in clock_w.items()}
        trained = np.concatenate([unit_data(u)[2] for u in units_of(target) if external or cohort_of[u] != k])
        null = dict(mae_null=(raw.age - trained.mean()).abs().mean(), age_offset=raw.age.mean() - trained.mean())
        for c, w in weights.items():
            test = sps.pearsonr(std[c], std.age)
            mae = (raw[c] - raw.age).abs().mean() if c in VARIANTS else np.nan
            rows.append(dict(set=scored, species=species, cohort=k, weights=c, n=len(std), r=test.statistic))
            rows[-1] |= dict(p=test.pvalue, mae=mae, r_scores=r_scores[c]) | null | dict(age_sd=raw.age.std())
            n_genes = ((w != 0) & quantified).sum()
            benchmark.append(dict(set=scored, species=species, cohort=k, weights=c, n_genes=n_genes, r=test.statistic))
    table = pd.DataFrame(rows)
    table["r_plain"] = table.cohort.map(table[table.weights == "elastic_net"].set_index("cohort").r)
    cohorts.append(table)
    for c, g in table.groupby("weights", sort=False):
        fit = utils.random_effects(np.arctanh(g.r), 1 / np.sqrt(g.n - 3))
        row = dict(set=scored, species=species, weights=c, r=np.tanh(fit["est"]), ci_low=np.tanh(fit["ci_low"]))
        row |= dict(ci_high=np.tanh(fit["ci_high"]), p=fit["p"], i2=fit["i2"], k=fit["k"])
        pooled.append(row)
        benchmark.append(
            dict(set=scored, species=species, cohort="pooled", weights=c)
            | {x: row[x] for x in ["r", "ci_low", "ci_high"]}
        )
        if c in VARIANTS[1:]:
            g = g[g.r_scores < 1]  # cohorts where the two scores agree exactly carry no difference
            d = np.arctanh(g.r) - np.arctanh(g.r_plain)
            diff = utils.random_effects(d, meng_se(g.r, g.r_plain, g.r_scores, g.n))
            compared.append(
                dict(set=scored, species=species, weights=c, control="elastic_net", difference=diff["est"]) | diff
            )

cohorts, pooled, compared = pd.concat(cohorts), pd.DataFrame(pooled), pd.DataFrame(compared).drop(columns="est")
utils.write(cohorts, OUT / "cohorts.tsv")
utils.write(pd.concat(scores), OUT / "scores.tsv.gz")
utils.write(pooled, OUT / "pooled.tsv")
utils.write(compared, OUT / "compared.tsv")
new = pd.DataFrame(benchmark).assign(group=lambda d: d.weights.map(GROUPS).fillna("published (trained elsewhere)"))
benchmark = pd.concat([utils.read(S52 / "benchmark.tsv"), new]).astype({"n_genes": "Int64"})
benchmark["tier"] = benchmark.weights.map(tiers).fillna("protein clocks")
utils.write(benchmark, OUT / "benchmark.tsv")
ranked = benchmark.query("cohort == 'pooled'").drop(columns=["cohort", "n_genes"])
ranked["rank"] = ranked.groupby(["set", "species"]).r.rank(ascending=False).astype(int)
ranked["rank_in_tier"] = ranked.groupby(["set", "species", "tier"]).r.rank(ascending=False).astype(int)
utils.write(ranked.sort_values(["set", "species", "tier", "rank_in_tier"]), OUT / "benchmark_pooled.tsv")
print(cohorts.query("weights in @VARIANTS").round(3).to_string(index=False))
print(pooled.round(4).to_string(index=False))
print(compared.round(4).to_string(index=False))

# paired native modality: the RNA clocks on the RNA and the protein models on the protein of the same samples, the
# difference of the m6Brain clock to each tested on the shared donors (dependent correlations, Meng's SE)
rna_clocks = [c for c, t in tiers.items() if t == "RNA clocks applied to protein"]
protein_models = ["m6brain_clock", "elastic_net", *[c for c, t in tiers.items() if t == "protein clocks"]]
modality = []
for unit, (species, target) in PAIRED.items():
    cohort = cohort_of[unit]
    genes, fits = FITS[target, cohort]
    protein = donor_scores([unit], genes, fits, species)[1]
    x, rna, samples, _ = utils.load_unit(unit)
    donor = samples.donor_id.astype(str).loc[x.columns].to_numpy()
    on_rna = utils.clock_scores(x, species, rna=rna, genes=genes)[rna_clocks].set_axis(donor)
    scores = protein[[*protein_models, "age"]].join(on_rna.groupby(level=0).mean(), how="inner")
    for c in [*protein_models, *rna_clocks]:
        test = sps.pearsonr(scores[c], scores.age)
        modality.append(dict(cohort=cohort, species=species, modality="RNA" if c in rna_clocks else "protein"))
        modality[-1] |= dict(weights=c, tier=tiers.get(c, "protein clocks"), n=len(scores), r=test.statistic)
        modality[-1] |= dict(p=test.pvalue, r_scores=sps.pearsonr(scores[c], scores.m6brain_clock).statistic)
modality = pd.DataFrame(modality)
modality["r_m6brain"] = modality.cohort.map(modality.query("weights == 'm6brain_clock'").set_index("cohort").r)
other = modality.weights != "m6brain_clock"
d = np.arctanh(modality.r_m6brain) - np.arctanh(modality.r)
modality.loc[other, "difference"] = d[other]
modality.loc[other, "se"] = meng_se(modality.r_m6brain, modality.r, modality.r_scores, modality.n)[other]
modality["p_difference"] = 2 * sps.norm.sf(abs(modality.difference / modality.se))
utils.write(modality, OUT / "native_modality.tsv")
print(modality.round(3).to_string(index=False))

# within-cohort nested CV: absolute age in real units and seed variability; the plain net first, the slowest last
data, jobs = {}, []
for cohort, target in WITHIN.items():
    genes = utils.read(S52 / "runs" / f"{target}__{cohort}.parquet").index
    X, y, donor, _ = stack([u for u in units_of(target) if cohort_of[u] == cohort], genes)
    data[cohort] = X, y, donor, priors_of(target, cohort, genes)
    for s in range(N_SEEDS):
        for i, (tr, te) in enumerate(group_folds(donor, N_OUTER, utils.SEED + s)):
            jobs += [(cohort, s, i, tr, te, v) for v in VARIANTS]
jobs.sort(key=lambda j: (VARIANTS.index(j[5]), list(WITHIN).index(j[0])))
runs = parallel(
    delayed(outer_fold)(c, *data[c][:3], tr, te, data[c][3][v], v, utils.SEED + s, i) for c, s, i, tr, te, v in jobs
)
seeds, paired = [], []
for cohort, target in WITHIN.items():
    _, y, donor, _ = data[cohort]
    pred, null = np.zeros((N_SEEDS, len(y), len(VARIANTS))), np.zeros((N_SEEDS, len(y)))
    for (c, s, _, tr, te, v), (p, hyper) in zip(jobs, runs, strict=True):
        if c == cohort:
            pred[s, te, VARIANTS.index(v)] = p
            null[s, te] = y[tr].mean()  # the null predicts the mean age of the training fold
            hyper_rows.append(hyper | dict(target=target, cohort=cohort))
    error = {v: [] for v in VARIANTS}
    for s in range(N_SEEDS):
        by_donor = pd.DataFrame(pred[s], index=donor, columns=VARIANTS).assign(age=y, null=null[s])
        by_donor = by_donor.groupby(level=0).mean()
        for v in VARIANTS:
            error[v].append((by_donor[v] - by_donor.age).abs())
            seeds.append(dict(target=target, cohort=cohort, variant=v, seed=utils.SEED + s, n_donors=len(by_donor)))
            seeds[-1] |= dict(mae=error[v][-1].mean(), r=sps.pearsonr(by_donor[v], by_donor.age).statistic)
            seeds[-1] |= dict(mae_null=(by_donor.null - by_donor.age).abs().mean())
    error = {v: pd.concat(e, axis=1).mean(axis=1) for v, e in error.items()}  # per donor, over the seeds
    for v in VARIANTS[1:]:
        paired.append(dict(cohort=cohort, variant=v, mae_vs_unguided=(error[v] - error["elastic_net"]).mean()))
        paired[-1] |= dict(p_wilcoxon=sps.wilcoxon(error[v], error["elastic_net"]).pvalue)
seeds = pd.DataFrame(seeds)
within = (
    seeds.groupby(["target", "cohort", "variant"], sort=False)
    .agg(
        n_donors=("n_donors", "first"),
        mae_mean=("mae", "mean"),
        mae_sd=("mae", "std"),
        r_mean=("r", "mean"),
        r_sd=("r", "std"),
        mae_null_mean=("mae_null", "mean"),
        mae_null_sd=("mae_null", "std"),
    )
    .reset_index()
    .merge(pd.DataFrame(paired), how="left", on=["cohort", "variant"])
)
utils.write(seeds, OUT / "within_cohort_seeds.tsv")
utils.write(within, OUT / "within_cohort.tsv")
utils.write(pd.DataFrame(hyper_rows), OUT / "hyperparameters.tsv")
print(within.round(3).to_string(index=False))

# where the m6Brain clock (and its label-free mode) falls short: the r of a cohort with a CI that includes 0, an MAE
# not below that of the training-mean null, a published clock of higher r
key = ["set", "species", "cohort"]
modes = ["m6brain_clock", "model"]
by_cohort = pd.concat([cohorts, utils.read(S52 / "cohorts.tsv").query("weights == 'model'")], ignore_index=True)
others = by_cohort[
    by_cohort.weights.isin([c for c in clock_w if tiers[c] in ("protein clocks", "RNA clocks applied to protein")])
]
best = others.loc[others.groupby(key).r.idxmax()].set_index(key)
failures = []
for _, row in by_cohort[by_cohort.weights.isin(modes)].iterrows():
    half = utils.Z95 / np.sqrt(row.n - 3)
    lo, hi = np.tanh(np.arctanh(row.r) - half), np.tanh(np.arctanh(row.r) + half)
    where = dict(set=row["set"], species=row.species, cohort=row.cohort, weights=row.weights)
    if lo <= 0 <= hi:
        failures.append(
            where | dict(criterion="r CI includes 0", reason=f"{row.n} donors: r {row.r:.2f}, CI {lo:.2f} to {hi:.2f}")
        )
    if row.weights == "m6brain_clock" and row.mae >= row.mae_null:
        reason = f"MAE {row.mae:.2f} against {row.mae_null:.2f}; mean age {row.age_offset:+.1f} from the training mean"
        reason += f", SD of age {row.age_sd:.1f}"
        failures.append(where | dict(criterion="MAE not below the training-mean null", reason=reason))
    top = best.loc[(row["set"], row.species, row.cohort)]
    if top.r > row.r:
        n_above = (others[(others[key] == row[key]).all(axis=1)].r > row.r).sum()
        failures.append(
            where
            | dict(
                criterion="below a published clock",
                reason=f"{top.weights} r {top.r:.2f} against {row.r:.2f}; {n_above} published clocks higher",
            )
        )
for _, r in within.query("variant == 'm6brain_clock' and mae_mean >= mae_null_mean").iterrows():
    reason = f"MAE {r.mae_mean:.2f} against {r.mae_null_mean:.2f}"
    failures.append(
        dict(
            set="within-cohort",
            species=r.target.split("_")[0],
            cohort=r.cohort,
            weights=r.variant,
            criterion="MAE not below the training-mean null",
            reason=reason,
        )
    )
failures = pd.DataFrame(failures, columns=[*key, "weights", "criterion", "reason"])
utils.write(failures, OUT / "failures.tsv")
print(failures.to_string(index=False))

# theory: does the guidance act as prior information for few samples? Luo only, outer 10-fold by donor x 5 seeds, the
# training fold subsampled to n_train donors (nested subsets, the held-out fold in full); unguided, m6Brain-guided and
# the same prior values permuted across the genes (a new permutation per seed and fold); l1_ratio and gamma fixed at
# the values chosen most often in the nested run (unguided: l1_ratio 0.9; guided and shuffled: l1_ratio 0.9, gamma 2),
# alpha tuned on the same path by inner 5-fold CV by donor
N_TRAIN, THEORY_SEEDS, SHUFFLED = [25, 50, 100, None], 5, "shuffled_prior"
cohort = "luo2024_syn32136022"
X, y, donor, priors = data[cohort]
chosen = pd.DataFrame(hyper_rows).query("stage == 'nested' and cohort == @cohort")
chosen = chosen.groupby("variant")[["l1_ratio", "gamma"]].agg(lambda x: x.mode()[0])
setting = {v: chosen.loc[v] for v in VARIANTS[:2]} | {SHUFFLED: chosen.loc["m6brain_clock"]}
jobs = []
for s in range(THEORY_SEEDS):
    for i, (tr, te) in enumerate(group_folds(donor, N_OUTER, utils.SEED + s)):
        order = np.random.default_rng([utils.SEED, s, i, 0]).permutation(tr)
        prior = {"elastic_net": priors["elastic_net"], "m6brain_clock": priors["m6brain_clock"]}
        prior[SHUFFLED] = np.random.default_rng([utils.SEED, s, i, 1]).permutation(priors["m6brain_clock"])
        for n in N_TRAIN:
            jobs += [(v, n, s, i, order[:n], te, prior[v]) for v in setting]
runs = parallel(
    delayed(theory_fold)(f"theory__{s}__{i}__{n}__{v}", X, y, donor, sub, te, p, *setting[v][::-1], utils.SEED + s, i)
    for v, n, s, i, sub, te, p in jobs
)
pred, nonzero = {}, {}
for (v, n, s, _, _, te, _), (p_te, n_nonzero) in zip(jobs, runs, strict=True):
    pred.setdefault((v, n), np.zeros((THEORY_SEEDS, len(y))))[s, te] = p_te
    nonzero.setdefault((v, n), []).append(n_nonzero)
against = {"m6brain_clock": {"unguided": "elastic_net", "shuffled": SHUFFLED}, SHUFFLED: {"unguided": "elastic_net"}}
theory = []
for n in N_TRAIN:
    error = {v: np.abs(pred[v, n] - y) for v in setting}  # seeds x donors
    for v in setting:
        row = dict(n_train=n or "all", variant=v, mae_mean=error[v].mean(1).mean(), mae_sd=error[v].mean(1).std(ddof=1))
        r = [sps.pearsonr(p, y).statistic for p in pred[v, n]]
        row |= dict(r_mean=np.mean(r), r_sd=np.std(r, ddof=1), n_nonzero=np.mean(nonzero[v, n]))
        for name, other in against.get(v, {}).items():
            donors = error[v].mean(0), error[other].mean(0)  # per donor over the seeds
            row |= {
                f"mae_vs_{name}": (donors[0] - donors[1]).mean(),
                f"p_wilcoxon_{name}": sps.wilcoxon(*donors).pvalue,
            }
        theory.append(row)
theory = pd.DataFrame(theory)
utils.write(theory, OUT / "theory.tsv")
print(theory.round(3).to_string(index=False))
