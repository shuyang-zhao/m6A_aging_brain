# Zero-shot proteome clock: the sample score is the sum over genes of the TabPFN-predicted age z (script 50 model,
# refitted without the scored cohort) times the protein's standardised abundance in the sample; the Spearman rho of
# the score with age in each held-out cohort, pooled by species, against control weights (gate G1 of
# dev/gates/criteria.md); the mouse external validation set (utils.EXTERNAL) scored with the weights of the model of
# all cohorts
#   52_proteome_clock.py inputs   leave-one-cohort-out gene tables (CPU)
#   52_proteome_clock.py TABLE    TabPFN predictions of one table from the full, length-only and no-m6A features (GPU)
#   52_proteome_clock.py          sample scores, correlations with age, pooled and compared (CPU)
# reads scripts 10 and 50, the prepared units and ref/gene_sets/ageing_signatures.tsv; writes output/52_proteome_clock/
import sys

import numpy as np
import pandas as pd
from scipy import stats as sps

import utils

OUT = utils.OUT / "52_proteome_clock"
INPUTS, RUNS = OUT / "inputs", OUT / "runs"
MODEL_INPUTS = utils.OUT / "50_transcript_ageing_model" / "inputs"
TARGETS = {"mouse": "mouse_protein", "human": "human_protein"}  # context of the protein_only layer
N_OUTER, N_RANDOM = 5, 1000
LENGTHS = ["log_tx_len", "log_gene_len", "log1p_intron_len", "log1p_utr5_len", "log_cds_len", "log1p_utr3_len"]
LABELS = ["symbol", "chrom", "fold", "target", "shuffled"]
MODELS = ["model", "length_only", "without_m6a"]  # TabPFN weights; the others are made here
# compared with the model per cohort (the no-m6A scores equal the model's in rank in most cohorts: not compared)
CONTROLS = ["length_only", "m6a_load", "DEMAGALHAES_AGING", "LU_AGING_BRAIN", "observed"]


def units_of(target):
    return next(units.split() for name, layer, units in utils.CONTEXTS if (name, layer) == (target, "protein_only"))


cohort_of = utils.read(utils.PREPARED / "units.tsv", index_col="unit").dataset

if sys.argv[1:] == ["inputs"]:
    # the target and expression of script 50 recomputed from the other cohorts; features and folds of script 50
    effects = utils.read(
        utils.OUT / "10_protein_age_effects" / "protein_only.tsv.gz",
        usecols=["unit", "gene_id", "n", "t", "p", "status"],
    ).query("status == 'ok'")
    effects = effects.assign(z=np.sign(effects.t) * sps.norm.isf(effects.p / 2), cohort=effects.unit.map(cohort_of))
    INPUTS.mkdir(parents=True, exist_ok=True)
    for target in TARGETS.values():
        base = pd.read_parquet(MODEL_INPUTS / f"{target}.parquet")
        units = units_of(target)
        scores = {}
        for unit in units:
            level = utils.load_unit(unit)[0].mean(axis=1)
            scores[unit] = pd.Series(sps.norm.ppf((level.rank() - 0.5) / level.count()), index=level.index)
        table = effects[effects.unit.isin(units)]
        base.to_parquet(INPUTS / f"{target}__all.parquet")  # every cohort: the table of script 50
        for k in sorted({cohort_of[u] for u in units}):
            rest = table[table.cohort != k]
            cohorts = rest.groupby(["gene_id", "cohort"])[["z", "n"]].mean()
            z = cohorts.groupby("gene_id").apply(lambda d: utils.stouffer(d.z, np.sqrt(d.n))["z"])
            expression = pd.concat([s for u, s in scores.items() if cohort_of[u] != k], axis=1).mean(axis=1)
            loco = base.drop(columns=["target", "expression"]).join(z.rename("target"), how="inner")
            loco = loco.assign(expression=expression.reindex(loco.index))[base.columns]
            loco.to_parquet(INPUTS / f"{target}__{k}.parquet")
            agreement = loco.target.corr(base.target, "spearman")
            print(target, k, loco.shape, "rho with the all-cohort target", round(agreement, 3))

elif sys.argv[1:]:
    # TabPFN-3.5 as script 50: each chromosome group predicted in context of the other groups
    from tabpfn import TabPFNRegressor
    from tabpfn.constants import ModelVersion

    name = sys.argv[1]
    table = pd.read_parquet(INPUTS / f"{name}.parquet")
    X, y, fold = table.drop(columns=[c for c in LABELS if c in table]), table.target.to_numpy(), table.fold.to_numpy()
    m6a = [c for c in X if c.startswith(("m6A_", "bound_"))]
    features = {"model": list(X), "length_only": LENGTHS, "without_m6a": [c for c in X if c not in m6a]}
    runs = table[["symbol", "fold", "target"]].copy()
    for model, columns in features.items():
        runs[model] = np.nan
        for k in range(N_OUTER):
            fit = TabPFNRegressor.create_default_for_version(ModelVersion.V3_5)
            fit.set_params(device="cuda", random_state=utils.SEED, ignore_pretraining_limits=True)
            fit.fit(X.loc[fold != k, columns], y[fold != k])
            runs.loc[fold == k, model] = fit.predict(X.loc[fold == k, columns])
        print(name, model, "held-out rho", round(sps.spearmanr(y, runs[model]).statistic, 3), flush=True)
    utils.write(runs, RUNS / f"{name}.parquet")

else:
    signatures = utils.read(utils.ROOT / "ref/gene_sets/ageing_signatures.tsv")  # MSigDB, downloaded 2026-10-17
    signatures = signatures.assign(sign=np.where(signatures.gene_set.str.endswith("_UP"), 1, -1))
    signatures["gene_set"] = signatures.gene_set.str.replace(r"_(UP|DN)$", "", regex=True)
    to_human = utils.orthologs().set_index("mouse").human
    rng = np.random.default_rng(utils.SEED)

    def standardised(x):
        return (x - x.mean()) / x.std()

    def weights(target, k):
        # genes x weights of the held-out cohort k: TabPFN, linear load, signatures, observed z of the other cohorts
        runs = utils.read(RUNS / f"{target}__{k}.parquet")
        load = pd.read_parquet(INPUTS / f"{target}__{k}.parquet").m6A_all
        w = runs[MODELS].assign(m6a_load=np.polyfit(load, runs.target, 1)[0] * load, observed=runs.target)
        symbol = runs.symbol.str.upper() if target.startswith("human") else runs.symbol.str.upper().map(to_human)
        for name, rows in signatures.groupby("gene_set"):
            w[name] = symbol.map(rows.drop_duplicates("symbol").set_index("symbol").sign).fillna(0)
        return w

    def donor_scores(units, w):
        # unit scores (proteins standardised in the unit, missing = 0; score standardised), donor mean over units
        out = []
        for unit in units:
            protein, _, samples, _ = utils.load_unit(unit)
            z = protein.reindex(w.index).T.pipe(standardised).fillna(0)
            score = (z @ w).pipe(standardised)
            age = samples.age_months if "age_months" in samples else samples.age_years
            out.append(score.assign(donor=samples.donor_id.astype(str), age=age))
        return pd.concat(out).groupby("donor").mean()

    def rho(scores, age):
        return np.array([sps.spearmanr(s, age).statistic for s in scores.T])

    def meng_se(r1, r2, r12, n):
        # SE of atanh r1 - atanh r2 for two correlations with a shared variable (Meng, Rosenthal & Rubin 1992)
        mean_sq = (r1**2 + r2**2) / 2
        f = np.minimum((1 - r12) / (2 * (1 - mean_sq)), 1)
        h = (1 - f * mean_sq) / (1 - mean_sq)
        return np.sqrt(2 * (1 - r12) * h / (n - 3))

    # held-out: each cohort of the context with the weights fitted without it; external: the external cohorts with the
    # weights of all cohorts
    SETS = [("held-out", species, target, units_of(target), False) for species, target in TARGETS.items()]
    SETS.append(("external", "mouse", "mouse_protein", utils.EXTERNAL.split(), True))
    cohorts, pooled, compared, nulls, scores_out = [], [], [], [], []
    for scored, species, target, units, external in SETS:
        universe = pd.read_parquet(MODEL_INPUTS / f"{target}.parquet").index
        keys = rng.random((N_RANDOM, len(universe)))  # one random gene order per draw, shared by the cohorts
        rows, null = [], []
        for k in sorted({cohort_of[u] for u in units}):
            w = weights(target, "all" if external else k)
            ranks = keys[:, universe.get_indexer(w.index)].argsort(axis=1).argsort(axis=1)
            random = pd.DataFrame(np.sort(w.model.to_numpy())[ranks].T, index=w.index)
            scores = donor_scores([u for u in units if cohort_of[u] == k], w.join(random.add_prefix("random_")))
            n = len(scores)
            observed = rho(scores[list(w)].to_numpy(), scores.age)
            null.append(rho(scores.filter(like="random_").to_numpy(), scores.age))
            r12 = scores[list(w)].corr("spearman").model
            rows += [
                dict(set=scored, species=species, cohort=k, weights=c, n=n, rho=r)
                | dict(p=sps.spearmanr(scores[c], scores.age).pvalue)
                | dict(r_model=r12[c], rho_model=observed[0])
                for c, r in zip(w, observed, strict=True)
            ]
            scores_out.append(scores[["age", *w]].assign(set=scored, species=species, cohort=k))
        table = pd.DataFrame(rows)
        cohorts.append(table)
        # pooled rho of each weight set on the Fisher scale, and of the random draws
        null_fisher = np.arctanh(np.array(null))  # cohorts x draws
        n_cohort = table.drop_duplicates("cohort").n.to_numpy()
        se = 1 / np.sqrt(n_cohort - 3)
        null_pooled = [utils.random_effects(null_fisher[:, d], se)["est"] for d in range(N_RANDOM)]
        draws = np.arange(1, N_RANDOM + 1)
        nulls.append(pd.DataFrame(dict(set=scored, species=species, draw=draws, rho=np.tanh(null_pooled))))
        for c, g in table.groupby("weights", sort=False):
            fit = utils.random_effects(np.arctanh(g.rho), 1 / np.sqrt(g.n - 3))
            row = dict(set=scored, species=species, weights=c, rho=np.tanh(fit["est"]), ci_low=np.tanh(fit["ci_low"]))
            row |= dict(ci_high=np.tanh(fit["ci_high"]), p=fit["p"], i2=fit["i2"], k=fit["k"])
            if c == "model":
                row["p_random"] = utils.perm_pvalue(fit["est"], null_pooled, "greater")
            pooled.append(row)
            if c in CONTROLS:
                g = g[g.r_model < 1]  # cohorts where the two scores rank the donors identically carry no difference
                d = np.arctanh(g.rho_model) - np.arctanh(g.rho)
                diff = utils.random_effects(d, meng_se(g.rho_model, g.rho, g.r_model, g.n))
                compared.append(dict(set=scored, species=species, control=c, difference=diff["est"]) | diff)

    OUT.mkdir(parents=True, exist_ok=True)
    utils.write(pd.concat(cohorts), OUT / "cohorts.tsv")
    utils.write(pd.DataFrame(pooled), OUT / "pooled.tsv")
    utils.write(pd.DataFrame(compared).drop(columns="est"), OUT / "compared.tsv")
    utils.write(pd.concat(nulls), OUT / "random_null.tsv")
    utils.write(pd.concat(scores_out).rename_axis("donor").reset_index(), OUT / "scores.tsv.gz")
    print(pd.concat(cohorts).query("weights == 'model'").round(3).to_string(index=False))
    print(pd.DataFrame(pooled).round(4).to_string(index=False))
    print(pd.DataFrame(compared).round(4).to_string(index=False))
