# Individual m6A load index and reader titration. The index of a sample is the OLS slope over genes of the standardised
# protein abundance on the m6A load of script 21 with LOAD_COVARIATES (genes complete in the unit, as script 33).
# (1) Pearson r of the index with age in each cohort (donor mean of the unit-standardised index), pooled by species,
# and in the external mouse proteomes (utils.EXTERNAL; index only).
# (2) Per unit, OLS of the index (reader genes left out; standardised in the unit, so the coefficient is SD of index per
# SD of reader) on the reader level, age and the unit's design; the reader coefficient pooled over the cohorts of a
# species (random effects; BH over the four readers of a species) and ranked among the control RBPs of script 33.
# Reader level: z of the protein in the unit (a module: mean of its members' z), a missing value at the unit mean;
# proteins quantified in fewer than half of the samples are skipped.
# (3) Bridge to the clock: per cohort, OLS of the age acceleration of the m6Brain clock (residual of the donor's
# predicted age on age, months in mouse and years in human; script 55 weights without the cohort, external proteomes:
# all cohorts) on the index residual on age, standardised: acceleration per SD of the index, pooled by species (random
# effects).
# reads the prepared units (08), scripts 09, 33 and 55; writes output/54_reader_titration/
import numpy as np
import pandas as pd
import statsmodels.api as sm
from joblib import Parallel, delayed
from scipy import stats as sps

import model
import utils

OUT = utils.OUT / "54_reader_titration"
CLOCK = utils.OUT / "55_m6brain_clock" / "weights"
READERS = {
    "FMR1": ["FMR1"],
    "FMRP family": ["FMR1", "FXR1", "FXR2"],
    "YTHDF module": ["YTHDF1", "YTHDF2", "YTHDF3"],
    "YTHDF3": ["YTHDF3"],
}
CONTROLS = utils.read(utils.OUT / "33_regulator_network" / "control_rbps.tsv").query("selected").rbp.tolist()
REGULATORS = READERS | {c: [c] for c in CONTROLS}
UNITS = {
    species: next(u.split() for n, layer, u in utils.CONTEXTS if (n, layer) == (f"{species}_protein", "protein_only"))
    for species in ("mouse", "human")
}
UNITS["external"] = utils.EXTERNAL.split()
cohort_of = utils.read(utils.PREPARED / "units.tsv", index_col="unit").dataset
features = {
    s: utils.read(utils.OUT / "09_transcript_features" / f"{s}_features.tsv.gz").set_index("gene_id")
    for s in ("mouse", "human")
}
gene_id = {
    s: pd.Series(f.index, index=f.symbol.str.upper()).pipe(lambda v: v[~v.index.duplicated(keep=False)])
    for s, f in features.items()
}
to_mouse = utils.orthologs().set_index("human").mouse


def z(v):
    return (v - v.mean()) / v.std()


def pool_r(g):
    # Fisher scale, random effects over the cohorts of g
    fit = utils.random_effects(np.arctanh(g.r), 1 / np.sqrt(g.n - 3))
    return dict(r=np.tanh(fit["est"]), ci_low=np.tanh(fit["ci_low"]), ci_high=np.tanh(fit["ci_high"])) | dict(
        p=fit["p"], i2=fit["i2"], k=fit["k"], n_negative=int((g.r < 0).sum())
    )


def unit_fits(unit, group):
    species = "human" if group == "human" else "mouse"
    protein, _, samples, covariates = utils.load_unit(unit)
    genes, X, Y = utils.load_design(protein, features[species])
    age = samples.age_months if species == "mouse" else samples.age_years
    slope = z(np.linalg.pinv(X)[1] @ Y)
    index = pd.DataFrame(dict(unit=unit, donor=samples.donor_id.astype(str), age=age, index=slope))
    if group == "external":
        return index, []
    D = model.design(samples, covariates)
    rows = []
    for name, members in REGULATORS.items():
        ids = [gene_id[species].get(m if species == "human" else to_mouse.get(m)) for m in members]
        quantified = [i for i in ids if i in protein.index and protein.loc[i].notna().mean() >= 0.5]
        if not quantified:
            continue
        level = np.mean([z(protein.loc[i].fillna(protein.loc[i].mean()).to_numpy(float)) for i in quantified], axis=0)
        keep = ~genes.index.isin(ids)
        fit = sm.OLS(z(np.linalg.pinv(X[keep])[1] @ Y[keep]), D.assign(reader=level)).fit()
        kind = "reader" if name in READERS else "control"
        rows.append(dict(unit=unit, species=species, regulator=name, kind=kind, n_samples=len(samples))
                    | dict(n_genes=keep.sum(), est=fit.params["reader"], se=fit.bse["reader"], p=fit.pvalues["reader"]))
    return index, rows


jobs = [delayed(unit_fits)(u, g) for g, units in UNITS.items() for u in units]
results = Parallel(n_jobs=16)(jobs)
group_of = {u: g for g, units in UNITS.items() for u in units}

# (1) index against age
index = pd.concat([r[0] for r in results]).assign(cohort=lambda d: d.unit.map(cohort_of))
index["species"] = index.unit.map(group_of)  # mouse, human or external (mouse)
donors = index.groupby(["species", "cohort", "donor"])[["index", "age"]].mean().reset_index()
cohorts = pd.DataFrame(
    [
        dict(species=s, cohort=k, n=len(d), r=sps.pearsonr(d["index"], d.age).statistic)
        for (s, k), d in donors.groupby(["species", "cohort"])
    ]
)
pooled = [dict(species=species) | pool_r(g) for species, g in cohorts.groupby("species")]
utils.write(donors, OUT / "index_donors.tsv.gz")
utils.write(cohorts, OUT / "index_cohorts.tsv")
utils.write(pd.DataFrame(pooled), OUT / "index_pooled.tsv")

# (2) reader titration
coefficients = pd.DataFrame([row for r in results for row in r[1]])
# units of a cohort averaged, cohorts combined by random effects (as the load model of script 21)
pooled = []
for (species, name, kind), g in coefficients.groupby(["species", "regulator", "kind"]):
    cohorts = g.assign(cohort=g.unit.map(cohort_of)).groupby("cohort")[["est", "se"]].mean()
    fit = utils.random_effects(cohorts.est, cohorts.se)
    pooled.append(dict(species=species, regulator=name, kind=kind, n_units=len(g), n_cohorts=len(cohorts)) | {
        k: fit[k] for k in ("est", "se", "ci_low", "ci_high", "p", "i2")
    })
pooled = pd.DataFrame(pooled)
readers = pooled[pooled.kind == "reader"].copy()
controls = {s: g.est for s, g in pooled[pooled.kind == "control"].groupby("species")}
readers["q"] = readers.groupby("species").p.transform(utils.bh)
readers["control_pct"] = [100 * (controls[r.species] < r.est).mean() for r in readers.itertuples()]
readers["n_controls"] = readers.species.map(lambda s: len(controls[s]))
utils.write(coefficients, OUT / "titration_units.tsv")
utils.write(pooled[pooled.kind == "control"], OUT / "titration_controls.tsv")
utils.write(readers, OUT / "titration.tsv")

# (3) bridge to the clock
TARGET = dict(mouse="mouse_protein", human="human_protein", external="mouse_protein")


def residual(y, age):
    return y - np.polyval(np.polyfit(age, y, 1), age)


def predicted_age(units, w):
    genes, out = w.drop("(intercept)"), []
    for unit in units:
        protein, _, samples, _ = utils.load_unit(unit)
        predicted = w["(intercept)"] + protein.reindex(genes.index).T.pipe(z).fillna(0) @ genes
        out.append(pd.DataFrame(dict(predicted=predicted, donor=samples.donor_id.astype(str))))
    return pd.concat(out).groupby("donor").predicted.mean()


rows = []
for species, units in UNITS.items():
    for k in sorted({cohort_of[u] for u in units}):
        fit = "all" if species == "external" else k
        w = utils.read(CLOCK / f"{TARGET[species]}__{fit}__m6brain_clock.tsv", index_col="gene_id").weight
        d = donors.query("species == @species and cohort == @k").set_index("donor")
        d = d.join(predicted_age([u for u in units if cohort_of[u] == k], w), how="inner")
        X = sm.add_constant(z(residual(d["index"], d.age)).rename("index"))
        fit = sm.OLS(residual(d.predicted, d.age), X).fit()
        unit = "years" if species == "human" else "months"
        rows.append(dict(set=species, cohort=k, n=len(d), unit=unit, est=fit.params["index"], se=fit.bse["index"]))
        rows[-1] |= dict(p=fit.pvalues["index"])
bridge = pd.DataFrame(rows)
bridge_pooled = []
for species, g in bridge.groupby("set"):
    fit = utils.random_effects(g.est, g.se)
    row = dict(set=species, cohort="pooled", n=g.n.sum(), unit=g.unit.iloc[0], n_negative=(g.est < 0).sum())
    bridge_pooled.append(row | {k: fit[k] for k in ("est", "se", "ci_low", "ci_high", "p", "i2", "k")})
bridge = pd.concat([bridge, pd.DataFrame(bridge_pooled)])
utils.write(bridge, OUT / "bridge.tsv")
print(cohorts.round(3).to_string(index=False))
print(pd.read_csv(OUT / "index_pooled.tsv", sep="\t").round(4).to_string(index=False))
print(readers.round(4).to_string(index=False))
print(bridge.round(3).to_string(index=False))
