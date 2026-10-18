# Individual m6A load index and reader titration. The index of a sample is the OLS slope over genes of the standardised
# protein abundance on the m6A load of script 21 with LOAD_COVARIATES (genes complete in the unit, as script 33).
# (1) Spearman rho of the index with age in each cohort (donor mean of the unit-standardised index), pooled by species.
# (2) Per unit, OLS of the index (reader genes left out) on the reader level, age and the unit's design; the reader
# coefficient pooled over the cohorts of a species (random effects; BH over the four readers of a species) and ranked
# among the control RBPs of script 33. Reader level: z of the protein in the unit (a module: mean of its members' z),
# a missing value at the unit mean; proteins quantified in fewer than half of the samples are skipped.
# reads the prepared units (08), scripts 09 and 33; writes output/54_reader_titration/
import numpy as np
import pandas as pd
import statsmodels.api as sm
from joblib import Parallel, delayed
from scipy import stats as sps

import model
import utils

OUT = utils.OUT / "54_reader_titration"
READERS = {
    "FMR1": ["FMR1"],
    "FMRP family": ["FMR1", "FXR1", "FXR2"],
    "YTHDF module": ["YTHDF1", "YTHDF2", "YTHDF3"],
    "YTHDF3": ["YTHDF3"],
}
LOAD = {"mouse": "n_sites_m6A", "human": "n_sites_m6A"}  # human: brain peaks (script 09)
CONTROLS = utils.read(utils.OUT / "33_regulator_network" / "control_rbps.tsv").query("selected").rbp.tolist()
REGULATORS = READERS | {c: [c] for c in CONTROLS}
UNITS = {
    species: next(u.split() for n, layer, u in utils.CONTEXTS if (n, layer) == (f"{species}_protein", "protein_only"))
    for species in LOAD
}
cohort_of = utils.read(utils.PREPARED / "units.tsv", index_col="unit").dataset
features = {
    s: utils.read(utils.OUT / "09_transcript_features" / f"{s}_features.tsv.gz").set_index("gene_id") for s in LOAD
}
gene_id = {
    s: pd.Series(f.index, index=f.symbol.str.upper()).pipe(lambda v: v[~v.index.duplicated(keep=False)])
    for s, f in features.items()
}
to_mouse = utils.orthologs().set_index("human").mouse


def z(v):
    return (v - v.mean()) / v.std()


def unit_fits(unit, species):
    protein, _, samples, covariates = utils.load_unit(unit)
    complete = protein.index[protein.notna().all(axis=1)].intersection(features[species].index)
    genes = features[species].loc[complete].join(protein.mean(axis=1).rename("level"))
    genes = genes.dropna(subset=[LOAD[species], *utils.LOAD_COVARIATES])
    Y = protein.loc[genes.index].apply(z, axis=1).to_numpy(float)
    X = np.column_stack([np.ones(len(genes)), np.log1p(genes[LOAD[species]]), genes[utils.LOAD_COVARIATES]])
    X = X.astype(float)
    age = samples.age_months if species == "mouse" else samples.age_years
    slope = z(np.linalg.pinv(X)[1] @ Y)
    index = pd.DataFrame(dict(unit=unit, donor=samples.donor_id.astype(str), age=age, index=slope))
    D = model.design(samples, covariates)
    rows = []
    for name, members in REGULATORS.items():
        ids = [gene_id[species].get(m if species == "human" else to_mouse.get(m)) for m in members]
        quantified = [i for i in ids if i in protein.index and protein.loc[i].notna().mean() >= 0.5]
        if not quantified:
            continue
        level = np.mean([z(protein.loc[i].fillna(protein.loc[i].mean()).to_numpy(float)) for i in quantified], axis=0)
        keep = ~genes.index.isin(ids)
        fit = sm.OLS(np.linalg.pinv(X[keep])[1] @ Y[keep], D.assign(reader=level)).fit()
        kind = "reader" if name in READERS else "control"
        rows.append(dict(unit=unit, species=species, regulator=name, kind=kind, n_samples=len(samples))
                    | dict(n_genes=keep.sum(), est=fit.params["reader"], se=fit.bse["reader"], p=fit.pvalues["reader"]))
    return index, rows


jobs = [delayed(unit_fits)(u, s) for s, units in UNITS.items() for u in units]
results = Parallel(n_jobs=16)(jobs)
species_of = {u: s for s, units in UNITS.items() for u in units}

# (1) index against age
index = pd.concat([r[0] for r in results]).assign(cohort=lambda d: d.unit.map(cohort_of))
index["species"] = index.unit.map(species_of)
donors = index.groupby(["species", "cohort", "donor"])[["index", "age"]].mean().reset_index()
cohorts = pd.DataFrame(
    [
        dict(species=s, cohort=k, n=len(d), rho=sps.spearmanr(d["index"], d.age).statistic)
        for (s, k), d in donors.groupby(["species", "cohort"])
    ]
)
pooled = []
for species, g in cohorts.groupby("species"):
    fit = utils.random_effects(np.arctanh(g.rho), 1 / np.sqrt(g.n - 3))
    pooled.append(dict(species=species, rho=np.tanh(fit["est"]), ci_low=np.tanh(fit["ci_low"])) | dict(
        ci_high=np.tanh(fit["ci_high"]), p=fit["p"], i2=fit["i2"], k=fit["k"], n_negative=int((g.rho < 0).sum())
    ))
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
print(cohorts.round(3).to_string(index=False))
print(pd.read_csv(OUT / "index_pooled.tsv", sep="\t").round(4).to_string(index=False))
print(readers.round(4).to_string(index=False))
