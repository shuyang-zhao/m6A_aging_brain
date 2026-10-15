# Per-gene age models: design matrix, protected PCs, OLS with limma variance moderation, batch BLUP,
# age effect against a gene property, its heterogeneity (causal forest), cameraPR set test, mixed-model age slopes
import subprocess
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats as sps

from utils import SEED, bh, pct

HERE = Path(__file__).resolve().parent
MIN_RESIDUAL_DF = 3  # fewer residual df: gene is unestimable
PC_NORM_FLOOR = 1e-8  # a protected PC that vanishes after projection is dropped
R_CORES = 12  # cores of the R mixed-model fits


def design(samples, covariates):
    """[Intercept, covariates, age]: numeric columns z-scored, categorical ones dummy coded.

    Constant or aliased covariate columns are dropped; age is the last column.
    """
    columns = {"Intercept": np.ones(len(samples))}
    for name in covariates.get("numeric", []):
        values = pd.to_numeric(samples[name]).to_numpy(float)
        if np.std(values) > 0:
            columns[name] = (values - values.mean()) / values.std()
    for name in covariates.get("categorical", []):
        dummies = pd.get_dummies(samples[name].astype(str), prefix=name, drop_first=True, dtype=float)
        columns.update({column: values.to_numpy() for column, values in dummies.items()})
    kept = {}
    for name, values in columns.items():
        trial = np.column_stack([*kept.values(), values])
        if np.linalg.matrix_rank(trial) == trial.shape[1]:
            kept[name] = values
    X = pd.DataFrame(kept, index=samples.index)
    X["age"] = samples["age"].to_numpy(float)
    return X


def protected_pca(matrix, X, k):
    """Top k principal components (samples x k) of row-standardised residuals after X."""
    X = np.asarray(X, float)
    values = np.asarray(matrix, float) @ (np.eye(len(X)) - X @ np.linalg.pinv(X))
    values = values[np.isfinite(values).all(axis=1)]
    # rows without variation after the projection drop out
    spread = values.std(axis=1)
    values = values[spread > 1e-8] / spread[spread > 1e-8, None]
    vt = np.linalg.svd(values, full_matrices=False)[2]
    return pd.DataFrame(vt[:k].T, index=matrix.columns, columns=[f"PC{i + 1}" for i in range(k)])


def _run_r(script, folder, *args):
    run = subprocess.run(["Rscript", str(HERE / script), str(folder), *map(str, args)], capture_output=True, text=True)
    if run.returncode:
        raise RuntimeError(f"{script} failed:\n{run.stderr}")


def _fit_gene(y, mrna, X, pcs):
    ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
    if mrna is not None:
        ok &= np.isfinite(mrna)
    y, X = y[ok], X[ok]
    covars = X[:, :-1] if mrna is None else np.column_stack([X[:, :-1], mrna[ok]])
    full = np.column_stack([covars, X[:, -1]])
    if pcs is not None and pcs.shape[1]:
        # PCs made orthogonal to this gene's own model so they cannot absorb the age effect
        pc = pcs[ok] - full @ np.linalg.pinv(full) @ pcs[ok]
        full = np.column_stack([covars, pc[:, np.linalg.norm(pc, axis=0) > PC_NORM_FLOOR], X[:, -1]])
    rank = np.linalg.matrix_rank(full)
    df = len(y) - rank
    if rank != full.shape[1] or df < MIN_RESIDUAL_DF:
        return dict(n=len(y), df=df, status="unestimable")
    inv = np.linalg.inv(full.T @ full)
    beta = inv @ full.T @ y
    resid = y - full @ beta
    mrna_coef = np.nan if mrna is None else beta[X.shape[1] - 1]
    return dict(
        n=len(y),
        df=df,
        estimate=beta[-1],
        variance=resid @ resid / df,
        unscaled=inv[-1, -1],
        mrna_coef=mrna_coef,
        status="ok",
    )


def _moderate(table):
    # limma squeezeVar (robust) on the per-gene variances; lmFit cannot take a per-gene mRNA column
    with tempfile.TemporaryDirectory() as folder:
        table[["estimate", "variance", "df", "unscaled"]].to_csv(
            f"{folder}/fit.tsv", sep="\t", index_label="feature_id"
        )
        _run_r("moderate.R", folder)
        return pd.read_csv(f"{folder}/moderated.tsv", sep="\t", index_col="feature_id", dtype={"feature_id": str})


def fit_age(y, samples, covariates, mrna=None, pca=None):
    """Age effect on each row of y (log2 per age unit), given the gene's own mRNA when mrna is set.

    One OLS per gene with limma variance moderation and BH q-values; pca adds protected PCs.
    """
    assert y.columns.equals(samples.index)
    X = design(samples, covariates).to_numpy()
    pcs = None if pca is None else np.asarray(pca, float)
    own = None if mrna is None else mrna.loc[y.index].to_numpy(float)
    rows = [_fit_gene(row, None if own is None else own[i], X, pcs) for i, row in enumerate(y.to_numpy(float))]
    table = pd.DataFrame(rows, index=pd.Index(y.index, name="feature_id"))
    moderated = _moderate(table[table.status == "ok"]).rename(columns={"df_total": "df"})
    table = table.drop(columns="df").join(moderated)
    table["q"] = bh(table.p)
    for column, source in [
        ("effect_pct", "estimate"),
        ("effect_pct_ci_low", "ci_low"),
        ("effect_pct_ci_high", "ci_high"),
    ]:
        table[column] = pct(table[source])
    columns = ["n", "estimate", "se", "t", "df", "p", "ci_low", "ci_high", "q"]
    return table[columns + ["effect_pct", "effect_pct_ci_low", "effect_pct_ci_high", "mrna_coef", "status"]]


def fit_protein_age(protein, mrna, samples, covariates, pca=None):
    """mRNA-adjusted protein age effect: fit_age of protein with the gene's own mRNA."""
    return fit_age(protein, samples, covariates, mrna=mrna, pca=pca)


def fit_layers(protein, rna, samples, covariates, pca=None):
    """Age fits of the layers a unit has: mRNA-adjusted protein, protein and mRNA, with the same protected PCs."""
    fits = {}
    if protein is not None and rna is not None:
        fits["mrna_adjusted"] = fit_protein_age(protein, rna, samples, covariates, pca)
    if protein is not None:
        fits["protein_only"] = fit_age(protein, samples, covariates, pca=pca)
    if rna is not None:
        fits["mrna_only"] = fit_age(rna, samples, covariates, pca=pca)
    return fits


def property_effect(data, prop, covariates):
    """Age effect (data.estimate) against a gene property: OLS with HC3 standard errors, effect also in percent."""
    X = sm.add_constant(data[[prop, *covariates]].astype(float))
    fit = sm.OLS(data.estimate.astype(float), X).fit(cov_type="HC3")
    low, high = fit.conf_int().loc[prop]
    estimate = fit.params[prop]
    return dict(
        n=len(data),
        estimate=estimate,
        se=fit.bse[prop],
        ci_low=low,
        ci_high=high,
        z=fit.tvalues[prop],
        p=fit.pvalues[prop],
        effect_pct=pct(estimate),
        effect_pct_ci_low=pct(low),
        effect_pct_ci_high=pct(high),
    )


def matched_camera(genes, marked, match, n_bins, name):
    """cameraPR of the marked genes against unmarked genes in the same joint quantile bins of the matching columns.

    In every bin the larger group is sampled down (fixed seed) to the size of the smaller one; returns the test row and
    the two gene lists.
    """
    bins = 0
    for column in match:
        bins = bins * n_bins + pd.qcut(genes[column].rank(method="first"), n_bins, labels=False)
    rng = np.random.default_rng(SEED)
    carriers, background = [], []
    for _, group in genes.groupby(bins):
        has_site = marked[group.index].to_numpy()
        size = min(has_site.sum(), (~has_site).sum())
        carriers += rng.choice(group.index[has_site], size, replace=False).tolist()
        background += rng.choice(group.index[~has_site], size, replace=False).tolist()
    return camera(genes.t[carriers + background], {name: carriers}).loc[name], carriers, background


def camera(t, sets):
    """Competitive set test (limma cameraPR) of a per-gene statistic; sets maps name -> gene ids."""
    with tempfile.TemporaryDirectory() as folder:
        t.rename("t").rename_axis("gene_id").reset_index().to_csv(f"{folder}/stat.tsv", sep="\t", index=False)
        rows = [(name, gene) for name, genes in sets.items() for gene in genes]
        pd.DataFrame(rows, columns=["set", "gene_id"]).to_csv(f"{folder}/sets.tsv", sep="\t", index=False)
        _run_r("camera.R", folder)
        out = pd.read_csv(f"{folder}/camera.tsv", sep="\t", dtype={"set": str}, keep_default_na=False).set_index("set")
    # signed z for combining cohorts
    out["z"] = sps.norm.isf(out.p / 2) * np.where(out.direction == "Up", 1, -1)
    return out


def causal_forest(data, classes, threads):
    """Causal forest (grf) of y on w: average partial effect overall, in and out of each gene class, and the difference.

    data has the columns y, w and x_*; classes maps a class name to 1, 0 or missing for each gene.
    """
    with tempfile.TemporaryDirectory() as folder:
        data = data.assign(**{f"g_{name}": member for name, member in classes.items()})
        data.to_csv(f"{folder}/data.tsv", sep="\t", index=False, na_rep="NA")
        _run_r("causal_forest.R", folder, threads)
        return pd.read_csv(f"{folder}/effects.tsv", sep="\t")


def batch_blup(protein, mrna, samples, covariates):
    """Remove a random effect of the protein batch (lme4 BLUP) from protein, keeping age, mRNA and covariates."""
    with tempfile.TemporaryDirectory() as folder:
        protein.to_csv(f"{folder}/protein.tsv", sep="\t", index_label="feature_id")
        mrna.loc[protein.index].to_csv(f"{folder}/rna.tsv", sep="\t", index_label="feature_id")
        design(samples, covariates).to_csv(f"{folder}/design.tsv", sep="\t", index=False)
        pd.DataFrame({"sample_id": samples.index, "batch": samples.protein_batch.to_numpy()}).to_csv(
            f"{folder}/samples.tsv", sep="\t", index=False
        )
        _run_r("batch_blup.R", folder, R_CORES)
        corrected = pd.read_csv(f"{folder}/corrected.tsv", sep="\t", index_col=0, dtype={"feature_id": str})
        status = pd.read_csv(f"{folder}/batch_audit.tsv", sep="\t", index_col=0).status
    assert status.isin(["ok", "insufficient"]).all(), "batch correction failed"
    return corrected


def fit_slopes(y, X, age, donor, focal):
    """Mixed model with a donor intercept: age slopes in the focal and the other rows (log2 per decade), difference."""
    X = X.assign(age_focal=age * focal, age_other=age * ~focal)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # fits at a zero donor variance warn
        fit = sm.MixedLM(y, X, groups=donor).fit(reml=True, method="powell")
    names = ["age_focal", "age_other"]
    slope, cov = fit.params[names], fit.cov_params().loc[names, names].to_numpy()
    diff, var = slope.iloc[0] - slope.iloc[1], cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
    return dict(
        slope_focal=slope.iloc[0],
        se_focal=np.sqrt(cov[0, 0]),
        slope_other=slope.iloc[1],
        se_other=np.sqrt(cov[1, 1]),
        estimate=diff,
        se=np.sqrt(var),
        p=2 * sps.norm.sf(abs(diff) / np.sqrt(var)),
        converged=fit.converged,
        n_donors=len(set(donor)),
    )
