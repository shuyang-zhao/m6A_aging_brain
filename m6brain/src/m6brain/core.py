"""Age prediction (m6Brain clock) and label-free scoring of a proteome with the m6Brain weights; see docs/method.md."""

from __future__ import annotations

import json
import re
import warnings
from functools import cache
from importlib.resources import files
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SPECIES = ("mouse", "human")
ENSEMBL = re.compile(r"^ENS([A-Z]{3})?G\d{6,}", re.IGNORECASE)
MIN_MATCH = 0.5  # match rate below which the species / id type is questioned


@cache
def _table(species):
    with (files("m6brain") / "data" / f"{species}_weights.tsv").open() as f:
        return pd.read_csv(f, sep="\t")


def weights_version():
    """Date, git commit and size of the weight tables shipped with this package."""
    return json.loads((files("m6brain") / "data" / "version.json").read_text())


def weights(species):
    """Weight table of a species: gene_id, symbol, clock_weight, prediction, observed_z, m6A_sites."""
    return _table(species).copy()


def gene(symbol_or_id, species):
    """Row of one gene: clock_weight, prediction (label-free weight, predicted age z), observed_z and m6A_sites."""
    table = _table(species)
    query = str(symbol_or_id).strip()
    hit = table[(table.symbol == query) | (table.gene_id == re.sub(r"\.\d+$", "", query).upper())]
    if hit.empty:  # symbols are case-sensitive in predict_age(); a single-gene lookup is forgiving
        hit = table[table.symbol.str.upper() == query.upper()]
    if hit.empty:
        raise KeyError(f"{query!r} has no weight in the {species} table")
    return hit.iloc[0]


def _read(table):
    """Gene ids (first column or index) and the numeric sample matrix of a DataFrame or a CSV/TSV path."""
    if isinstance(table, (str, Path)):
        tab = str(table).lower().removesuffix(".gz").endswith((".tsv", ".txt", ".tab"))
        table = pd.read_csv(table, sep="\t" if tab else ",")
    if isinstance(table.index, pd.RangeIndex):
        ids, table = table.iloc[:, 0], table.iloc[:, 1:]
    else:
        ids = table.index.to_series()
    values = table.astype(float)
    values.columns = values.columns.astype(str)
    return ids.astype(str).str.strip().reset_index(drop=True), values.reset_index(drop=True)


def _keys(ids, kind):
    """Ids as in the weight table: first entry of 'A;B' protein groups; Ensembl ids without version, upper-case."""
    ids = ids.str.split(";").str[0].str.strip()
    return ids.str.replace(r"\.\d+$", "", regex=True).str.upper() if kind == "ensembl" else ids


def _rate(ids, species, kind):
    table = _table(species)
    return _keys(ids, kind).isin(table.gene_id if kind == "ensembl" else table.symbol).mean()


def _check_species(ids, species, kind):
    rate = _rate(ids, species, kind)
    if rate >= MIN_MATCH:
        return
    other = SPECIES[SPECIES.index(species) - 1]
    msg = f"only {rate:.0%} of the {kind} ids match the {species} table"
    if kind == "symbol":
        table = _table(species)
        folded = ids.str.split(";").str[0].str.strip().str.upper().isin(table.symbol.str.upper()).mean()
        if folded > rate:
            msg += f" ({folded:.0%} ignoring case; symbols are case-sensitive: mouse Gnai3, human GNAI3)"
    if _rate(ids, other, kind) > rate:
        msg += f"; {_rate(ids, other, kind):.0%} match {other}: is species={species!r} wrong?"
    warnings.warn(msg, stacklevel=3)


def _by_sample(x, samples):
    if isinstance(x, pd.Series):
        x = x.set_axis(x.index.astype(str)).reindex(samples)
    return np.asarray(x, float)


def _pearson(x, y):
    ok = ~(np.isnan(x) | np.isnan(y))
    r = stats.pearsonr(x[ok], y[ok])
    return {"pearson_r": float(r.statistic), "pearson_p": float(r.pvalue)}


def _prepare(table, species, ids, age, column):
    """Proteins of the weighted genes standardised across the samples (missing stays NaN), their weights, the ages."""
    gene_ids, values = _read(table)
    is_age = gene_ids.str.lower() == "age"
    if is_age.any() and age is None:
        age = values[is_age].iloc[0]
    gene_ids, values = gene_ids[~is_age], values[~is_age]
    if values.shape[1] < 3:
        raise ValueError("at least 3 samples are needed")

    if ids == "auto":
        ids = "ensembl" if gene_ids.str.match(ENSEMBL).mean() > 0.5 else "symbol"
    _check_species(gene_ids, species, ids)
    keys = _keys(gene_ids, ids)
    weight = keys.map(_table(species).set_index("gene_id" if ids == "ensembl" else "symbol")[column]).replace(0, np.nan)
    keep = weight.notna() & ~keys.duplicated()  # a gene in several rows counts once (the first)
    if not keep.any():
        raise ValueError(f"no {ids} id matches a weighted gene of the {species} table")

    x = values[keep]
    z = x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1), axis=0)
    attrs = {"species": species, "ids": ids, "genes_matched": int(keep.sum()), "genes_input": len(gene_ids)}
    return z, weight[keep], values.columns, None if age is None else _by_sample(age, values.columns), attrs


def predict_age(table, species, age=None, controls=None, ids="auto"):
    """m6Brain clock: predicted age (months in mouse, years in human) of every sample of a proteome.

    table: DataFrame (gene ids in the first column, or in a non-default index) or a CSV/TSV path (.tsv/.txt = tab) of
        log abundances, genes x samples; a row whose id is "age" holds the ages of the samples.
    species: "mouse" or "human". ids: "auto", "symbol" or "ensembl". age: per-sample ages (Series by sample name or a
        list in column order); overrides an "age" row. controls: names of the samples on which the age
        acceleration is fitted (default all samples).
    Predicted age = intercept + sum of weight x protein, each protein standardised across the samples of the table
    (missing = 0), so the mean over the samples is the training mean age (exactly without missing values):
    differences between samples carry the information, the level of the table does not.
    Returns a DataFrame (sample, predicted_age, n_genes_matched = clock proteins quantified in the sample; with ages
    also age and age_acceleration = residual of predicted_age on age in a linear fit over the controls). With ages,
    result.attrs["pearson_r"] and ["pearson_p"] hold the Pearson correlation of predicted_age and age.
    """
    z, w, samples, age, attrs = _prepare(table, species, ids, age, "clock_weight")
    version = weights_version()[species]
    result = pd.DataFrame(
        {
            "sample": samples,
            "predicted_age": (z.fillna(0).T @ w).to_numpy() + version["intercept"],
            "n_genes_matched": z.notna().sum().to_numpy(),
        }
    )
    result.attrs |= attrs | {"clock_genes": version["clock_genes"], "unit": version["unit"]}
    if age is not None:
        result["age"] = age
        result.attrs |= _pearson(result.predicted_age.to_numpy(), age)
        fit = ~np.isnan(age)
        if controls is not None:
            fit &= samples.isin(controls)
        slope, offset = np.polyfit(age[fit], result.predicted_age[fit], 1)
        result["age_acceleration"] = result.predicted_age - (offset + slope * age)
    return result


def score(table, species, ids="auto", age=None):
    """Label-free score of every sample of a proteome (genes x samples, log abundances).

    The weights are the age z of each gene predicted from its transcript (column prediction), not fitted on any
    proteome. Arguments as in predict_age.
    Returns a DataFrame (sample, score, n_genes_matched; age if given). With ages, result.attrs["pearson_r"] and
    ["pearson_p"] hold the Pearson correlation of score and age.
    """
    z, w, samples, age, attrs = _prepare(table, species, ids, age, "prediction")
    raw = z.fillna(0).T @ w  # score standardised across the samples
    result = pd.DataFrame(
        {
            "sample": samples,
            "score": ((raw - raw.mean()) / raw.std()).to_numpy(),
            "n_genes_matched": z.notna().sum().to_numpy(),
        }
    )
    result.attrs |= attrs
    if age is not None:
        result["age"] = age
        result.attrs |= _pearson(result.score.to_numpy(), age)
    return result
