"""Scoring of a proteome with the m6Brain weights (predicted age z per gene); see docs/method.md."""

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
    """Weight table of a species: gene_id, symbol, weight (predicted age z), observed_z, m6A_sites."""
    return _table(species).copy()


def gene(symbol_or_id, species):
    """Row of one gene: weight (predicted age z), observed_z (combined age z of the cohorts) and m6A_sites."""
    table = _table(species)
    query = str(symbol_or_id).strip()
    hit = table[(table.symbol == query) | (table.gene_id == re.sub(r"\.\d+$", "", query).upper())]
    if hit.empty:  # symbols are case-sensitive in score(); a single-gene lookup is forgiving
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


def score(table, species, ids="auto", age=None):
    """Clock score of every sample of a proteome (genes x samples, log abundances).

    table: DataFrame (gene ids in the first column, or in a non-default index) or a CSV/TSV path (.tsv/.txt = tab);
        a row whose id is "age" holds the ages of the samples.
    species: "mouse" or "human". ids: "auto", "symbol" or "ensembl". age: per-sample ages (Series by sample name or a
        list in column order); overrides an "age" row.
    Returns a DataFrame (sample, score, n_genes_matched; age if given). With ages, result.attrs["spearman_rho"] and
    ["spearman_p"] hold the Spearman correlation of score and age.
    """
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
    weight = keys.map(_table(species).set_index("gene_id" if ids == "ensembl" else "symbol").weight)
    keep = weight.notna() & ~keys.duplicated()  # a gene in several rows counts once (the first)
    if not keep.any():
        raise ValueError(f"no {ids} id matches the {species} table")

    # each protein standardised across the samples (missing -> 0), weighted sum, score standardised across samples
    x, w = values[keep], weight[keep]
    z = x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1), axis=0)
    raw = z.fillna(0).T @ w
    result = pd.DataFrame(
        {
            "sample": values.columns,
            "score": ((raw - raw.mean()) / raw.std()).to_numpy(),
            "n_genes_matched": z.notna().sum().to_numpy(),
        }
    )
    result.attrs |= {"species": species, "ids": ids, "genes_matched": int(keep.sum()), "genes_input": len(gene_ids)}
    if age is not None:
        if isinstance(age, pd.Series):
            age = age.set_axis(age.index.astype(str)).reindex(values.columns)
        result["age"] = np.asarray(age, float)
        rho = stats.spearmanr(result.score, result.age, nan_policy="omit")
        result.attrs |= {"spearman_rho": float(rho.statistic), "spearman_p": float(rho.pvalue)}
    return result
