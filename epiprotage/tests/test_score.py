"""Scoring equals the reference implementation (pipeline/52_proteome_clock.py donor_scores) and the pipeline's rho."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import epiprotage

HERE = Path(__file__).resolve().parent
EXAMPLE = HERE.parent / "examples" / "data" / "scifo2025_whole_brain.csv"
COHORTS = HERE.parents[1] / "output" / "52_proteome_clock" / "cohorts.tsv"


def standardised(x):
    return (x - x.mean()) / x.std()


def reference(protein, weight):
    """donor_scores of the pipeline for one unit: standardised proteins (missing = 0), weighted sum, standardised."""
    z = protein.reindex(weight.index).T.pipe(standardised).fillna(0)
    return (z @ weight).pipe(standardised)


@pytest.fixture(scope="module")
def mouse():
    return epiprotage.weights("mouse").set_index("gene_id")


@pytest.fixture(scope="module")
def synthetic(mouse):
    rng = np.random.default_rng(1)
    genes = list(mouse.index[:400]) + ["ENSMUSG00000099999"]  # one gene without a weight
    x = pd.DataFrame(rng.normal(size=(len(genes), 15)), index=genes, columns=[f"s{i}" for i in range(15)])
    x = x.mask(rng.random(x.shape) < 0.1)
    x.iloc[3] = 5.0  # constant protein
    x.iloc[4, 1:] = np.nan  # protein quantified once
    x.index.name = "gene_id"
    return x


def test_score_equals_reference(synthetic, mouse):
    result = epiprotage.score(synthetic.reset_index(), "mouse")
    expected = reference(synthetic, mouse.weight)
    assert result["sample"].tolist() == list(synthetic.columns)
    np.testing.assert_allclose(result.score, expected, atol=1e-12)
    assert result.score.mean() == pytest.approx(0, abs=1e-12) and result.score.std() == pytest.approx(1)


def test_score_symbols_equal_ids(synthetic, mouse):
    by_symbol = synthetic.rename(index=mouse.symbol)
    np.testing.assert_allclose(epiprotage.score(by_symbol, "mouse").score, epiprotage.score(synthetic, "mouse").score)


def test_n_genes_matched(synthetic):
    result = epiprotage.score(synthetic, "mouse")
    quantified = (
        synthetic.iloc[:400].drop(synthetic.index[[3, 4]]).notna().sum()
    )  # without the constant and the single-value protein
    assert result.n_genes_matched.tolist() == quantified.tolist()


def test_example_equals_pipeline_rho():
    result = epiprotage.score(EXAMPLE, "mouse")
    assert result.attrs["spearman_rho"] == pytest.approx(stats.spearmanr(result.score, result.age).statistic)
    if not COHORTS.exists():
        pytest.skip("analysis outputs not available")
    row = pd.read_csv(COHORTS, sep="\t").query(
        "`set` == 'external' and weights == 'model' and cohort == 'scifo2025_pxd058684'"
    )
    assert len(result) == row.n.item()
    assert result.attrs["spearman_rho"] == pytest.approx(row.rho.item(), abs=1e-6)


def test_age_argument(synthetic):
    age = pd.Series(np.arange(15.0), index=synthetic.columns)
    rho = epiprotage.score(synthetic, "mouse", age=age).attrs["spearman_rho"]
    assert epiprotage.score(synthetic, "mouse", age=age.tolist()).attrs["spearman_rho"] == rho
    with_row = pd.concat([synthetic, age.rename("age").to_frame().T])
    assert epiprotage.score(with_row, "mouse").attrs["spearman_rho"] == rho
