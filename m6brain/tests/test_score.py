"""Clock and label-free score equal the reference implementation (donor_scores of pipeline/52 and 55) and its r."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import m6brain

HERE = Path(__file__).resolve().parent
EXAMPLE = HERE.parent / "examples" / "data" / "scifo2025_whole_brain.csv"
OUTPUT = HERE.parents[1] / "output"


def standardised(x):
    return (x - x.mean()) / x.std()


def reference(protein, weight, intercept=None):
    """donor_scores of the pipeline for one unit: standardised proteins (missing = 0), weighted sum; the label-free
    score is standardised again, the clock adds the intercept."""
    z = protein.reindex(weight.index).T.pipe(standardised).fillna(0)
    return (z @ weight).pipe(standardised) if intercept is None else z @ weight + intercept


@pytest.fixture(scope="module")
def mouse():
    return m6brain.weights("mouse").set_index("gene_id")


@pytest.fixture(scope="module")
def synthetic(mouse):
    rng = np.random.default_rng(1)
    clock = mouse.index[mouse.clock_weight != 0]
    other = mouse.index[mouse.clock_weight == 0][:300]
    genes = [*clock[:100], *other, "ENSMUSG00000099999"]  # the last gene has no weight
    x = pd.DataFrame(rng.normal(size=(len(genes), 15)), index=genes, columns=[f"s{i}" for i in range(15)])
    x = x.mask(rng.random(x.shape) < 0.1)
    x.iloc[3] = 5.0  # constant protein
    x.iloc[4, 1:] = np.nan  # protein quantified once
    x.index.name = "gene_id"
    return x


def test_score_equals_reference(synthetic, mouse):
    result = m6brain.score(synthetic.reset_index(), "mouse")
    expected = reference(synthetic, mouse.prediction)
    assert result["sample"].tolist() == list(synthetic.columns)
    np.testing.assert_allclose(result.score, expected, atol=1e-12)
    assert result.score.mean() == pytest.approx(0, abs=1e-12) and result.score.std() == pytest.approx(1)


def test_predict_age_equals_reference(synthetic, mouse):
    result = m6brain.predict_age(synthetic.reset_index(), "mouse")
    clock = mouse.clock_weight[mouse.clock_weight != 0]
    expected = reference(synthetic, clock, m6brain.weights_version()["mouse"]["intercept"])
    np.testing.assert_allclose(result.predicted_age, expected, atol=1e-12)
    assert result.attrs["genes_matched"] == 100 and result.attrs["clock_genes"] == len(clock)


def test_predicted_age_centred_on_intercept(synthetic):
    complete = synthetic.fillna(0.0) + np.random.default_rng(2).normal(size=synthetic.shape)
    mean = m6brain.predict_age(complete, "mouse").predicted_age.mean()
    assert mean == pytest.approx(m6brain.weights_version()["mouse"]["intercept"])


def test_symbols_equal_ids(synthetic, mouse):
    by_symbol = synthetic.rename(index=mouse.symbol)
    for fn in (m6brain.score, m6brain.predict_age):
        np.testing.assert_allclose(fn(by_symbol, "mouse").iloc[:, 1], fn(synthetic, "mouse").iloc[:, 1])


def test_n_genes_matched(synthetic):
    usable = synthetic.drop(synthetic.index[[3, 4]])  # without the constant and the single-value protein
    quantified = usable.iloc[:98].notna().sum()
    assert m6brain.predict_age(synthetic, "mouse").n_genes_matched.tolist() == quantified.tolist()
    both = usable.iloc[:398].notna().sum()
    assert m6brain.score(synthetic, "mouse").n_genes_matched.tolist() == both.tolist()


def test_example_equals_pipeline_r():
    clock, label_free = m6brain.predict_age(EXAMPLE, "mouse"), m6brain.score(EXAMPLE, "mouse")
    assert clock.attrs["pearson_r"] == pytest.approx(stats.pearsonr(clock.predicted_age, clock.age).statistic)
    assert label_free.attrs["pearson_r"] == pytest.approx(stats.pearsonr(label_free.score, label_free.age).statistic)
    if not OUTPUT.exists():
        pytest.skip("analysis outputs not available")
    # the example file has 5 significant digits and 4,714 of the 4,724 proteins of the unit: r agrees to ~1e-4
    modes = [(clock, "55_m6brain_clock", "m6brain_clock"), (label_free, "52_proteome_clock", "model")]
    for result, file, weights in modes:
        row = pd.read_csv(OUTPUT / file / "cohorts.tsv", sep="\t").query(
            "`set` == 'external' and weights == @weights and cohort == 'scifo2025_pxd058684'"
        )
        assert len(result) == row.n.item()
        assert result.attrs["pearson_r"] == pytest.approx(row.r.item(), abs=2e-4)


def test_age_acceleration(synthetic):
    age = pd.Series(np.arange(15.0) * 2 + 3, index=synthetic.columns)
    result = m6brain.predict_age(synthetic, "mouse", age=age)
    assert np.corrcoef(result.age_acceleration, age)[0, 1] == pytest.approx(0, abs=1e-12)
    assert result.age_acceleration.mean() == pytest.approx(0, abs=1e-12)
    controls = synthetic.columns[:8]
    fit = m6brain.predict_age(synthetic, "mouse", age=age, controls=controls)
    slope, offset = np.polyfit(age[:8], result.predicted_age[:8], 1)
    np.testing.assert_allclose(fit.age_acceleration, result.predicted_age - offset - slope * age.to_numpy(), atol=1e-10)
    assert fit.age_acceleration[:8].mean() == pytest.approx(0, abs=1e-10)
    assert fit.attrs["pearson_r"] == result.attrs["pearson_r"]


def test_age_argument(synthetic):
    age = pd.Series(np.arange(15.0), index=synthetic.columns)
    for fn in (m6brain.score, m6brain.predict_age):
        r = fn(synthetic, "mouse", age=age).attrs["pearson_r"]
        assert fn(synthetic, "mouse", age=age.tolist()).attrs["pearson_r"] == r
        with_row = pd.concat([synthetic, age.rename("age").to_frame().T])
        assert fn(with_row, "mouse").attrs["pearson_r"] == r
