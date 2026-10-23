"""Id detection, species warnings, single-gene lookup and the command line."""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import m6brain
from m6brain.cli import main

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "data" / "scifo2025_whole_brain.csv"
MODES = [(m6brain.score, "prediction"), (m6brain.predict_age, "clock_weight")]  # function, weight column it uses


def matrix(index, n=8):
    x = np.random.default_rng(0).normal(size=(len(index), n))
    return pd.DataFrame(x, index=pd.Index(index, name="id"), columns=[f"s{i}" for i in range(n)])


def weighted(species, column):
    t = m6brain.weights(species)
    return t[t[column] != 0].iloc[:200]


@pytest.mark.parametrize(("fn", "column"), MODES)
def test_ids_auto_detected(fn, column):
    t = weighted("mouse", column)
    by_id, by_symbol = fn(matrix(t.gene_id), "mouse"), fn(matrix(t.symbol), "mouse")
    assert by_id.attrs["ids"] == "ensembl" and by_symbol.attrs["ids"] == "symbol"
    np.testing.assert_allclose(by_id.iloc[:, 1], by_symbol.iloc[:, 1])


@pytest.mark.parametrize(("fn", "column"), MODES)
def test_versions_and_protein_groups_stripped(fn, column):
    t = weighted("human", column)
    plain = fn(matrix(t.gene_id), "human")
    versioned = fn(matrix(t.gene_id.str.lower() + ".7"), "human")
    grouped = fn(matrix(t.symbol + ";" + t.symbol.shift(-1).fillna("X")), "human")
    np.testing.assert_allclose(versioned.iloc[:, 1], plain.iloc[:, 1])
    np.testing.assert_allclose(grouped.iloc[:, 1], fn(matrix(t.symbol), "human").iloc[:, 1])


def test_wrong_species_warnings():
    mouse, human = (m6brain.weights(s) for s in ("mouse", "human"))
    with pytest.warns(UserWarning, match="match mouse"):  # 9% human symbols among mouse symbols
        m6brain.score(matrix([*mouse.symbol[:300], *human.symbol[:30]]), "human")
    with pytest.warns(UserWarning, match="match human"):  # Ensembl ids
        m6brain.score(matrix([*human.gene_id[:300], *mouse.gene_id[:30]]), "mouse")
    with pytest.warns(UserWarning, match="ignoring case"):  # symbols are case-sensitive
        m6brain.score(matrix([*mouse.symbol[:30], *mouse.symbol[30:300].str.upper()]), "mouse")
    with pytest.warns(UserWarning), pytest.raises(ValueError, match="no symbol id matches"):
        m6brain.predict_age(matrix(mouse.symbol[:300]), "human")


def test_errors():
    for fn, _ in MODES:
        with pytest.raises(ValueError, match="3 samples"):
            fn(matrix(["Gnai3", "Apoh"], n=2), "mouse")
    with pytest.raises(ValueError, match="no symbol id matches"):  # genes in the table, none in the clock
        m6brain.predict_age(matrix(m6brain.weights("mouse").query("clock_weight == 0").symbol[:100]), "mouse")


@pytest.mark.parametrize(("fn", "column"), MODES)
def test_duplicate_gene_uses_first_row(fn, column):
    x = matrix(weighted("mouse", column).symbol)
    result = fn(pd.concat([x, x.iloc[:1]]), "mouse")
    np.testing.assert_allclose(result.iloc[:, 1], fn(x, "mouse").iloc[:, 1])


def test_gene_lookup():
    row = m6brain.gene("Gnai3", "mouse")
    assert row.gene_id == "ENSMUSG00000000001" and row.m6A_sites == 12
    assert m6brain.gene("ENSMUSG00000000001.4", "mouse").symbol == "Gnai3"
    assert m6brain.gene("gnai3", "mouse").prediction == row.prediction
    assert m6brain.gene("DPM1", "human").gene_id == "ENSG00000000419"
    with pytest.raises(KeyError):
        m6brain.gene("NotAGene", "human")


def test_weights_tables():
    version = m6brain.weights_version()
    for species, n, clock in [("mouse", 10606, 279), ("human", 9990, 421)]:
        t = m6brain.weights(species)
        columns = ["gene_id", "symbol", "clock_weight", "prediction", "observed_z", "m6A_sites"]
        assert list(t) == columns and len(t) == n and t.gene_id.is_unique
        assert (t.clock_weight != 0).sum() == clock == version[species]["clock_genes"]
        assert version[species]["unit"] == {"mouse": "months", "human": "years"}[species]
    assert version["mouse"]["intercept"] == pytest.approx(13.2165, abs=1e-4)
    assert version["human"]["intercept"] == pytest.approx(69.3683, abs=1e-4)
    assert set(version) >= {"date", "git_commit"}


def test_cli(tmp_path, capsys):
    out = tmp_path / "ages.csv"
    main(["age", str(EXAMPLE), "--species", "mouse", "--out", str(out)])
    ages = pd.read_csv(out)
    assert list(ages) == ["sample", "predicted_age", "n_genes_matched", "age", "age_acceleration"] and len(ages) == 45
    err = capsys.readouterr().err
    assert "190 of 279 clock proteins matched" in err and "Pearson r with age = 0.811" in err
    np.testing.assert_allclose(ages.predicted_age, m6brain.predict_age(EXAMPLE, "mouse").predicted_age, atol=1e-6)

    controls = ",".join(ages["sample"][ages.age < 10])
    main(["age", str(EXAMPLE), "--species", "mouse", "--controls", controls, "--out", str(out)])
    young = pd.read_csv(out).query("age < 10")
    assert young.age_acceleration.mean() == pytest.approx(0, abs=1e-6)

    main(["score", str(EXAMPLE), "--species", "mouse", "--out", str(out)])
    scores = pd.read_csv(out)
    assert list(scores) == ["sample", "score", "n_genes_matched", "age"] and len(scores) == 45
    assert "Pearson r with age = 0.742" in capsys.readouterr().err
    np.testing.assert_allclose(scores.score, m6brain.score(EXAMPLE, "mouse").score, atol=1e-6)


def test_cli_entry_point():
    exe = Path(sys.executable).parent / "m6brain"
    if not exe.exists():
        pytest.skip("package not installed with its script")
    done = subprocess.run([exe, "age", EXAMPLE, "--species", "mouse"], capture_output=True, text=True, check=True)
    assert done.stdout.splitlines()[0] == "sample,predicted_age,n_genes_matched,age,age_acceleration"
