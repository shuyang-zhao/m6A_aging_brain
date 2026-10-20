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


def matrix(index, n=8):
    x = np.random.default_rng(0).normal(size=(len(index), n))
    return pd.DataFrame(x, index=pd.Index(index, name="id"), columns=[f"s{i}" for i in range(n)])


def test_ids_auto_detected():
    t = m6brain.weights("mouse").iloc[:200]
    by_id = m6brain.score(matrix(t.gene_id), "mouse")
    by_symbol = m6brain.score(matrix(t.symbol), "mouse")
    assert by_id.attrs["ids"] == "ensembl" and by_symbol.attrs["ids"] == "symbol"
    np.testing.assert_allclose(by_id.score, by_symbol.score)


def test_versions_and_protein_groups_stripped():
    t = m6brain.weights("human").iloc[:200]
    plain = m6brain.score(matrix(t.gene_id), "human")
    versioned = m6brain.score(matrix(t.gene_id.str.lower() + ".7"), "human")
    grouped = m6brain.score(matrix(t.symbol + ";" + t.symbol.shift(-1).fillna("X")), "human")
    np.testing.assert_allclose(versioned.score, plain.score)
    np.testing.assert_allclose(grouped.score, m6brain.score(matrix(t.symbol), "human").score)


def test_wrong_species_warnings():
    mouse, human = (m6brain.weights(s) for s in ("mouse", "human"))
    with pytest.warns(UserWarning, match="match mouse"):  # 9% human symbols among mouse symbols
        m6brain.score(matrix([*mouse.symbol[:300], *human.symbol[:30]]), "human")
    with pytest.warns(UserWarning, match="match human"):  # Ensembl ids
        m6brain.score(matrix([*human.gene_id[:300], *mouse.gene_id[:30]]), "mouse")
    with pytest.warns(UserWarning, match="ignoring case"):  # symbols are case-sensitive
        m6brain.score(matrix([*mouse.symbol[:30], *mouse.symbol[30:300].str.upper()]), "mouse")
    with pytest.warns(UserWarning), pytest.raises(ValueError, match="no symbol id matches"):
        m6brain.score(matrix(mouse.symbol[:300]), "human")


def test_errors():
    with pytest.raises(ValueError, match="3 samples"):
        m6brain.score(matrix(["Gnai3", "Apoh"], n=2), "mouse")


def test_duplicate_gene_uses_first_row():
    x = matrix(m6brain.weights("mouse").symbol.iloc[:100])
    result = m6brain.score(pd.concat([x, x.iloc[:1]]), "mouse")
    np.testing.assert_allclose(result.score, m6brain.score(x, "mouse").score)


def test_gene_lookup():
    row = m6brain.gene("Gnai3", "mouse")
    assert row.gene_id == "ENSMUSG00000000001" and row.m6A_sites == 12
    assert m6brain.gene("ENSMUSG00000000001.4", "mouse").symbol == "Gnai3"
    assert m6brain.gene("gnai3", "mouse").weight == row.weight
    assert m6brain.gene("DPM1", "human").gene_id == "ENSG00000000419"
    with pytest.raises(KeyError):
        m6brain.gene("NotAGene", "human")


def test_weights_tables():
    for species, n in [("mouse", 10606), ("human", 9990)]:
        t = m6brain.weights(species)
        assert (
            list(t) == ["gene_id", "symbol", "weight", "observed_z", "m6A_sites"]
            and len(t) == n
            and t.gene_id.is_unique
        )
    assert set(m6brain.weights_version()) >= {"date", "git_commit"}


def test_cli(tmp_path, capsys):
    out = tmp_path / "scores.csv"
    main(["score", str(EXAMPLE), "--species", "mouse", "--out", str(out)])
    scores = pd.read_csv(out)
    assert list(scores) == ["sample", "score", "n_genes_matched", "age"] and len(scores) == 45
    assert "Spearman rho with age = 0.753" in capsys.readouterr().err
    expected = m6brain.score(EXAMPLE, "mouse").score
    np.testing.assert_allclose(scores.score, expected, atol=1e-6)


def test_cli_entry_point(tmp_path):
    exe = Path(sys.executable).parent / "m6brain"
    if not exe.exists():
        pytest.skip("package not installed with its script")
    done = subprocess.run([exe, "score", EXAMPLE, "--species", "mouse"], capture_output=True, text=True, check=True)
    assert done.stdout.splitlines()[0] == "sample,score,n_genes_matched,age"
