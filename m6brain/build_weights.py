"""Build the weight tables of m6Brain (src/m6brain/data/) from the clock runs of the analysis pipeline.

Run from anywhere after pipeline/50, 52 and 55 have been (re)run:  env/bin/python m6brain/build_weights.py
clock_weight = coefficient of the m6Brain clock, an elastic net of age on the standardised proteins guided by m6Brain
(pipeline/55_m6brain_clock.py, fitted on all training cohorts; 0 for genes outside the clock); its intercept is in
version.json. prediction = label-free weight, the TabPFN prediction of the age z of the gene (pipeline/52); observed_z =
combined age z of the gene over the cohorts (the target); m6A_sites = m6A sites of the gene (script 50).
"""

import hashlib
import json
import subprocess
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
DATA = Path(__file__).resolve().parent / "src" / "m6brain" / "data"
UNIT = {"mouse": "months", "human": "years"}


def git(*args):
    return subprocess.run(["git", "-C", ROOT, *args], capture_output=True, text=True, check=True).stdout.strip()


version = {"date": datetime.now().astimezone().date().isoformat(), "git_commit": git("rev-parse", "--short=8", "HEAD")}

DATA.mkdir(parents=True, exist_ok=True)
for species in ["mouse", "human"]:
    runs = pd.read_parquet(OUT / "52_proteome_clock" / "runs" / f"{species}_protein__all.parquet")
    load = pd.read_parquet(
        OUT / "50_transcript_ageing_model" / "inputs" / f"{species}_protein.parquet", columns=["m6A_all"]
    )
    file = OUT / "55_m6brain_clock" / "weights" / f"{species}_protein__all__m6brain_clock.tsv"
    clock = pd.read_csv(file, sep="\t", index_col="gene_id").weight
    intercept = clock.pop("(intercept)")
    table = pd.DataFrame(
        {
            "symbol": runs.symbol,
            "clock_weight": clock.reindex(runs.index).fillna(0),
            "prediction": runs.model,
            "observed_z": runs.target,
            "m6A_sites": np.rint(np.expm1(load.m6A_all.loc[runs.index])).astype(int),  # m6A_all = log1p(count)
        }
    )
    path = DATA / f"{species}_weights.tsv"
    table.to_csv(path, sep="\t", float_format="%.10g")
    version[species] = {
        "genes": len(table),
        "clock_genes": int((table.clock_weight != 0).sum()),
        "intercept": float(intercept),
        "unit": UNIT[species],
        "weights_computed": datetime.fromtimestamp(file.stat().st_mtime).astimezone().date().isoformat(),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest()[:16],
    }
    print(species, table.shape, version[species]["clock_genes"], len(clock))

(DATA / "version.json").write_text(json.dumps(version, indent=1) + "\n")
print(version)
