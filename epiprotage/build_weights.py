"""Build the weight tables of EpiProtAge (src/epiprotage/data/) from the clock run of the analysis pipeline.

Run from anywhere after pipeline/50 and 52 have been (re)run:  env/bin/python epiprotage/build_weights.py
weight = TabPFN prediction of the age z of the gene (model of all cohorts, pipeline/52_proteome_clock.py),
observed_z = combined age z of the gene over the cohorts (the target); m6A_sites = m6A sites of the gene (script 50).
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
DATA = Path(__file__).resolve().parent / "src" / "epiprotage" / "data"


def git(*args):
    return subprocess.run(["git", "-C", ROOT, *args], capture_output=True, text=True, check=True).stdout.strip()


version = {"date": datetime.now().astimezone().date().isoformat(), "git_commit": git("rev-parse", "--short=8", "HEAD")}

DATA.mkdir(parents=True, exist_ok=True)
for species in ["mouse", "human"]:
    run = OUT / "52_proteome_clock" / "runs" / f"{species}_protein__all.parquet"
    runs = pd.read_parquet(run)
    load = pd.read_parquet(
        OUT / "50_transcript_ageing_model" / "inputs" / f"{species}_protein.parquet", columns=["m6A_all"]
    )
    table = pd.DataFrame(
        {
            "symbol": runs.symbol,
            "weight": runs.model,
            "observed_z": runs.target,
            "m6A_sites": np.rint(np.expm1(load.m6A_all.loc[runs.index])).astype(int),  # m6A_all = log1p(count)
        }
    )
    path = DATA / f"{species}_weights.tsv"
    table.to_csv(path, sep="\t", float_format="%.8g")
    version[species] = {
        "genes": len(table),
        "weights_computed": datetime.fromtimestamp(run.stat().st_mtime).astimezone().date().isoformat(),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest()[:16],
    }
    print(species, table.shape, version[species]["weights_computed"])

(DATA / "version.json").write_text(json.dumps(version, indent=1) + "\n")
print(version)
