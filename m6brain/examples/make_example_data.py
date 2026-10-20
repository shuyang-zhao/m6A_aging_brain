"""Export the example proteome (Scifo et al. 2025, PXD058684, mouse whole brain, 45 samples) from the analysis pipeline.

Run from the repository root:  env/bin/python m6brain/examples/make_example_data.py
"""

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
unit = ROOT / "output" / "prepared" / "scifo2025_whole_brain"

protein = pd.read_parquet(unit / "protein.parquet")  # Ensembl gene ids x samples, log intensities
samples = pd.read_csv(unit / "samples.tsv", sep="\t", dtype={"sample_id": str}).set_index("sample_id")
symbol = (
    pd.read_csv(
        ROOT / "output" / "09_transcript_features" / "mouse_features.tsv.gz", sep="\t", usecols=["gene_id", "symbol"]
    )
    .set_index("gene_id")
    .symbol
)

table = protein.rename(index=symbol)
table = table[table.index.isin(symbol)].astype(float)
table.index.name = "symbol"
table.loc["age"] = samples.age_months.reindex(table.columns)
table.to_csv(Path(__file__).parent / "data" / "scifo2025_whole_brain.csv", float_format="%.5g")
print(table.shape)
