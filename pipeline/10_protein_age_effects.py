# Age effects per gene and unit: mRNA-adjusted protein, protein-only and mRNA-only fits
# reads the prepared units; writes output/10_protein_age_effects/
import pandas as pd
from joblib import Parallel, delayed

import model
import utils

LAYERS = ["mrna_adjusted", "protein_only", "mrna_only"]
WORKERS = 40
OUT = utils.OUT / "10_protein_age_effects"


def fit_unit(unit):
    protein, rna, samples, design = utils.load_unit(unit)
    # protected PCs of the mRNA residuals, shared by the layers of the unit
    k = design["hidden_factors"]
    pcs = model.protected_pca(rna, model.design(samples, design), k) if k else None
    if unit in utils.BATCH_CORRECTED:
        protein = model.batch_blup(protein, rna, samples, design)
    return unit, model.fit_layers(protein, rna, samples, design, pcs)


# all units in parallel
units = utils.read(utils.PREPARED / "units.tsv", index_col="unit").index
results = Parallel(n_jobs=WORKERS)(delayed(fit_unit)(unit) for unit in units)

# one table per layer, units stacked
symbols = pd.concat([utils.ensg_to_symbol("mouse"), utils.ensg_to_symbol("human")])
for layer in LAYERS:
    parts = [
        fits[layer].rename_axis("gene_id").reset_index().assign(unit=unit) for unit, fits in results if layer in fits
    ]
    table = pd.concat(parts)
    table.insert(0, "unit", table.pop("unit"))
    table.insert(2, "symbol", table.gene_id.map(symbols))
    assert not table.duplicated(["unit", "gene_id"]).any()
    utils.write(table, OUT / f"{layer}.tsv.gz")
    print(f"{layer}: {len(parts)} units, {(table.status == 'ok').sum()} fits, {(table.q < utils.Q).sum()} hits")
