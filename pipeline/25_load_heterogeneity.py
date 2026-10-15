# Where the m6A load effect on the age effect concentrates: one causal forest per unit and layer (outcome the gene age
# effect, treatment the m6A load, covariates those of the load model and the protein lifetime) gives the average
# partial effect inside and outside three gene classes (neuron-enriched: neuron share 0.5 or more; synaptic: GO
# synapse; long-lived: protein lifetime above the median) and the difference; cohorts combined by random effects, BH
# over the differences
# reads scripts 10 and 09 and the unit levels; writes output/25_load_heterogeneity/
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import model
import utils

OUT = utils.OUT / "25_load_heterogeneity"
FAMILIES = [("mouse_protein", "protein_only"), ("mouse_protein", "mrna_adjusted"), ("human_protein", "protein_only")]
COVARIATES = [*utils.LOAD_COVARIATES, "log2_protein_lifetime"]  # lifetime may be missing
WORKERS, THREADS = 14, 4

features = {
    species: utils.read(utils.OUT / "09_transcript_features" / f"{species}_features.tsv.gz").set_index("gene_id")
    for species in ("mouse", "human")
}
units = {(context, layer): listed.split() for context, layer, listed in utils.CONTEXTS if (context, layer) in FAMILIES}
effects = {
    layer: utils.read(
        utils.OUT / "10_protein_age_effects" / f"{layer}.tsv.gz", usecols=["unit", "gene_id", "estimate", "status"]
    ).query("status == 'ok'")
    for layer in {layer for _, layer in FAMILIES}
}
cohort_of = utils.read(utils.PREPARED / "units.tsv", index_col="unit").dataset


def unit_effects(unit, layer):
    # m6A load as in script 21 (human: peaks of the brain m6A-seq map) and the three gene classes of the unit's genes
    table = effects[layer][effects[layer].unit == unit]
    species = "mouse" if table.gene_id.iloc[0].startswith("ENSMUSG") else "human"
    protein, rna, _, _ = utils.load_unit(unit)
    level = (protein if layer == "protein_only" else rna).mean(axis=1).rename("level")
    genes = table.set_index("gene_id")[["estimate"]].join(features[species]).join(level)
    genes = genes.dropna(subset=["estimate", *utils.LOAD_COVARIATES])
    count = genes["n_peaks_m6A_brain" if species == "human" else "n_sites_m6A"]
    lifetime = genes.log2_protein_lifetime
    classes = {
        "neuron-enriched": (genes.ct_frac_neuron >= 0.5).astype(float),
        "synaptic": genes.synapse.astype(float),
        "long-lived": (lifetime > lifetime.median()).astype(float).where(lifetime.notna()),
    }
    data = pd.DataFrame({"y": genes.estimate, "w": np.log1p(count)}).join(genes[COVARIATES].add_prefix("x_"))
    return model.causal_forest(data, classes, THREADS).assign(unit=unit, layer=layer, n_genes=len(genes))


# every unit and layer in parallel
jobs = sorted({(unit, layer) for (_, layer), listed in units.items() for unit in listed})
results = Parallel(n_jobs=WORKERS)(delayed(unit_effects)(unit, layer) for unit, layer in jobs)
unit_table = pd.concat(results).assign(cohort=lambda d: d.unit.map(cohort_of))
utils.write(unit_table, OUT / "unit_effects.tsv")

# contexts: units of a cohort averaged, cohorts by random effects, in percent; BH over the differences
rows = []
for (context, layer), listed in units.items():
    part = unit_table[(unit_table.layer == layer) & unit_table.unit.isin(listed)]
    for (gene_class, quantity), group in part.groupby(["class", "quantity"], sort=False):
        cohorts = group.groupby("cohort")[["estimate", "se"]].mean()
        effect = utils.random_effects(cohorts.estimate, cohorts.se)
        row = dict(context=context, layer=layer, gene_class=gene_class, quantity=quantity, n_cohorts=len(cohorts))
        row |= dict(n_negative=(cohorts.estimate < 0).sum(), z=effect["z"], p=effect["p"], i2=effect["i2"])
        row |= dict(effect_pct=utils.pct(effect["est"]), effect_pct_ci_low=utils.pct(effect["ci_low"]))
        rows.append(row | dict(effect_pct_ci_high=utils.pct(effect["ci_high"])))
combined = pd.DataFrame(rows)
difference = combined.quantity == "difference"
combined.loc[difference, "q"] = combined[difference].groupby(["context", "layer"]).p.transform(utils.bh)
utils.write(combined, OUT / "combined_effects.tsv")
shown = ["context", "layer", "gene_class", "quantity", "n_negative", "n_cohorts", "effect_pct", "z", "p", "q"]
print(combined[shown].round(4).to_string(index=False))
