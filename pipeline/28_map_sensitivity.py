# Sensitivity of the m6A results to the site map that defines the marked genes: the matched cameraPR signature of
# script 20 and the OLS/HC3 load effect of script 21 (model without protein lifetime) with the genes marked by one map
# or the union; human: brain peaks vs cell-line (GLORI) sites, mouse: each of the four maps vs the union
# reads scripts 10 and 09, the site atlas and the unit levels; writes output/28_map_sensitivity/
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import model
import utils

OUT = utils.OUT / "28_map_sensitivity"
WORKERS = 16
N_BINS, MATCH = 10, ["level", "log_tx_len"]
CONTEXTS = ["human_striatum", "human_cortex", "human_hippocampus", "human_cerebellum", "human_neurons", "human_protein"]
CONTEXTS += ["mouse_mrna", "mouse_neurons", "mouse_protein"]
MAPS = {
    "human": {"brain peaks": utils.BRAIN_M6A, "cell line": "m6a_hek293t_hela"},
    "mouse": {
        "hippocampus cyto": "m6a_hippocampus_cyto",
        "hippocampus syn": "m6a_hippocampus_syn",
        "cortex scDART": "m6a_cortex_scdart_tegowski2024",
        "neuron DART": "m6a_neuron_dart_tegowski2024",
    },
}

features = {
    species: utils.read(utils.OUT / "09_transcript_features" / f"{species}_features.tsv.gz").set_index("gene_id")
    for species in ("mouse", "human")
}
cohort_of = utils.read(utils.PREPARED / "units.tsv", index_col="unit").dataset
sites = utils.read(utils.RESOURCES / "site_atlas" / "sites.parquet")
sites = sites[(sites.modification == "m6A") & sites.gene_id.notna() & sites.tx_offset.notna() & sites.ref_ok]
sites = sites.drop_duplicates(["map", "chrom", "start", "end", "strand"])
per_gene = sites.groupby(["map", "gene_id"]).size()

# sites per gene and definition; the human union marks a gene carried by either map, the mouse union is the atlas count
counts = {}
for species, maps in MAPS.items():
    counts[species] = {name: per_gene[key].reindex(features[species].index).fillna(0) for name, key in maps.items()}
counts["human"]["union"] = ((counts["human"]["brain peaks"] > 0) | (counts["human"]["cell line"] > 0)).astype(int)
counts["mouse"]["union"] = features["mouse"].n_sites_m6A
LOAD_DEFINITIONS = {"human": ["brain peaks", "cell line"], "mouse": [*MAPS["mouse"], "union"]}

families = [(context, layer, units.split()) for context, layer, units in utils.CONTEXTS if context in CONTEXTS]
jobs = sorted({(unit, layer) for _, layer, units in families for unit in units})


def matched_test(genes):
    # as script 20: joint deciles of level and length
    test = model.matched_camera(genes, genes.count_sites > 0, MATCH, N_BINS, "m6A")[0]
    return dict(n_set=test.n_genes, z=test.z, p=test.p)


def unit_tests(unit, layer):
    table = utils.read(
        utils.OUT / "10_protein_age_effects" / f"{layer}.tsv.gz", usecols=["unit", "gene_id", "estimate", "t", "status"]
    )
    table = table[(table.unit == unit) & (table.status == "ok")].set_index("gene_id")
    protein, rna, samples, _ = utils.load_unit(unit)
    species = "mouse" if table.index[0].startswith("ENSMUSG") else "human"
    # signature: mean level of the layer; load: mean protein level, otherwise mean mRNA level (as script 21)
    level = (rna if layer == "mrna_only" else protein).mean(axis=1)
    load_level = (protein if layer == "protein_only" else rna).mean(axis=1)
    genes = table.join(features[species], how="inner").assign(level=level)
    signature, load = [], []
    for name, count in counts[species].items():
        data = genes.assign(count_sites=count.reindex(genes.index).to_numpy())
        key = dict(unit=unit, layer=layer, species=species, definition=name, n_samples=len(samples))
        signature.append(key | dict(n_marked=(data.count_sites > 0).sum()) | matched_test(data))
        if name in LOAD_DEFINITIONS[species]:
            data = data.assign(level=load_level.reindex(data.index))
            data = data.dropna(subset=["estimate", *utils.LOAD_COVARIATES])
            data = data.assign(load=np.log1p(data.count_sites))
            fit = model.property_effect(data, "load", utils.LOAD_COVARIATES)
            load.append(key | dict(n_marked=(data.load > 0).sum()) | fit)
    return signature, load


results = Parallel(n_jobs=WORKERS)(delayed(unit_tests)(unit, layer) for unit, layer in jobs)
signature = pd.DataFrame([row for rows, _ in results for row in rows]).assign(cohort=lambda d: d.unit.map(cohort_of))
load = pd.DataFrame([row for _, rows in results for row in rows]).assign(cohort=lambda d: d.unit.map(cohort_of))
utils.write(signature, OUT / "unit_signatures.tsv")
utils.write(load, OUT / "unit_load_effects.tsv")

# contexts of utils.CONTEXTS: signature z of a cohort's units averaged, Stouffer with the square root of the samples;
# load coefficients of a cohort averaged, cohorts combined by random effects
sig_rows, load_rows = [], []
for context, layer, units in families:
    family = signature[(signature.layer == layer) & signature.unit.isin(units)]
    for definition, group in family.groupby("definition", sort=False):
        cohorts = group.groupby("cohort").agg(z=("z", "mean"), n=("n_samples", "mean"))
        pooled = utils.stouffer(cohorts.z, np.sqrt(cohorts.n))
        row = dict(context=context, layer=layer, definition=definition, n_cohorts=len(cohorts))
        sig_rows.append(row | dict(n_marked=int(group.n_marked.mean()), n_set=int(group.n_set.mean())) | pooled)
    family = load[(load.layer == layer) & load.unit.isin(units)]
    for definition, group in family.groupby("definition", sort=False):
        cohorts = group.groupby("cohort")[["estimate", "se"]].mean()
        effect = utils.random_effects(cohorts.estimate, cohorts.se)
        row = dict(context=context, layer=layer, definition=definition, n_cohorts=len(cohorts))
        row |= dict(n_marked=int(group.n_marked.mean()), effect_pct=utils.pct(effect["est"]))
        row |= dict(ci_low_pct=utils.pct(effect["ci_low"]), ci_high_pct=utils.pct(effect["ci_high"]))
        load_rows.append(row | dict(z=effect["z"], p=effect["p"], i2=effect["i2"]))
combined_signature, combined_load = pd.DataFrame(sig_rows), pd.DataFrame(load_rows)
# BH over the definitions of a context and layer
for table in (combined_signature, combined_load):
    table["q"] = table.groupby(["context", "layer"]).p.transform(utils.bh)
utils.write(combined_signature, OUT / "combined_signatures.tsv")
utils.write(combined_load, OUT / "combined_load_effects.tsv")
print(combined_signature.round(4).to_string(index=False))
print(combined_load.round(4).to_string(index=False))
