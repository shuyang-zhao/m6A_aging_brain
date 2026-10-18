# Sensitivity of the m6A results to the site map that defines the marked genes: the matched cameraPR signature of
# script 20 and the OLS/HC3 load effect of script 21 (model without protein lifetime) with the genes marked by one map
# or the union of the mouse GLORI 3.0 hippocampus maps (cytoplasm, synaptosome); also the DRACH share of the m6A maps
# (mouse: sites in a DRACH motif; human: brain peaks containing a DRACH)
# reads scripts 10 and 09, the site atlas, the human genome and the unit levels; writes output/28_map_sensitivity/
import re

import numpy as np
import pandas as pd
import pysam
from joblib import Parallel, delayed

import model
import utils

OUT = utils.OUT / "28_map_sensitivity"
WORKERS = 16
N_BINS, MATCH = 10, ["level", "log_tx_len"]
CONTEXTS = ["mouse_mrna", "mouse_neurons", "mouse_protein"]
MAPS = {"hippocampus cyto": "m6a_hippocampus_cyto", "hippocampus syn": "m6a_hippocampus_syn"}
DRACH = re.compile(r"[AGT][AG]AC[ACT]")

features = utils.read(utils.OUT / "09_transcript_features" / "mouse_features.tsv.gz").set_index("gene_id")
cohort_of = utils.read(utils.PREPARED / "units.tsv", index_col="unit").dataset
atlas = utils.read(utils.RESOURCES / "site_atlas" / "sites.parquet")

# DRACH share of the m6A maps: mouse, sites with a DRACH motif around the site; human, peaks that contain a DRACH
genome = pysam.FastaFile(str(utils.RESOURCES / "gencode_human" / "genome.fa"))
peaks = atlas[atlas["map"] == utils.BRAIN_M6A]
has_drach = []
for chrom, start, end, strand in zip(peaks.chrom, peaks.start, peaks.end, peaks.strand, strict=True):
    text = genome.fetch(chrom, start - 1, end).upper()
    has_drach.append(bool(DRACH.search(text.translate(str.maketrans("ACGT", "TGCA"))[::-1] if strand == "-" else text)))
mouse = atlas[atlas["map"].isin(MAPS.values())]
drach = [
    dict(map=name, species="Mus musculus", n=len(rows), n_drach=(rows.motif_class == "DRACH").sum(), unit="sites")
    for name, rows in mouse.groupby("map")
]
drach.append(dict(map=utils.BRAIN_M6A, species="Homo sapiens", n=len(peaks), n_drach=sum(has_drach), unit="peaks"))
drach = pd.DataFrame(drach).assign(drach_fraction=lambda d: d.n_drach / d.n)
utils.write(drach, OUT / "map_drach.tsv")
print(drach.round(4).to_string(index=False))

sites = atlas[(atlas.modification == "m6A") & atlas.gene_id.notna() & atlas.tx_offset.notna() & atlas.ref_ok]
sites = sites[sites["map"].isin(MAPS.values())].drop_duplicates(["map", "chrom", "start", "end", "strand"])
per_gene = sites.groupby(["map", "gene_id"]).size()

# sites per gene and definition; the union is the atlas count of script 09
counts = {name: per_gene[key].reindex(features.index).fillna(0) for name, key in MAPS.items()}
counts["union"] = features.n_sites_m6A

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
    # signature: mean level of the layer; load: mean protein level, otherwise mean mRNA level (as script 21)
    level = (rna if layer == "mrna_only" else protein).mean(axis=1)
    load_level = (protein if layer == "protein_only" else rna).mean(axis=1)
    genes = table.join(features, how="inner").assign(level=level)
    signature, load = [], []
    for name, count in counts.items():
        data = genes.assign(count_sites=count.reindex(genes.index).to_numpy())
        key = dict(unit=unit, layer=layer, species="mouse", definition=name, n_samples=len(samples))
        signature.append(key | dict(n_marked=(data.count_sites > 0).sum()) | matched_test(data))
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
