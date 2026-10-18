# Age signature of the transcripts that carry each RNA modification (a site of its maps, single-base except the m7G
# miCLIP windows and METTL1-dependent regions; m6A: mouse GLORI sites, human brain m6A-seq peaks): limma cameraPR on the
# per-gene age t against unmarked transcripts matched on mean level (protein in both protein layers) and length, per
# unit, then combined across cohorts; robustness_* columns repeat the test with the matched background also built on
# the neuron share and genomic length
# reads scripts 10 and 09 and the unit levels; writes output/20_modification_signatures/ (matched_genes.tsv.gz lists
# the marked genes and the matched background of every primary test)
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import model
import utils

LAYERS = ["mrna_only", "protein_only", "mrna_adjusted"]
N_BINS = {2: 10, 4: 5}  # deciles of the two primary matching variables, quintiles of the four robustness ones
WORKERS = 40
OUT = utils.OUT / "20_modification_signatures"
MATCH = ["level", "log_tx_len"]  # mean level in the layer, transcript length
MATCH_ROBUST = [*MATCH, "ct_frac_neuron", "log_gene_len"]  # robustness: plus the neuron share and the genomic length

cohort_of = utils.read(utils.PREPARED / "units.tsv", index_col="unit").dataset
feature_columns = ["gene_id", "log_tx_len", "log_gene_len", "ct_frac_neuron"]
feature_columns += [f"n_sites_{mod}" for mod in utils.MODIFICATIONS]
features = {
    species: utils.read(
        utils.OUT / "09_transcript_features" / f"{species}_features.tsv.gz", usecols=feature_columns
    ).set_index("gene_id")
    for species in ("mouse", "human")
}
effects = {
    layer: utils.read(
        utils.OUT / "10_protein_age_effects" / f"{layer}.tsv.gz", usecols=["unit", "gene_id", "estimate", "t", "status"]
    ).query("status == 'ok'")
    for layer in LAYERS
}


def matched_test(genes, match, modification):
    # carriers have at least one site of the modification, against non-carriers matched on the quantiles of match
    marked = genes[f"n_sites_{modification}"] > 0
    test, carriers, background = model.matched_camera(genes, marked, match, N_BINS[len(match)], modification)
    gap = genes.estimate[carriers].mean() - genes.estimate[background].mean()
    stats = dict(n_set=test.n_genes, direction=test.direction, z=test.z, p=test.p, diff_pct=utils.pct(gap))
    return stats, carriers, background


def test_unit(unit, layer, table):
    # genes with a canonical transcript: age effect and t, mean level in the layer, length, neuron share, site counts
    protein, rna, _, _ = utils.load_unit(unit)
    matrix = rna if layer == "mrna_only" else protein
    species = "mouse" if table.gene_id.iloc[0].startswith("ENSMUSG") else "human"
    genes = table.set_index("gene_id").join(features[species], how="inner")
    genes["level"] = matrix.mean(axis=1).reindex(genes.index)
    assert genes.level.notna().all()
    # genes without a neuron share cannot be matched on it and stay out of the robustness test only
    with_share = genes.dropna(subset=["ct_frac_neuron"])
    rows, sets = [], []
    for modification in utils.MODIFICATIONS:
        robust, _, _ = matched_test(with_share, MATCH_ROBUST, modification)
        primary, carriers, background = matched_test(genes, MATCH, modification)
        key = dict(unit=unit, layer=layer, modification=modification)
        marked = [True] * len(carriers) + [False] * len(background)
        sets.append(pd.DataFrame(key | dict(gene_id=carriers + background, marked=marked)))
        n_carriers = (genes[f"n_sites_{modification}"] > 0).sum()
        rows.append(
            key
            | dict(species=species, n_samples=matrix.shape[1], n_carriers=n_carriers)
            | primary
            | {f"robustness_{name}": value for name, value in robust.items()}
        )
    return rows, pd.concat(sets)


# every unit and layer in parallel
jobs = [delayed(test_unit)(unit, layer, table) for layer in LAYERS for unit, table in effects[layer].groupby("unit")]
results = Parallel(n_jobs=WORKERS)(jobs)
signatures = pd.DataFrame([row for rows, _ in results for row in rows])
utils.write(pd.concat([sets for _, sets in results]), OUT / "matched_genes.tsv.gz")
assert np.isfinite(signatures[["z", "robustness_z"]]).all(axis=None)
context_of = {
    (unit, layer): context
    for context, layer, units in utils.CONTEXTS
    for unit in units.split()
}
signatures["cohort"] = signatures.unit.map(cohort_of)
signatures["context"] = [context_of.get(key, "none") for key in zip(signatures.unit, signatures.layer, strict=True)]
signatures["age_scale"] = signatures.species.map(utils.AGE_SCALE)
order = ["unit", "cohort", "species", "context", "layer", "modification", "n_samples", "n_carriers", "n_set"]
unit_stats = ["direction", "z", "p", "diff_pct"]
signatures = signatures[order + unit_stats + ["age_scale"] + [f"robustness_{c}" for c in ["n_set", *unit_stats]]]
utils.write(signatures, OUT / "unit_signatures.tsv")

# contexts, then all units of a species and layer
families = [(context, layer, units.split()) for context, layer, units in utils.CONTEXTS]
for (species, layer), group in signatures.groupby(["species", "layer"]):
    families.append((f"{species}_all", layer, group.unit.unique()))

# cohorts: units averaged first, weighted Stouffer with the square root of the samples
rows = []
for context, layer, units in families:
    family = signatures[(signatures.layer == layer) & signatures.unit.isin(units)]
    assert set(units) == set(family.unit), f"{context} {layer}: units without this layer"
    for modification, group in family.groupby("modification", sort=False):
        row = dict(context=context, layer=layer, modification=modification, n_units=len(group))
        row["n_cohorts"] = group.cohort.nunique()
        for prefix in ["", "robustness_"]:
            cohorts = group.groupby("cohort").agg(
                z=(f"{prefix}z", "mean"), n=("n_samples", "mean"), diff_pct=(f"{prefix}diff_pct", "mean")
            )
            pooled = utils.stouffer(cohorts.z, np.sqrt(cohorts.n))
            row[f"{prefix}n_agree"] = (np.sign(group[f"{prefix}z"]) == np.sign(pooled["z"])).sum()
            row[f"{prefix}direction"] = "Up" if pooled["z"] > 0 else "Down"
            row[f"{prefix}z"] = pooled["z"]
            row[f"{prefix}p"] = pooled["p"]
            row[f"{prefix}diff_pct"] = cohorts.diff_pct.mean()
        row["age_scale"] = group.age_scale.iloc[0]
        rows.append(row)
combined = pd.DataFrame(rows)
# BH over the four modifications of a context and layer, for the primary and for the robustness test separately
for prefix in ["", "robustness_"]:
    combined[f"{prefix}q"] = combined.groupby(["context", "layer"])[f"{prefix}p"].transform(utils.bh)
combined_stats = ["n_agree", "direction", "z", "p", "q", "diff_pct"]
key = ["context", "layer", "modification"]
combined = combined[
    [*key, "n_units", "n_cohorts", *combined_stats, "age_scale", *[f"robustness_{c}" for c in combined_stats]]
]
utils.write(combined, OUT / "combined_signatures.tsv")
print(combined[[*key, "n_units", "n_cohorts", *combined_stats]].round(3).to_string(index=False))
print("robustness: matching also on neuron share")
print(combined[[*key, *[f"robustness_{c}" for c in combined_stats]]].round(3).to_string(index=False))
