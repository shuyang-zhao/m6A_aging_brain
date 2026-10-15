# Age effect of every RNA-modification regulator in each context: cohort means, random-effects meta-analysis, BH
# reads output/10_protein_age_effects and ref/gene_sets/rna_mod_regulators.tsv; writes output/30_regulator_screen/
import pandas as pd

import utils

OUT = utils.OUT / "30_regulator_screen"

# regulators without the ADAR family (A-to-I) and the genes listed as not acting on RNA
regulators = utils.read(utils.ROOT / "ref/gene_sets/rna_mod_regulators.tsv")
regulators = regulators[(regulators.tier != "excluded") & (regulators.family != "ADAR")]
assert regulators.human_symbol.is_unique
cohort_of = utils.read(utils.PREPARED / "units.tsv", index_col="unit").dataset
effects = {
    layer: utils.read(
        utils.OUT / "10_protein_age_effects" / f"{layer}.tsv.gz",
        usecols=["unit", "gene_id", "estimate", "se", "status"],
    )
    for layer in ("mrna_only", "protein_only", "mrna_adjusted")
}


def regulator_genes(species):
    # regulator by gene id; mouse Mettl2 is shared by METTL2A and METTL2B and goes to the first
    symbols = utils.ensg_to_symbol(species)
    ids = pd.DataFrame({"gene_id": symbols.index, "key": symbols.str.upper().to_numpy()})
    keys = pd.DataFrame({"regulator": regulators.human_symbol, "key": regulators[f"{species}_symbol"].str.upper()})
    return keys.merge(ids, on="key").drop_duplicates("gene_id").set_index("gene_id").regulator


genes = {species: regulator_genes(species) for species in utils.AGE_SCALE}
names = pd.concat(genes.values())
print({species: len(found) for species, found in genes.items()}, "regulators with a gene id")

# unit effects of the regulators, one cohort mean per context, layer, gene and dataset
parts = []
for context, layer, listed in utils.CONTEXTS:
    table, listed = effects[layer], listed.split()
    assert set(listed) <= set(table.unit), f"{context} {layer}: units without this layer"
    rows = table[
        table.unit.isin(listed) & table.gene_id.isin(genes[context.split("_")[0]].index) & (table.status == "ok")
    ]
    assert set(listed) <= set(rows.unit), f"{context} {layer}: units without any regulator"
    parts.append(rows.assign(context=context, layer=layer, cohort=rows.unit.map(cohort_of)))
units = pd.concat(parts)
keys = ["context", "layer", "gene_id", "cohort"]
cohorts = units.groupby(keys).agg(
    n_units=("unit", "size"), units=("unit", ",".join), estimate=("estimate", "mean"), se=("se", "mean")
)
cohorts = cohorts.reset_index()
cohorts.insert(2, "regulator", cohorts.gene_id.map(names))
utils.write(cohorts, OUT / "regulator_cohorts.tsv")

# independent cohorts: DerSimonian-Laird random effects, then BH over the regulators of each context and layer
rows = []
for key, group in cohorts.groupby(["context", "layer", "gene_id"], sort=False):
    pooled = utils.random_effects(group.estimate, group.se)
    rows.append(
        dict(zip(["context", "layer", "gene_id"], key, strict=True)) | pooled | dict(cohorts=",".join(group.cohort))
    )
screen = pd.DataFrame(rows).rename(columns={"est": "estimate", "k": "n_cohorts"})
screen.insert(2, "regulator", screen.gene_id.map(names))
screen["age_scale"] = screen.context.str.split("_").str[0].map(utils.AGE_SCALE)
screen["q"] = screen.groupby(["context", "layer"]).p.transform(utils.bh)
for column, source in [("effect_pct", "estimate"), ("effect_pct_ci_low", "ci_low"), ("effect_pct_ci_high", "ci_high")]:
    screen[column] = utils.pct(screen[source])
columns = ["context", "layer", "regulator", "gene_id", "n_cohorts", "estimate", "se", "ci_low", "ci_high", "effect_pct"]
columns += ["effect_pct_ci_low", "effect_pct_ci_high", "p", "q", "tau2", "i2", "age_scale", "cohorts"]
screen = screen[columns]
utils.write(screen, OUT / "regulator_effects.tsv")

# regulators per family and the hits
for (context, layer), family in screen.groupby(["context", "layer"], sort=False):
    hits = family[family.q < utils.Q].sort_values("p")
    print(f"{context} {layer}: {len(family)} regulators, {len(hits)} with q < {utils.Q}")
    print(
        hits[["regulator", "n_cohorts", "effect_pct", "effect_pct_ci_low", "effect_pct_ci_high", "q"]]
        .round(3)
        .to_string(index=False)
    )
