# Concordance of the age effects of the human units (every layer they have: protein, and for BrainGVEX also
# mRNA-adjusted protein and mRNA) with those of mouse protein, mRNA-adjusted protein and mRNA on one-to-one orthologs,
# with permutation P values from the age-permuted mouse fits
# reads script 10 age effects and script 12 permuted t; writes output/13_species_concordance/
import numpy as np
import pandas as pd
from scipy import stats as sps

import utils

AGE_EFFECTS = utils.OUT / "10_protein_age_effects"
PERMUTED = utils.OUT / "12_permutation_fdr"
OUT = utils.OUT / "13_species_concordance"
LAYERS = ["mrna_adjusted", "protein_only", "mrna_only"]
MOUSE = ["chdi2016_expanded", "takasugi2024"]
HUMAN = [*utils.HUMAN_PROTEIN_UNITS, "luo2024_bankadj"]


def age_effects(unit, layer):
    # genes a unit estimated in one layer
    return utils.by_symbol(tables[layer].query("unit == @unit and status == 'ok'"))


def to_mouse(frame):
    # human genes re-indexed by mouse ortholog
    frame = frame[frame.index.isin(human_to_mouse.index)]
    return frame.set_axis(frame.index.map(human_to_mouse))


def compare(x, y, q=None, null=None):
    # rank correlation on shared genes, permutation p from null draws
    genes = x.index.intersection(y.index)
    x, y = x[genes], y[genes]
    rho = sps.spearmanr(x, y)[0]
    half_width = utils.Z95 * np.sqrt(1.06 / (len(genes) - 3))  # 95% interval on the Fisher z scale, Spearman variance
    same_sign = np.sign(x) == np.sign(y)
    row = dict(
        n=len(genes),
        rho=rho,
        ci_low=np.tanh(np.arctanh(rho) - half_width),
        ci_high=np.tanh(np.arctanh(rho) + half_width),
        sign_agreement=same_sign.mean(),
    )
    if q is not None:
        hits = (q[genes] < utils.Q).to_numpy()
        row |= dict(n_hits=hits.sum(), hits_same_sign=same_sign[hits].sum())
    if null is not None:
        draws, rhos = null.reindex(genes), []
        for draw in draws:
            ok = draws[draw].notna().to_numpy()
            rhos.append(sps.spearmanr(x[ok], draws.loc[ok, draw])[0])
        row |= dict(null_mean=np.mean(rhos), null_sd=np.std(rhos, ddof=1), p=utils.perm_pvalue(rho, rhos))
        row["null_rho"] = pd.Series(rhos, index=draws.columns)
    return row


tables = {layer: utils.read(AGE_EFFECTS / f"{layer}.tsv.gz") for layer in LAYERS}
orthologs = utils.orthologs()
human_to_mouse = orthologs.set_index("human").mouse
assert orthologs.human.is_unique and orthologs.mouse.is_unique
units = {layer: set(table.unit) for layer, table in tables.items()}
assert set(MOUSE) <= units["mrna_adjusted"] and set(HUMAN) <= units["protein_only"]

# mouse layers on common genes, permuted t of each layer
mouse, null = {}, {}
for unit in MOUSE:
    mouse[unit] = pd.concat(
        {layer: age_effects(unit, layer)[["estimate", "t", "q"]] for layer in LAYERS}, axis=1, join="inner"
    ).dropna()
    for layer in LAYERS:
        estimated = age_effects(unit, layer)
        permuted = utils.read(PERMUTED / f"null_t_{unit}_{layer}.parquet")
        null[unit, layer] = permuted.loc[estimated.gene_id].set_axis(estimated.index)
        assert null[unit, layer].shape[1] == utils.N_PERMUTATIONS
    assert mouse[unit].index.is_unique

# human layers per unit
human = {}
for unit in HUMAN:
    for layer in [layer for layer in LAYERS if unit in units[layer]]:
        human[unit, layer] = to_mouse(age_effects(unit, layer)[["estimate", "t"]])

# human-mouse pairs: key columns, human vector and mouse unit
pairs = [
    (dict(human_unit=unit, human_layer=layer, mouse_unit=mouse_unit), h, m)
    for (unit, layer), h in human.items()
    for mouse_unit, m in mouse.items()
]

# human t against each mouse layer
rows, draws = [], []
for key, h, m in pairs:
    for layer in LAYERS:
        row = compare(h.t, m[layer].t, m[layer].q, null[key["mouse_unit"], layer])
        key_layer = key | dict(mouse_layer=layer)
        draws.append(row.pop("null_rho").rename("rho").rename_axis("draw").reset_index().assign(**key_layer))
        rows.append(key_layer | row)
species = pd.DataFrame(rows)
utils.write(species, OUT / "species_concordance.tsv")
utils.write(pd.concat(draws), OUT / "null_correlations.tsv.gz")
shown = ["human_unit", "human_layer", "mouse_unit", "mouse_layer", "n", "rho", "sign_agreement", "p"]
print(species[shown].round(3).to_string(index=False))
