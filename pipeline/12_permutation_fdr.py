# Age-label permutations within blocks: BH discoveries per draw and the permuted estimate and t of every layer
# reads the prepared units; writes output/12_permutation_fdr/ discoveries.tsv and
# null_<statistic>_<unit>_<layer>.parquet (estimate or t, genes x draws); units with all files are skipped
import zlib

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import model
import utils

AGE_GROUPS = ("middle_", "early_", "late_")  # columns that code age groups and move with age
SALT = "m2a"  # arbitrary label hashed into the seed of every draw; kept to reproduce the stored permutations
NULL_STATS = ["estimate", "t"]  # saved for every draw
WORKERS = 40
OUT = utils.OUT / "12_permutation_fdr"


def rng_for(*keys):
    # seed of a draw from the global seed and a key
    return np.random.default_rng([utils.SEED, zlib.crc32("|".join(map(str, keys)).encode())])


def permute_within(values, blocks, rng):
    # shuffle values within each block
    values = np.asarray(values)
    out = values.copy()
    blocks = np.asarray(blocks).astype(str)
    for block in np.unique(blocks):
        rows = np.where(blocks == block)[0]
        out[rows] = values[rng.permutation(rows)]
    return out


def age_columns(design):
    # age and the columns that code its groups
    return ["age"] + [c for c in design["numeric"] if c.startswith(AGE_GROUPS)]


def donor_pool(members):
    # donors of the units of a dataset: block and age coding, first unit first
    columns = age_columns(data[members[0]][3])
    assert all(age_columns(data[unit][3]) == columns for unit in members)
    pool = pd.concat([data[unit][2][["permutation_block", *columns]] for unit in members])
    # a donor in several units has the same age coding in each, up to rounding
    assert np.allclose(pool[columns] - pool.groupby(level=0)[columns].transform("first"), 0, atol=1e-3)
    return pool[~pool.index.duplicated()]


def donor_groups(members):
    # units that share donors form a group, first unit first; a unit without shared donors stands alone
    groups = []
    for unit in members:
        donors = set(data[unit][2].index)
        linked = [group for group in groups if donors & group[1]]
        merged = (sum([group[0] for group in linked], []) + [unit], donors.union(*[group[1] for group in linked]))
        groups = [group for group in groups if group not in linked] + [merged]
    return [group[0] for group in groups]


def permuted_age(unit, draw):
    # age coding each sample of the unit receives; units that share donors share the donor permutation
    pool = pools[source[unit]]
    seed = rng_for(SALT, source[unit], f"L{draw:03d}")
    order = permute_within(np.arange(len(pool)), pool.permutation_block, seed)
    return pool.drop(columns="permutation_block").to_numpy()[order[pool.index.get_indexer(data[unit][2].index)]]


def permuted_fit(unit, draw, protein, rna, samples, design, age):
    # shuffled age coding, then the fits of script 10
    shuffled = samples.copy()
    shuffled[age_columns(design)] = age
    k = design["hidden_factors"]
    pcs = model.protected_pca(rna, model.design(shuffled, design), k) if k else None
    fits = model.fit_layers(protein, rna, shuffled, design, pcs)
    # discoveries of the first layer
    table = next(iter(fits.values()))
    ok = table.status == "ok"
    stats = dict(unit=unit, draw=draw, age_correlation=np.corrcoef(shuffled.age, samples.age)[0, 1])
    stats.update(n_tested=ok.sum(), n_discoveries=(table.q < utils.Q).sum())
    stats["share_p_below_0_05"] = (table.p[ok] < 0.05).mean()
    return stats, {(layer, name): fit[name].to_numpy(np.float32) for layer, fit in fits.items() for name in NULL_STATS}


def layers_of(unit):
    protein, rna = data[unit][:2]
    return (
        (["mrna_adjusted"] if protein is not None and rna is not None else [])
        + (["protein_only"] if protein is not None else [])
        + (["mrna_only"] if rna is not None else [])
    )


def null_path(statistic, unit, layer):
    return OUT / f"null_{statistic}_{unit}_{layer}.parquet"


# prepared units
units = utils.read(utils.PREPARED / "units.tsv", index_col="unit")
data = {unit: utils.load_unit(unit) for unit in units.index}
for _, _, samples, _ in data.values():
    assert (samples.donor_id.astype(str) == samples.index).all()

# units of a dataset that share donors share the age permutation of the first of them, drawn over all their donors
groups = [group for _, members in units.groupby("dataset", sort=False) for group in donor_groups(list(members.index))]
source = pd.Series({unit: group[0] for group in groups for unit in group})
pools = {group[0]: donor_pool(group) for group in groups}

# units without all their files; a finished unit of the group must already hold all donors of the pool
pending = [
    unit
    for unit in units.index
    if not all(null_path(name, unit, layer).exists() for name in NULL_STATS for layer in layers_of(unit))
]
for unit in pending:
    finished = [u for u in next(group for group in groups if unit in group) if u not in pending]
    covered = set().union(*[data[u][2].index for u in finished])
    assert not finished or pools[source[unit]].index.isin(covered).all(), f"{unit} adds donors to a finished group"
print(f"{len(pending)} of {len(units)} units to permute, {WORKERS} workers")

columns = [f"L{draw:04d}" for draw in range(1, utils.N_PERMUTATIONS + 1)]
for unit in pending:
    protein, rna, samples, design = data[unit]
    if unit in utils.BATCH_CORRECTED:
        protein = model.batch_blup(protein, rna, samples, design)
    draws = range(1, utils.N_PERMUTATIONS + 1)
    jobs = (delayed(permuted_fit)(unit, d, protein, rna, samples, design, permuted_age(unit, d)) for d in draws)
    results = Parallel(n_jobs=WORKERS, batch_size=5)(jobs)
    # permuted estimate and t, genes x draws, one file per layer
    genes = (protein if protein is not None else rna).index
    for layer, statistic in results[0][1]:
        values = np.column_stack([fits[layer, statistic] for _, fits in results])
        utils.write(
            pd.DataFrame(values, index=genes, columns=columns).rename_axis("gene_id"), null_path(statistic, unit, layer)
        )
    # discoveries per draw, the unit added to the table of the earlier units
    stats = pd.DataFrame([row for row, _ in results])
    earlier = utils.read(OUT / "discoveries.tsv") if (OUT / "discoveries.tsv").exists() else stats.iloc[:0]
    utils.write(pd.concat([earlier[earlier.unit != unit], stats], ignore_index=True), OUT / "discoveries.tsv")
    print(f"{unit}: {stats.n_discoveries.mean():.2f} discoveries per draw", flush=True)

discoveries = utils.read(OUT / "discoveries.tsv")
discoveries = discoveries[discoveries.unit.isin(units.index)]  # units of utils.EXCLUDED may remain from earlier runs
utils.write(discoveries, OUT / "discoveries.tsv")
print(discoveries.groupby("unit", sort=False)[["n_discoveries", "share_p_below_0_05"]].agg(["mean", "max"]).round(3))
