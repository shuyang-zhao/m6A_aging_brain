# TabPFN model of the combined age z of each gene from its transcript features, on held-out chromosomes, against
# LightGBM, an elastic net and a shuffled target; m6A contribution as the held-out R2 gain from m6A sites and reader
# binding and as the in silico deletion of the sites, both against 100 copies of the two groups permuted jointly across
# genes; transfer of the mouse protein model to human genes
#   50_transcript_ageing_model.py inputs        gene tables of the targets and the feature tables of all genes (CPU)
#   50_transcript_ageing_model.py TABLE [A-B]   TabPFN on one gene table, or on its permuted copies A to B (one GPU)
#   50_transcript_ageing_model.py               LightGBM and the elastic net, and the result tables (CPU)
# reads scripts 09, 10 and 40, the prepared units and the site atlas; writes output/50_transcript_ageing_model/
import sys

import numpy as np
import pandas as pd
from scipy import stats as sps
from sklearn.metrics import r2_score
from sklearn.model_selection import GroupKFold

import utils

OUT = utils.OUT / "50_transcript_ageing_model"
INPUTS, RUNS = OUT / "inputs", OUT / "runs"
TABPFN = "V3_5"  # TabPFN-3.5 regressor, weights tabpfn-v3.5-20260909 (non-commercial licence), default settings
N_OUTER, N_INNER = 5, 3  # chromosome groups held out, and splits of the training genes for tuning (LightGBM, EN)
N_SEARCH = 50  # LightGBM random-search budget
N_PERMUTED = 100  # permuted copies of the m6A block
N_BOOTSTRAP = 1000
MAX_BIN = 63  # histogram bins per feature, the LightGBM setting for speed
WORKERS = -1  # every available core
SEARCH_SPACE = {
    "n_estimators": sps.randint(100, 500),
    "learning_rate": sps.loguniform(0.01, 0.1),
    "num_leaves": sps.randint(8, 64),
    "min_child_samples": sps.randint(10, 100),
    "subsample": sps.uniform(0.5, 0.5),
    "colsample_bytree": sps.uniform(0.3, 0.7),
    "reg_lambda": sps.loguniform(1e-3, 10),
}
L1_RATIOS = [0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1]
# species, layer and context of the units (None: every unit of the species with the layer)
TARGETS = {
    "mouse_protein": ("mouse", "protein_only", "mouse_protein"),
    "human_mrna": ("human", "mrna_only", None),
    "human_neurons": ("human", "mrna_only", "human_neurons"),
    "mouse_mrna_adjusted": ("mouse", "mrna_adjusted", "mouse_protein"),
    "human_protein": ("human", "protein_only", "human_protein"),
}
regulators = utils.read(utils.ROOT / "ref/gene_sets/rna_mod_regulators.tsv")
READERS = regulators[(regulators.role == "reader") & regulators.modification.str.contains("m6A")].human_symbol
REGIONS = ["5UTR", "CDS", "3UTR"]
SITES = [f"m6A_{region}" for region in ["all", *REGIONS]]
GROUPS = {
    "Transcript structure": [
        "log_tx_len",
        "log_gene_len",
        "log1p_intron_len",
        "log1p_utr5_len",
        "log_cds_len",
        "log1p_utr3_len",
        "log_n_exons",
        "log1p_internal_exon_mean_len",
        "log_first_exon_len",
        "log_last_exon_len",
        "junction_density_per_kb",
        "gc_tx",
        "gc_utr5",
        "gc_cds",
        "gc_utr3",
        "kr_max10",
        "kr_frac",
        "net_charge_per_res",
    ],
    "Expression and turnover": ["expression", "mrna_half_life_pc1", "log2_protein_lifetime"],
    "Cell-type origin": ["ct_frac_neuron", "ct_frac_astro", "ct_frac_oligo", "ct_frac_micro"],
    "Functional program": ["myelin", "lysosome", "mitochondrion", "synapse"],
    "m6A sites": SITES,
    "m6A reader binding": [f"bound_{reader}" for reader in READERS],
    "Sequence representation": ["ribonn_te_mean", "ribonn_te_brain"] + [f"orthrus_pc{i}" for i in range(1, 65)],
}
GROUP_OF = {column: group for group, columns in GROUPS.items() for column in columns}
BLOCKS = {  # feature groups removed to measure their gain; the first block is also tested against permuted copies
    "m6A sites and binding": ["m6A sites", "m6A reader binding"],
    "m6A sites": ["m6A sites"],
    "m6A reader binding": ["m6A reader binding"],
    **{group: [group] for group in GROUPS if not group.startswith("m6A")},  # every other group, one at a time
}
TESTED = "m6A sites and binding"
# features of both species: no reader binding (CLIP of other readers in each species), no sequence representation
SHARED = [c for c in GROUP_OF if GROUP_OF[c] not in ("m6A reader binding", "Sequence representation")]
DELETED = ["all", *REGIONS]


def tabpfn(X, y):
    from tabpfn import TabPFNRegressor  # GPU steps only
    from tabpfn.constants import ModelVersion

    model = TabPFNRegressor.create_default_for_version(ModelVersion[TABPFN])
    # no check of the table size against the pretraining limits (about 15,000 training genes here)
    model.set_params(device="cuda", random_state=utils.SEED, ignore_pretraining_limits=True)
    model.fit(X, y)
    return model


def heldout(X, y, fold, variants=()):
    # each chromosome group predicted by TabPFN in context of the other groups (no tuning), also from variants of the
    # held-out features: test rows attend only to the training rows, so one call serves all
    pred = np.empty((1 + len(variants), len(y)))
    for k in range(N_OUTER):
        model = tabpfn(X[fold != k], y[fold != k])
        pred[:, fold == k] = model.predict(pd.concat([v[fold == k] for v in (X, *variants)])).reshape(len(pred), -1)
    return pred


def deletions(X):
    # the m6A sites of every gene deleted in silico, all of them or those of one region; reader binding is kept
    out = [X.assign(**dict.fromkeys(SITES, 0.0))]
    for region in REGIONS:
        rest = np.log1p(np.maximum(np.expm1(X.m6A_all) - np.expm1(X[f"m6A_{region}"]), 0))
        out.append(X.assign(m6A_all=rest, **{f"m6A_{region}": 0.0}))
    return out


def summary(y, pred):
    # Spearman rho and R2 of the held-out predictions, 95% intervals from resampling genes
    rng = np.random.default_rng(utils.SEED)
    draws = rng.integers(0, len(y), (N_BOOTSTRAP, len(y)))
    rho = [sps.spearmanr(y[i], pred[i]).statistic for i in draws]
    r2 = [r2_score(y[i], pred[i]) for i in draws]
    test = sps.spearmanr(y, pred)
    (rho_low, r2_low), (rho_high, r2_high) = np.percentile([rho, r2], [2.5, 97.5], axis=1)
    return dict(n=len(y), rho=test.statistic, rho_ci_low=rho_low, rho_ci_high=rho_high, p=test.pvalue) | dict(
        r2=r2_score(y, pred), r2_ci_low=r2_low, r2_ci_high=r2_high
    )


def features_of(table):
    return table[[c for c in GROUP_OF if c in table]]


def permuted_runs(name):
    # R2 and deletion effects of the permuted copies of a table
    return pd.concat(utils.read(f) for f in sorted(RUNS.glob(f"{name}_permuted_*.tsv"))).sort_values("draw")


if sys.argv[1:] == ["inputs"]:
    # m6A sites on the canonical transcript by region, the stop-codon window of the atlas split at the CDS end: in mouse
    # the GLORI sites of the two maps counted once across maps (as in script 09), in human the brain m6A-seq peaks at
    # their midpoints (the human load of script 21)
    position = ["chrom", "start", "end", "strand"]
    columns = ["map", "modification", "species", *position, "gene_id", "region", "tx_offset", "ref_ok"]
    atlas = utils.read(utils.RESOURCES / "site_atlas" / "sites.parquet", columns=columns)
    atlas = atlas[atlas.gene_id.notna() & atlas.tx_offset.notna() & atlas.ref_ok & (atlas.modification == "m6A")]
    sites = atlas.drop_duplicates(["species", *position]).astype({"tx_offset": int})
    transcripts = [utils.read(utils.RESOURCES / f"gencode_{s}" / "transcripts.parquet") for s in ("mouse", "human")]
    cds_end = pd.concat(transcripts).set_index("gene_id").cds_end
    stop = sites.region == "stop"
    sites.loc[stop, "region"] = np.where(sites.tx_offset[stop] < sites.gene_id[stop].map(cds_end), "CDS", "3UTR")
    sites = sites[sites.region.isin(REGIONS)]
    maps = {"mouse": sites[sites.species == "Mus musculus"], "human": sites[sites.species == "Homo sapiens"]}

    # reader binding: 1 for the transcripts with an m6A site within 50 nt of the reader's CLIP peaks (the bound
    # transcripts of script 40), one column per reader with CLIP in the species
    bound = utils.read(utils.OUT / "40_target_sets" / "target_genes.tsv.gz")
    bound = bound[(bound.window == utils.PEAK_WINDOW) & (bound.group == "target") & (bound.modification == "m6A")]
    binding = {
        species: pd.crosstab(rows.gene_id, "bound_" + rows.regulator) for species, rows in bound.groupby("species")
    }

    # transcript features of script 09 with log(1 + m6A sites) and m6A reader binding
    features = {}
    for species, rows in maps.items():
        table = utils.read(utils.OUT / "09_transcript_features" / f"{species}_features.tsv.gz").set_index("gene_id")
        counts = rows.groupby(["gene_id", "region"]).size().unstack(fill_value=0).reindex(table.index, columns=REGIONS)
        counts = counts.fillna(0).set_axis(SITES[1:], axis=1)
        counts.insert(0, "m6A_all", counts.sum(axis=1))
        features[species] = table.join(np.log1p(counts)).join(binding[species].reindex(table.index, fill_value=0))

    effects = {
        layer: utils.read(
            utils.OUT / "10_protein_age_effects" / f"{layer}.tsv.gz",
            usecols=["unit", "gene_id", "n", "t", "p", "status"],
        ).query("status == 'ok'")
        for layer in ("protein_only", "mrna_only", "mrna_adjusted")
    }
    cohort_of = utils.read(utils.PREPARED / "units.tsv", index_col="unit").dataset

    def target_table(species, layer, context):
        # z of each gene from its limma p and the sign of t; units of a cohort averaged, cohorts combined by Stouffer
        # weighted by the square root of the samples; expression: normal scores of the mean level, averaged over units
        table = effects[layer]
        if context is None:
            units = table.unit[table.gene_id.str.startswith("ENSMUSG" if species == "mouse" else "ENSG")].unique()
        else:
            units = next(listed.split() for name, lay, listed in utils.CONTEXTS if (name, lay) == (context, layer))
        table = table[table.unit.isin(units)]
        table = table.assign(z=np.sign(table.t) * sps.norm.isf(table.p / 2), cohort=table.unit.map(cohort_of))
        cohorts = table.groupby(["gene_id", "cohort"])[["z", "n"]].mean()
        z = cohorts.groupby("gene_id").apply(lambda d: utils.stouffer(d.z, np.sqrt(d.n))["z"])
        scores = []
        for unit in units:
            protein, rna, _, _ = utils.load_unit(unit)
            level = (protein if layer == "protein_only" else rna).mean(axis=1)
            scores.append(pd.Series(sps.norm.ppf((level.rank() - 0.5) / level.count()), index=level.index))
        expression = pd.concat(scores, axis=1).mean(axis=1).rename("expression")
        print(f"{species} {layer} {context}: {len(units)} units")
        return pd.concat([z.rename("target"), expression], axis=1, join="inner")

    # one table per target: chromosome folds, the target shuffled within each fold and the gene order of each permuted
    # copy of the m6A block
    INPUTS.mkdir(parents=True, exist_ok=True)
    for species, table in features.items():
        table.to_parquet(INPUTS / f"features_{species}.parquet")  # every gene, with or without a target (script 52)
    rng = np.random.default_rng(utils.SEED)
    for name, (species, layer, context) in TARGETS.items():
        targets = target_table(species, layer, context)
        table = pd.concat([features[species], targets], axis=1, join="inner").rename_axis("gene_id")
        fold = np.empty(len(table), int)
        for k, (_, test) in enumerate(GroupKFold(N_OUTER).split(table, groups=table.chrom)):
            fold[test] = k
        shuffled = table.target.to_numpy().copy()
        for k in range(N_OUTER):
            shuffled[fold == k] = rng.permutation(shuffled[fold == k])
        orders = np.array([rng.permutation(len(table)) for _ in range(N_PERMUTED)], dtype=np.int32)
        np.save(INPUTS / f"permutations_{name}.npy", orders)
        columns = ["symbol", "chrom", "fold", "target", "shuffled", *features_of(table)]
        frame = table.assign(fold=fold, shuffled=shuffled)[columns]
        frame.to_parquet(INPUTS / f"{name}.parquet")
        print(name, frame.shape, "genes with sites", (frame.m6A_all > 0).sum())

    # the human genes with a one-to-one mouse ortholog in the mouse protein table, with the fold of that ortholog, for
    # the transfer
    mouse = pd.read_parquet(INPUTS / "mouse_protein.parquet")
    ids = {}
    for species in ("mouse", "human"):
        keys = features[species].symbol.str.upper()
        ids[species] = pd.Series(keys.index, index=keys)[~keys.duplicated(keep=False).to_numpy()]  # unambiguous
    pairs = utils.orthologs()
    pairs = pairs.assign(mouse_id=pairs.mouse.map(ids["mouse"]), human_id=pairs.human.map(ids["human"]))
    pairs = pairs[pairs.mouse_id.isin(mouse.index)].dropna().set_index("human_id")
    pairs["mouse_fold"] = pairs.mouse_id.map(mouse.fold)
    utils.write(pairs[["mouse_id", "mouse_fold"]], INPUTS / "transfer_pairs.tsv")

elif sys.argv[1:]:
    # TabPFN on one gene table (GPU): held-out predictions with the sites deleted, on the shuffled target and without
    # each block; with draws A-B, R2 and mean deletion effects of the permuted copies
    name, draws = sys.argv[1], sys.argv[2:]
    table = pd.read_parquet(INPUTS / f"{name}.parquet")
    X, y, fold = features_of(table), table.target.to_numpy(), table.fold.to_numpy()
    RUNS.mkdir(parents=True, exist_ok=True)
    if draws:
        first, last = map(int, draws[0].split("-"))
        orders = np.load(INPUTS / f"permutations_{name}.npy")
        tested = [c for c in X if GROUP_OF[c] in BLOCKS[TESTED]]
        rows = []
        for draw in range(first, last + 1):
            permuted = X.copy()
            permuted[tested] = X[tested].to_numpy()[orders[draw - 1]]
            pred = heldout(permuted, y, fold, deletions(permuted))
            marked = permuted.m6A_all.to_numpy() > 0  # genes carrying moved sites
            effect = {f"delete_{r}": (pred[0] - p)[marked].mean() for r, p in zip(DELETED, pred[1:], strict=True)}
            rows.append(dict(target=name, draw=draw, r2=r2_score(y, pred[0])) | effect)
            print(rows[-1], flush=True)
        utils.write(pd.DataFrame(rows), RUNS / f"{name}_permuted_{first}-{last}.tsv")
        sys.exit()
    pred = heldout(X, y, fold, deletions(X))
    runs = table[["symbol", "chrom", "fold", "target", "shuffled"]].assign(tabpfn=pred[0])
    runs = runs.assign(**{f"deleted_{r}": p for r, p in zip(DELETED, pred[1:], strict=True)})
    runs["shuffled_target"] = heldout(X, table.shuffled.to_numpy(), fold)[0]
    for block, members in BLOCKS.items():
        runs[f"without {block}"] = heldout(X.drop(columns=[c for c in X if GROUP_OF[c] in members]), y, fold)[0]
    utils.write(runs, RUNS / f"{name}.parquet")
    print(name, "R2", r2_score(y, pred[0]))
    if name == "mouse_protein":
        # transfer: human genes predicted by the shared-feature mouse model fitted without the fold of their ortholog
        pairs = utils.read(INPUTS / "transfer_pairs.tsv", index_col="human_id")
        transfer = []
        for human in ("human_protein", "human_mrna"):
            target = pd.read_parquet(INPUTS / f"{human}.parquet")
            target = target[target.index.isin(pairs.index)]
            folds = pairs.mouse_fold[target.index].to_numpy()
            predicted = {}
            for model, train_y in (("shared", y), ("shared_shuffled", table.shuffled.to_numpy())):
                predicted[model] = np.empty(len(target))
                for k in range(N_OUTER):
                    mouse = tabpfn(X.loc[fold != k, SHARED], train_y[fold != k])
                    predicted[model][folds == k] = mouse.predict(target.loc[folds == k, SHARED])
            transfer.append(target[["symbol"]].assign(target=human, observed=target.target, **predicted))
        utils.write(pd.concat(transfer), RUNS / "transfer.parquet")

else:
    # LightGBM (tuned by random search on inner chromosome folds) and the elastic net on the same inputs, then the
    # tables: performance, gain from each block, deletion effects with their permutation nulls, transfer
    from joblib import Parallel, delayed
    from lightgbm import LGBMRegressor
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import ElasticNetCV
    from sklearn.model_selection import RandomizedSearchCV
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    def comparison(X, y, groups, linear):
        # each chromosome group predicted by a model tuned and fitted on the other groups
        pred, params = np.empty(len(y)), []
        for train, test in GroupKFold(N_OUTER).split(X, y, groups):
            inner = list(GroupKFold(N_INNER).split(X.iloc[train], groups=groups[train]))
            if linear:
                net = ElasticNetCV(l1_ratio=L1_RATIOS, cv=inner, max_iter=10000)
                model = make_pipeline(SimpleImputer(strategy="median"), StandardScaler(), net)
            else:
                trees = LGBMRegressor(max_bin=MAX_BIN, subsample_freq=1, n_jobs=1, random_state=utils.SEED, verbose=-1)
                model = RandomizedSearchCV(trees, SEARCH_SPACE, n_iter=N_SEARCH, cv=inner, random_state=utils.SEED)
            model.fit(X.iloc[train], y[train])
            pred[test] = model.predict(X.iloc[test])
            params.append({} if linear else dict(inner_r2=model.best_score_) | model.best_params_)
        return pred, params

    tables = {name: pd.read_parquet(INPUTS / f"{name}.parquet") for name in TARGETS}
    jobs = {
        (name, model): delayed(comparison)(
            features_of(t), t.target.to_numpy(), t.chrom.to_numpy(), model == "elastic_net"
        )
        for name, t in tables.items()
        for model in ("lightgbm", "elastic_net")
    }
    compared = dict(zip(jobs, Parallel(n_jobs=WORKERS, verbose=5)(jobs.values()), strict=True))

    rows, ablation, nulls, params, predictions = [], [], [], [], []
    for name, table in tables.items():
        runs = utils.read(RUNS / f"{name}.parquet")
        y, n_features = runs.target.to_numpy(), features_of(table).shape[1]
        rows.append(dict(target=name, model="tabpfn", n_features=n_features) | summary(y, runs.tabpfn.to_numpy()))
        for model in ("lightgbm", "elastic_net"):
            rows.append(dict(target=name, model=model, n_features=n_features) | summary(y, compared[name, model][0]))
        shuffled = runs.shuffled.to_numpy()
        rows.append(
            dict(target=name, model="shuffled_target", n_features=n_features)
            | summary(shuffled, runs.shuffled_target.to_numpy())
        )
        params += [dict(target=name, fold=k + 1) | p for k, p in enumerate(compared[name, "lightgbm"][1])]
        r2 = r2_score(y, runs.tabpfn)
        permuted = permuted_runs(name)
        assert permuted.draw.tolist() == list(range(1, N_PERMUTED + 1))
        for block in BLOCKS:
            without = r2_score(y, runs[f"without {block}"])
            row = dict(target=name, block=block, r2=r2, r2_without=without, delta_r2=r2 - without)
            if block == TESTED:
                null = permuted.r2.to_numpy() - without
                nulls.append(pd.DataFrame(dict(target=name, draw=permuted.draw, delta_r2=null)))
                low, high = np.percentile(null, [2.5, 97.5])
                p = utils.perm_pvalue(r2 - without, null, "greater")
                row |= dict(null_mean=null.mean(), null_low=low, null_high=high, p=p)
            ablation.append(row)
        predictions.append(
            runs[["symbol", "chrom", "fold", "target"]]
            .rename(columns={"target": "observed"})
            .assign(
                target=name,
                fold=runs.fold + 1,
                tabpfn=runs.tabpfn,
                lightgbm=compared[name, "lightgbm"][0],
                elastic_net=compared[name, "elastic_net"][0],
            )
        )
    performance = pd.DataFrame(rows)
    utils.write(performance, OUT / "performance.tsv")
    utils.write(pd.DataFrame(ablation), OUT / "block_ablation.tsv")
    utils.write(pd.concat(nulls), OUT / "block_ablation_permuted.tsv")
    utils.write(pd.DataFrame(params), OUT / "hyperparameters.tsv")
    utils.write(pd.concat(predictions), OUT / "predictions.tsv.gz")
    print(performance.round(3).to_string(index=False))
    print(pd.DataFrame(ablation).round(4).to_string(index=False))

    # in silico deletion: mean of prediction minus prediction without the sites over the genes with sites, against
    # the same mean in the permuted copies (two-sided)
    deletion, permuted_effects, genes = [], [], []
    for name in TARGETS:
        table, runs = pd.read_parquet(INPUTS / f"{name}.parquet"), utils.read(RUNS / f"{name}.parquet")
        marked = table.m6A_all.to_numpy() > 0
        effects = pd.DataFrame({r: runs.tabpfn - runs[f"deleted_{r}"] for r in DELETED})
        permuted = permuted_runs(name)
        for region in DELETED:
            row = dict(target=name, region=region, n_genes=marked.sum(), effect=effects[region][marked].mean())
            null = permuted[f"delete_{region}"].to_numpy()
            low, high = np.percentile(null, [2.5, 97.5])
            row |= dict(null_mean=null.mean(), null_low=low, null_high=high, p=utils.perm_pvalue(row["effect"], null))
            deletion.append(row)
        long = permuted.melt(["target", "draw"], [f"delete_{r}" for r in DELETED], "region", "effect")
        permuted_effects.append(long.assign(region=long.region.str.removeprefix("delete_")))
        gene_rows = table.loc[marked, ["symbol"]].assign(target=name, observed=runs.target, tabpfn=runs.tabpfn)
        genes.append(gene_rows.join(effects[marked].add_prefix("effect_")))
    utils.write(pd.DataFrame(deletion), OUT / "deletion.tsv")
    utils.write(pd.concat(permuted_effects), OUT / "deletion_permuted.tsv")
    utils.write(pd.concat(genes), OUT / "deletion_genes.tsv.gz")
    print(pd.DataFrame(deletion).round(4).to_string(index=False))

    # transfer: rank correlation only, as z scales differ between species
    transferred = utils.read(RUNS / "transfer.parquet")
    transfer = []
    for name, rows in transferred.groupby("target", sort=False):
        for model in ("shared", "shared_shuffled"):
            stats = summary(rows.observed.to_numpy(), rows[model].to_numpy())
            transfer.append(dict(target=name, model=model) | {k: v for k, v in stats.items() if not k.startswith("r2")})
    utils.write(pd.DataFrame(transfer), OUT / "transfer.tsv")
    utils.write(transferred, OUT / "transfer_predictions.tsv.gz")
    print(pd.DataFrame(transfer).round(3).to_string(index=False))
