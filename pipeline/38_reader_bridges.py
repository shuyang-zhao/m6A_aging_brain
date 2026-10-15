# Bridge tests of the IGF2BP2 and YTHDF readers: does a reader's change with age show in its m6A-bound targets? The
# score of a target set in a sample is the mean of the targets minus the mean of the comparison genes (unbound m6A
# genes, or bound genes without m6A) matched on expression level, genomic length (and m6A load); the control RBPs of
# script 33 with their ENCODE IDR targets give the null (those with enough targets in the test or quantified in the
# unit; n_controls in each row).
# IGF2BP2 (igf2bp2_*; targets: genes bound with m6A in at least 2 of the 4 neural sources; T2, T3 as in script 37):
#   striatum x age in GTEx donors (mixed model with a donor intercept); MSN and OPC against oligodendrocytes and
#   microglia in Green 2025 and Handsaker 2025 (nuclei), and the neuron-versus-glia statistic on the same targets;
#   coupling (script 33 model) and mediation in the MSN and OPC pseudobulks; age slope in striatal proteomes; targets
#   after IGF2BP2 loss.
# YTHDF1-3 (yth_*; target_type bound_set: the preferred binding sets of script 37; m6A_load: the m6A load of every gene
#   against the YTHDF level): age slope of the target scores in the Output, Translation and Stability arms, coupling,
#   mediation of the age slope, targets after YTHDF loss or gain, and the human units where YTHDF3 mRNA falls (mode
#   ythdf3_falls).
# BH families: striatum x age and nuclei over the primary tiers and comparisons; proteomes within species; IGF2BP2 loss
# over all tests; YTHDF target slopes within arm, coupling and mediation within target type and arm, YTHDF loss over
# all tests (the units where YTHDF3 mRNA falls apart from the main units).
# usage: 38_reader_bridges.py [igf2bp2|yth|wildtype]
# reads script 33 (unit data and control RBPs), the compendium of 37, the target sets of 40, the prepared units (08)
# and the pseudobulks of 03 (single_nucleus.py); writes output/38_reader_bridges/, one table per test family
import os

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"
import importlib
import multiprocessing as mp
import sys
from collections import namedtuple

import numpy as np
import pandas as pd
import statsmodels.api as sm
from joblib import Parallel, delayed

import model
import network as nw
import single_nucleus as sn
import utils

m33 = importlib.import_module("33_regulator_network")  # data of a unit, mediation and age-slope helpers

OUT = utils.OUT / "38_reader_bridges"
COMP = utils.OUT / "37_reader_compendium"
WORKERS = 24
NBOOT = 1000
MIN_TARGETS = 50  # genes in the target set and in its comparison set
COMPARISONS = ["unbound_m6A", "bound_unmarked"]
FOCUS = "IGF2BP2"
NEURAL = ["conway2016_H9_eCLIP", "degrauwe2016_MGG8_PARCLIP", "han2022_A172_eCLIP", "preitner2016_E14brain_CLIP"]
GTEX_FOCAL = ["gtex_caudate", "gtex_putamen", "gtex_nucleus_accumbens"]
PROTEOMES = {"mouse": ["keele2025_striatum", "chdi2016_expanded"],
             "human": ["blumenreich2024_striatum", "teng2023_nucleus_accumbens"]}
POPULATIONS = {"green2025_gse277313": ["D1_MSN", "D2_MSN", "MSN_unassigned", "OPC"],
               "handsaker2025_nemo_broadhdcaudate": ["D1_MSN", "D2_MSN", "OPC"]}
GLIA = ["Oligodendrocyte", "Microglia"]
READERS3 = ["YTHDF1", "YTHDF2", "YTHDF3"]
# gene ids of a target set, of the genes bound by the reader, and of the bound genes without m6A
Genes = namedtuple("Genes", "target bound unmarked")
G = {}  # data shared with the worker processes (fork)


# target scores
def weights(genes, target, pool, own, use_load):
    """Weights that turn expression into the target mean minus the comparison mean, matched over 5 level x 4 length
    [x 4 load] bins; returns the weights, both gene counts and the bins (for the bootstrap)."""
    idx = genes.index
    t = idx.isin(list(target)) & ~idx.isin(own)
    c = idx.isin(list(pool)) & ~t & ~idx.isin(own)
    code = np.zeros(len(idx), int)
    for column, q in [("level", 5), ("log_gene_len", 4)] + ([("load", 4)] if use_load else []):
        b = np.full(len(idx), -1)
        b[t | c] = pd.qcut(genes[column][t | c], q, labels=False, duplicates="drop")
        code = code * 10 + b
    n_t, n_c = pd.Series(code[t]).value_counts(), pd.Series(code[c]).value_counts()
    covered = n_t.index.intersection(n_c.index)
    t &= np.isin(code, covered)
    c &= np.isin(code, covered)
    w = np.zeros(len(idx))
    if t.sum() == 0:
        return w, 0, 0, None
    w[t] = 1 / t.sum()
    w[c] = -(n_t.reindex(code[c]).to_numpy() / n_c.reindex(code[c]).to_numpy()) / t.sum()
    return w, int(t.sum()), int(c.sum()), (t, c, code)


def score_weights(genes, sets, comparison, own):
    """genes: level, log_gene_len, load of the unit's genes. Comparison: unbound m6A genes (load > 0, not bound) or
    bound-unmarked genes."""
    if comparison == "unbound_m6A":
        return weights(genes, sets.target, set(genes.index[genes.load > 0]) - set(sets.bound), own, True)
    return weights(genes, sets.target, sets.unmarked, own, False)


def gene_table(species, level):
    f = nw.features[species][["log_gene_len", nw.LOAD[species]]].rename(columns={nw.LOAD[species]: "load"})
    return f.join(level.rename("level"), how="inner").dropna()


def own_of(species, name):
    gene = nw.symbols[species].get(name)
    return [gene] if gene is not None else []


def level_of(unit):
    rna = utils.load_unit(unit)[1]
    return rna[rna.notna().all(axis=1)].mean(axis=1)


def s_u(species, sets):
    """Links of a target set in the form of script 33's gene lists: binding x m6A load of the targets, and the
    bound-unmarked flag."""
    index = nw.features[species].index
    load = nw.features[species][nw.LOAD[species]].to_numpy(float)
    return (pd.Series(np.log1p(np.where(index.isin(list(sets.target)), load, 0)), index),
            pd.Series(index.isin(list(sets.unmarked)).astype(float), index))


def ols_slope(y, X, score_of_own_mrna=None):
    """Coefficient of the last column of X (age) in an OLS fit (HC3); the score of the genes' own mRNA is added as a
    covariate in the paired layer."""
    D = X if score_of_own_mrna is None else np.column_stack([X, score_of_own_mrna])
    fit = sm.OLS(y, D).fit(cov_type="HC3")
    k = X.shape[1] - 1
    return dict(estimate=fit.params[k], se=fit.bse[k], p=fit.pvalues[k])


def unit_design(unit):
    """Protein, mRNA, samples and the design of a unit with the protected hidden factors of script 10."""
    protein, rna, samples, design = utils.load_unit(unit)
    X = model.design(samples, design)
    k = design["hidden_factors"]
    if k:
        pcs = model.protected_pca(rna, X, k)
        X = pd.concat([X.drop(columns="age"), pcs.set_axis([f"PC{i + 1}" for i in range(k)], axis=1), X[["age"]]],
                      axis=1)
    return protein, rna, samples, X


# control RBPs
def control_sets():
    """Targets of the control RBPs: s and u of their ENCODE IDR peaks (human), and the target sets of both species
    (mouse by the bound orthologs)."""
    site = {s: nw.site_map(s) for s in ("human", "mouse")}
    idr = dict(Parallel(n_jobs=WORKERS)(delayed(nw.idr_links)(r, site["human"]) for r in m33.CONTROLS))
    sets = {"human": {}, "mouse": {}}
    for rbp, (s, u, bound) in idr.items():
        sets["human"][rbp] = Genes(set(s.index[s > 0]), bound, set(u.index[(u > 0) & ~(s > 0)]))
        mapped = nw.gene_ids(nw.features["human"].symbol.reindex(list(bound)).dropna(), "human", "mouse")
        s, u, bound = nw.links_from_genes(mapped, "mouse", site["mouse"])
        sets["mouse"][rbp] = Genes(set(s.index[s > 0]), bound, set(u.index[u > 0]))
    return idr, sets


def versus(table, ref_of, value="estimate"):
    """Percentile of each row's value among the control RBPs (ref_of(row): their values), empirical two-sided P."""
    out = []
    for _, row in table.iterrows():
        ref = pd.Series(ref_of(row), dtype=float).dropna()
        ok = len(ref) > 0 and pd.notna(row[value])
        out.append(dict(n_controls=len(ref),
                        pctile_signed=100 * (ref < row[value]).mean() if ok else np.nan,
                        pctile_abs=100 * (ref.abs() < abs(row[value])).mean() if ok else np.nan,
                        control_abs_p95=ref.abs().quantile(0.95) if len(ref) else np.nan,
                        p_control_two_sided=(1 + (ref.abs() >= abs(row[value])).sum()) / (1 + len(ref))
                        if ok else np.nan))
    return pd.DataFrame(out, index=table.index)


def pool(table, keys, label="unit"):
    """Random-effects pooling of estimate and se over the rows of each group of keys."""
    rows = [dict(zip(keys, key if isinstance(key, tuple) else (key,), strict=False))
            | utils.random_effects(g.estimate, g.se) | dict(units=",".join(g[label]))
            for key, g in table.dropna(subset=["estimate", "se"]).query("se > 0").groupby(keys)]
    return pd.DataFrame(rows).rename(columns={"est": "estimate", "k": "n_units"})


def percent(table):
    low, high = ((table.ci_low, table.ci_high) if "ci_low" in table
                 else (table.estimate - utils.Z95 * table.se, table.estimate + utils.Z95 * table.se))
    table = table.assign(estimate_pct=utils.pct(table.estimate), ci_low_pct=utils.pct(low), ci_high_pct=utils.pct(high))
    for column in ("slope_focal", "slope_other"):  # mixed-model slopes of the two groups
        if column in table:
            table[column + "_pct"] = utils.pct(table[column])
    return table


def perturbation_effect(lfc, genes, sets, comparison, own, rng, boot=True):
    """Targets minus matched comparison genes in log2 fold change of a perturbation; bootstrap within the bins for the
    CI."""
    w, nt, nc, detail = score_weights(genes, sets, comparison, own)
    if min(nt, nc) < MIN_TARGETS:
        return None
    est = w @ lfc
    if not boot:
        return dict(log2fc=est)
    t, c, code = detail
    draws = np.zeros(NBOOT)
    for k in np.unique(code[t | c]):
        tv, cv = lfc[t & (code == k)], lfc[c & (code == k)]
        if len(tv) and len(cv):
            draws += len(tv) / nt * (tv[rng.integers(0, len(tv), (NBOOT, len(tv)))].mean(1)
                                     - cv[rng.integers(0, len(cv), (NBOOT, len(cv)))].mean(1))
    se = draws.std(ddof=1)
    return dict(n_targets=nt, n_comparison=nc, log2fc=est, effect_pct=utils.pct(est), ci_low=np.percentile(draws, 2.5),
                ci_high=np.percentile(draws, 97.5), se=se, p=2 * utils.sps.norm.sf(abs(est) / se))


# IGF2BP2
def igf2bp2_sets():
    """Target sets per group: primary (4 neural sources), secondary (all sources) and each source alone; genes bound
    with m6A in at least 2 sources of the group."""
    bound = (utils.read(COMP / "bound_genes.tsv", low_memory=False).query("reader == @FOCUS")
             .dropna(subset=["human_gene_id"]))
    tiers = utils.read(COMP / "tiers_igf2bp2.tsv", low_memory=False)
    every = sorted(bound.source.unique())
    G["t3_by_dataset"] = {d: set(g.human_gene_id.dropna()) for d, g in tiers[tiers.tier == "T3"].groupby("source")}
    sets = {}
    groups = {"primary": NEURAL, "secondary": every} | {s: [s] for s in every}
    for name, sources in groups.items():
        rows = bound[bound.source.isin(sources)]
        marked = rows[rows.T1].groupby("human_gene_id").source.nunique()
        t1 = set(marked[marked >= min(2, len(sources))].index)
        targets = {"T1": t1}
        if name == "primary":
            targets |= {"T2": t1 & set(tiers[tiers.tier == "T2"].human_gene_id.dropna()),
                        "T3": t1 & set().union(*G["t3_by_dataset"].values())}
        ids = set(rows.human_gene_id)
        sets[name] = dict(targets=targets, bound=ids, unmarked=ids - set(rows[rows.T1].human_gene_id))
    return sets


def to_species(ids, species):
    return (set(ids) if species == "human"
            else nw.gene_ids(nw.features["human"].symbol.reindex(list(ids)).dropna(), "human", "mouse"))


def spec(name, tier, species, target=None):
    d = G["sets"][name]
    return Genes(to_species(d["targets"][tier] if target is None else target, species),
                 to_species(d["bound"], species), to_species(d["unmarked"], species))


def tests(species):
    """(source, tier, sets, own genes) of the IGF2BP2 target sets and of the control RBPs (T1) in a species."""
    out = [(f"{FOCUS}:{name}", tier, spec(name, tier, species), own_of(species, FOCUS))
           for name, d in G["sets"].items() for tier in d["targets"]]
    return out + [(f"CTRL:{r}", "T1", g, own_of(species, r)) for r, g in G["ctrl"][species].items()]


def striatum_job(job):
    label, sets, own = job
    parts, nt, nc = [], [], []
    for g in G["gtex"].values():
        w, a, b, _ = score_weights(g["genes"], sets, label["comparison"], own)
        parts.append(w @ g["rna"].loc[g["genes"].index, g["samples"].index].to_numpy())
        nt.append(a)
        nc.append(b)
    row = label | dict(n_targets=int(np.median(nt)), n_comparison=int(np.median(nc)))
    if min(nt) < MIN_TARGETS or min(nc) < MIN_TARGETS:
        return row | dict(status="too few genes")
    return row | model.fit_slopes(np.concatenate(parts), G["X"], G["age"], G["donor"], G["focal"]) | dict(status="ok")


def striatum_age():
    """Age slope of the target scores in the striatal regions against the other regions of GTEx donors (mixed model,
    donor intercept)."""
    gtex = {}
    for unit in [u for u in utils.read(utils.PREPARED / "units.tsv").unit if u.startswith("gtex_")]:
        _, rna, samples, design = utils.load_unit(unit)
        rna = rna[rna.notna().all(axis=1)]
        gtex[unit] = dict(rna=rna, samples=samples, genes=gene_table("human", rna.mean(axis=1)),
                          X=model.design(samples, design).drop(columns="age").add_prefix(f"{unit}:"))
    X = pd.concat([g["X"] for g in gtex.values()]).fillna(0).reset_index(drop=True)
    G.update(gtex=gtex, X=X.loc[:, (X != 0).any()],
             age=np.concatenate([g["samples"].age.to_numpy() for g in gtex.values()]),
             donor=np.concatenate([g["samples"].donor_id.to_numpy() for g in gtex.values()]),
             focal=np.concatenate([[u in GTEX_FOCAL] * len(g["samples"]) for u, g in gtex.items()]))
    jobs = [(dict(source=s, tier=t, comparison=c), sets, own)
            for s, t, sets, own in tests("human") for c in COMPARISONS]
    with mp.get_context("fork").Pool(WORKERS) as pool_:
        table = pd.DataFrame(pool_.map(striatum_job, jobs, chunksize=1))
    # IGF2BP2 mRNA in the units that quantify it
    gene = own_of("human", FOCUS)[0]
    have = [u for u, g in gtex.items() if gene in g["rna"].index]
    mask = np.concatenate([[u in have] * len(gtex[u]["samples"]) for u in gtex])
    y = np.concatenate([gtex[u]["rna"].loc[gene, gtex[u]["samples"].index].to_numpy() for u in have])
    X = G["X"][mask].reset_index(drop=True)
    own = (dict(source=FOCUS + " mRNA", tier="gene", comparison="-", status="ok", n_targets=len(have))
           | model.fit_slopes(y, X.loc[:, (X != 0).any()], G["age"][mask], G["donor"][mask], G["focal"][mask]))
    table = percent(pd.concat([table, pd.DataFrame([own])], ignore_index=True))
    table["is_control"] = table.source.str.startswith("CTRL:")
    ok = ~table.is_control & (table.tier != "gene") & (table.status == "ok")
    primary = ok & (table.source == f"{FOCUS}:primary")
    table.loc[primary, "q"] = utils.bh(table.loc[primary, "p"])
    table.loc[ok, "q_all"] = utils.bh(table.loc[ok, "p"])
    controls = table[table.is_control & (table.status == "ok")]
    table = table.join(versus(table[ok], lambda r: controls[controls.comparison == r.comparison].estimate))
    utils.write(table.drop(columns=["converged"]), OUT / "igf2bp2_striatum_x_age.tsv")
    print(table[~table.is_control][["source", "tier", "comparison", "n_targets", "slope_focal_pct", "slope_other_pct",
                                    "estimate_pct", "ci_low_pct", "ci_high_pct", "p", "q", "pctile_signed"]]
          .round(4).to_string(), flush=True)


def proteome_unit(unit, species):
    protein, _, samples, X = unit_design(unit)
    genes = gene_table(species, protein[protein.notna().all(axis=1)].mean(axis=1))
    Y, X = protein.loc[genes.index, samples.index].to_numpy(float), X.to_numpy(float)
    rows = []
    for source, tier, sets, own in tests(species):
        for comparison in COMPARISONS:
            w, nt, nc, _ = score_weights(genes, sets, comparison, own)
            row = dict(unit=unit, species=species, source=source, tier=tier, comparison=comparison, n_targets=nt,
                       n_comparison=nc, n=len(samples))
            rows.append(row | (ols_slope(w @ Y, X) if min(nt, nc) >= MIN_TARGETS else {}))
    return rows


def proteome_slopes():
    """Age slope of the target scores in the striatal proteomes, pooled by species (random effects)."""
    with mp.get_context("fork").Pool(4) as pool_:
        units = pd.DataFrame([r for rows in pool_.starmap(
            proteome_unit, [(u, sp) for sp, us in PROTEOMES.items() for u in us]) for r in rows])
    utils.write(percent(units), OUT / "igf2bp2_proteome_units.tsv")
    table = percent(pool(units, ["species", "source", "tier", "comparison"]))
    table["is_control"] = table.source.str.startswith("CTRL:")
    for species in ("human", "mouse"):
        m = (table.species == species) & ~table.is_control
        primary = m & (table.source == f"{FOCUS}:primary")
        table.loc[primary, "q"] = utils.bh(table.loc[primary, "p"])
        table.loc[m, "q_all"] = utils.bh(table.loc[m, "p"])
    controls = table[table.is_control]
    table = table.join(versus(table[~table.is_control],
                              lambda r: controls[(controls.species == r.species) & (controls.comparison == r.comparison)
                                                 & (controls.n_units == r.n_units)].estimate))
    utils.write(table, OUT / "igf2bp2_proteome_pooled.tsv")
    print(table[table.source.isin([f"{FOCUS}:primary", f"{FOCUS}:secondary"])]
          [["species", "source", "tier", "comparison", "n_units", "estimate_pct", "ci_low_pct", "ci_high_pct", "p", "q",
            "i2", "pctile_signed"]].round(4).to_string(), flush=True)


def igf2bp2_loss():
    """Log2 fold change of the targets minus matched non-targets in each IGF2BP2 knockdown or knockout; T3 without
    the dataset itself."""
    pert = (utils.read(COMP / "perturbation.tsv", low_memory=False)
            .query("reader == @FOCUS and direction in ['KD', 'KO']"))
    level = {"human": level_of("gtex_caudate"), "mouse": level_of("chdi2016_rna")}
    rng = np.random.default_rng(utils.SEED)
    rows = []
    for dataset, d in pert.groupby("dataset"):
        species = d.species.iloc[0]
        d = d.drop_duplicates("gene_id").set_index("gene_id")
        genes = gene_table(species, level[species])
        genes = genes.loc[genes.index.intersection(d.index)]
        lfc = d.log2fc.reindex(genes.index).to_numpy(float)
        for tier in ("T1", "T2", "T3"):
            target = (None if tier != "T3" else G["sets"]["primary"]["targets"]["T1"]
                      & set().union(*[v for k, v in G["t3_by_dataset"].items() if k != dataset]))
            for comparison in COMPARISONS:
                row = dict(dataset=dataset, species=species, tier=tier, comparison=comparison)
                rows.append(row | (perturbation_effect(lfc, genes, spec("primary", tier, species, target), comparison,
                                                       own_of(species, FOCUS), rng) or {}))
    table = pd.DataFrame(rows)
    ok = table.p.notna()
    table.loc[ok, "q"] = utils.bh(table.loc[ok, "p"])
    utils.write(table, OUT / "igf2bp2_perturbation.tsv")
    print(table[["dataset", "tier", "comparison", "n_targets", "n_comparison", "effect_pct", "ci_low", "ci_high", "p",
                 "q"]].round(4).to_string(), flush=True)


def nuclei_prepare(dataset):
    """Donor pseudobulks of the expressing populations (MSN, OPC; group 1) and the glia that do not express IGF2BP2
    (group 0), log2 CPM."""
    folder = utils.PREPARED / dataset
    samples = utils.read(folder / "pseudobulk_samples.tsv", low_memory=False, dtype={"sample_id": str})
    samples = samples[(samples.arm == "aging") & (samples.n_cells >= sn.MIN_NUCLEI)
                      & samples.cell_type.isin(POPULATIONS[dataset] + GLIA)].copy()
    samples["group"] = samples.cell_type.isin(POPULATIONS[dataset]).astype(float)
    samples["age"] = samples.age_years / 10
    counts = utils.read(folder / "pseudobulk.parquet", columns=samples.sample_id.tolist())
    all_log = np.log2((counts + 0.5).div(counts.sum() + 1, axis=1) * 1e6)
    cpm = counts.div(counts.sum(), axis=1) * 1e6
    expressed = pd.concat([(cpm.loc[:, samples.group.to_numpy() == g] >= 1).mean(axis=1) >= 0.5 for g in (0, 1)],
                          axis=1).all(axis=1)
    logcpm, all_log = sn.to_ensembl(all_log[expressed], folder), sn.to_ensembl(all_log, folder)
    frame = samples.rename(columns={"cell_type": "population"}).set_index("sample_id", drop=False)
    X = model.design(frame, {"numeric": [], "categorical": ["population", "sex"]
                             + (["batch"] if dataset.startswith("green") else [])})
    genes = gene_table("human", logcpm.mean(axis=1))
    return dict(frame=frame, X=X, age=X["age"].to_numpy(), X0=X.drop(columns="age"), donor=frame.donor_id.to_numpy(),
                focal=frame.group.to_numpy().astype(bool), genes=genes,
                Y=logcpm.loc[genes.index, frame.index].to_numpy(), all_log=all_log)


def nuclei_job(job):
    dataset, label, sets, own = job
    d = G["nuclei"][dataset]
    w, nt, nc, _ = score_weights(d["genes"], sets, label["comparison"], own)
    row = dict(cohort=dataset, n_targets=nt, n_comparison=nc, n_samples=len(d["frame"])) | label
    return row | (model.fit_slopes(w @ d["Y"], d["X0"], d["age"], d["donor"], d["focal"])
                  if min(nt, nc) >= MIN_TARGETS else {})


def nuclei_age():
    """Population group x age interaction of the target scores (mixed model with a donor intercept) in Green 2025
    and Handsaker 2025, combined by random effects."""
    G["nuclei"] = {d: nuclei_prepare(d) for d in POPULATIONS}
    sets = ([(f"{FOCUS}:primary", tier, spec("primary", tier, "human"), own_of("human", FOCUS))
             for tier in ("T1", "T2", "T3")]
            + [(f"CTRL:{r}", "T1", g, own_of("human", r)) for r, g in G["ctrl"]["human"].items()])
    jobs = [(d, dict(source=s, tier=t, comparison=c), g, own)
            for d in POPULATIONS for s, t, g, own in sets for c in COMPARISONS]
    with mp.get_context("fork").Pool(WORKERS) as pool_:
        cohorts = pd.DataFrame(pool_.map(nuclei_job, jobs, chunksize=1))
    gene = nw.symbols["human"][FOCUS]  # IGF2BP2 mRNA in the same model
    own = [dict(cohort=d, source=f"{FOCUS} mRNA", tier="gene", comparison="-", n_samples=len(v["frame"]))
           | model.fit_slopes(v["all_log"].loc[gene, v["frame"].index].to_numpy(), v["X0"], v["age"], v["donor"],
                              v["focal"])
           for d, v in G["nuclei"].items()]
    cohorts = percent(pd.concat([cohorts, pd.DataFrame(own)], ignore_index=True))
    utils.write(cohorts.drop(columns="converged"), OUT / "igf2bp2_nuclei_cohorts.tsv")
    table = percent(pool(cohorts, ["source", "tier", "comparison"], "cohort")).rename(columns={"n_units": "n_cohorts"})
    table["is_control"] = table.source.str.startswith("CTRL:")
    ok = ~table.is_control & (table.tier != "gene")
    table.loc[ok, "q"] = utils.bh(table.loc[ok, "p"])
    controls = table[table.is_control]
    table = table.join(versus(table[ok], lambda r: controls[(controls.comparison == r.comparison)
                                                            & (controls.n_cohorts == r.n_cohorts)].estimate))
    utils.write(table, OUT / "igf2bp2_nuclei_combined.tsv")
    print(table[~table.is_control][["source", "tier", "comparison", "n_cohorts", "estimate_pct", "ci_low_pct",
                                    "ci_high_pct", "p", "q", "i2", "pctile_signed"]].round(4).to_string(), flush=True)


def matched_score(logcpm, in_target, in_comparison):
    """Mean of the targets minus the comparison genes weighted to the targets over 5 x 5 bins of expression and
    length."""
    length = G["length"].reindex(logcpm.index).to_numpy()
    keep = (in_target | in_comparison) & ~np.isnan(length)
    values, target = logcpm[keep], in_target[keep]
    stratum = 5 * pd.qcut(values.mean(axis=1), 5, labels=False).to_numpy() + pd.qcut(length[keep], 5, labels=False)
    means = values[~target].groupby(stratum[~target]).mean()
    counts = pd.Series(stratum[target]).value_counts().reindex(means.index, fill_value=0)
    return values[target].mean().to_numpy() - (means.T @ counts).to_numpy() / counts.sum()


def neuron_glia_cohort(cohort):
    """Neurons against other nuclei (donor pseudobulks) for every set of G['neuron_glia_sets']."""
    spec_ = utils.SINGLE_NUCLEUS[cohort]
    folder = utils.PREPARED / cohort
    stem = "rna" if (folder / "rna.parquet").exists() else "pseudobulk"
    samples = utils.read(folder / f"{stem}_samples.tsv", dtype={"sample_id": str}).query(spec_["keep"])
    neurons, others = spec_["neuron"].split(","), spec_["other"].split(",")
    samples = samples[samples.cell_type.isin(neurons + others)]
    summed = (utils.read(folder / f"{stem}.parquet")[samples.sample_id].T
              .groupby([samples.donor_id.to_numpy(), samples.cell_type.isin(neurons).to_numpy()]).sum().T)
    library = summed.sum()
    expressed = ((summed.div(library, axis=1) * 1e6 >= 1).T.groupby(level=1).mean() >= 0.5).all()
    logcpm = np.log2((summed[expressed.to_numpy()] + 0.5).div(library + 1, axis=1) * 1e6)
    logcpm = sn.to_ensembl(logcpm, folder, stem)
    donors = samples.drop_duplicates("donor_id").set_index("donor_id")
    donor, neuron = (np.array(level) for level in zip(*logcpm.columns, strict=True))
    frame = pd.DataFrame({"neuron": neuron.astype(float), "sex": donors.sex.loc[donor].to_numpy()})
    frame["age"] = donors.age_years.loc[donor].to_numpy() / 10
    X = model.design(frame, {"numeric": ["neuron"], "categorical": ["sex"]})
    age = X.pop("age").to_numpy()
    rows = []
    for name, (target, comparison) in G["neuron_glia_sets"].items():
        in_t, in_c = logcpm.index.isin(target), logcpm.index.isin(comparison)
        if min(in_t.sum(), in_c.sum()) < MIN_TARGETS:
            continue
        score = {"all": logcpm[in_t].mean().to_numpy() - logcpm[in_c].mean().to_numpy(),
                 "matched": matched_score(logcpm, in_t, in_c)}
        for kind, values in score.items():
            rows.append(dict(cohort=cohort, set=name, comparison=kind, n_target=in_t.sum(), n_comparison=in_c.sum())
                        | model.fit_slopes(values, X, age, donor, neuron))
    return rows


def neuron_glia():
    """Neuron-versus-glia age slope of the IGF2BP2 targets, against the control RBPs' targets built the same way."""
    atlas = nw.atlas[(nw.atlas.species == "Homo sapiens")
                     & (nw.atlas.resolution == "single_base")].reset_index(drop=True)
    modified = set(atlas.gene_id)
    members = utils.read(utils.OUT / "40_target_sets" / "target_genes.tsv.gz").query(
        "species == 'human' and modification == 'm6A' and regulator == @FOCUS and window == @utils.PEAK_WINDOW")
    sets = {f"{FOCUS} (CLIP within 50 nt)": (set(members.gene_id[members.group == "target"]),
                                             set(members.gene_id[members.group == "comparison"]))}
    for tier in ("T1", "T2", "T3"):
        s = spec("primary", tier, "human")
        sets[f"{FOCUS} primary {tier} vs unbound m6A"] = (s.target, modified - set(s.bound))
        sets[f"{FOCUS} primary {tier} vs bound-unmarked"] = (s.target, s.unmarked)
    for rbp in m33.CONTROLS:
        target = set(atlas.gene_id[utils.near_peaks(atlas, nw.idr_peaks(rbp), utils.PEAK_WINDOW)])
        sets[f"CTRL:{rbp}"] = (target, modified - target)
    G.update(neuron_glia_sets=sets, length=nw.features["human"].log_gene_len)
    with mp.get_context("fork").Pool(len(utils.SINGLE_NUCLEUS)) as pool_:
        cohorts = pd.DataFrame([r for rows in pool_.map(neuron_glia_cohort, list(utils.SINGLE_NUCLEUS)) for r in rows])
    utils.write(cohorts.drop(columns="converged"), OUT / "igf2bp2_neuron_glia_cohorts.tsv")
    table = pd.DataFrame([dict(zip(["set", "comparison"], key, strict=True)) | utils.random_effects(g.estimate, g.se)
                          for key, g in cohorts.groupby(["set", "comparison"])]
                         ).rename(columns={"est": "estimate", "k": "n_cohorts"})
    table["estimate_pct"] = utils.pct(table.estimate)
    table["is_control"] = table["set"].str.startswith("CTRL:")
    controls = table[table.is_control]
    focus = table[~table.is_control]
    stats = pd.DataFrame([dict(n_controls=len(ref), pctile_z=100 * (ref.z < r.z).mean(),
                               p_control_two_sided=(1 + (ref.z.abs() >= abs(r.z)).sum()) / (1 + len(ref)))
                          for r in focus.itertuples()
                          for ref in [controls[(controls.comparison == r.comparison)
                                               & (controls.n_cohorts == r.n_cohorts)]]], index=focus.index)
    table = table.join(stats)
    utils.write(table, OUT / "igf2bp2_neuron_glia_combined.tsv")
    print(table[~table.is_control][["set", "comparison", "n_cohorts", "estimate_pct", "z", "p", "n_controls",
                                    "pctile_z", "p_control_two_sided"]].round(4).to_string(), flush=True)


def coupling_unit(unit):
    """Coupling of IGF2BP2 mRNA with its primary T1 targets (script 33 model, against the control RBPs) and mediation
    of the target score's age slope."""
    rng = np.random.default_rng(utils.SEED)
    controls = list(G["ctrl"]["human"])
    data = sn.load_unit(unit, controls)
    if FOCUS not in data["regulators"]:
        return pd.DataFrame(), []
    items = [nw.Item(FOCUS, "primary_T1", True, G["s_focus"], G["u_focus"], [FOCUS], "reader")]
    items += [nw.Item(r, "ENCODE_IDR", True, G["idr"][r][0], G["idr"][r][1], [r], "control_rbp")
              for r in controls if r in data["regulators"]]
    load = np.log1p(data["load"].reindex(data["genes"].index)).rename("load")
    perm = nw.permutations(data["strata"], rng)
    edges = nw.run_edges(data, nw.Stage1(data), perm, items, nw.stage2_fits(data, items, load))
    edges = edges.assign(population=unit["population"], dataset=unit["dataset"])
    genes = data["genes"].assign(load=data["genes"][nw.LOAD["human"]].astype(float))
    X0, age = data["X"][:, :-1], data["X"][:, -1]
    R = nw.regulator_matrix(data, items[:1])
    mediation = []
    for tier in ("T1", "T2", "T3"):
        for comparison in COMPARISONS:
            w, nt, nc, _ = score_weights(genes, spec("primary", tier, "human"), comparison, own_of("human", FOCUS))
            row = dict(unit=unit["unit"], dataset=unit["dataset"], population=unit["population"], tier=tier,
                       comparison=comparison, n_targets=nt, n_comparison=nc, n=len(age))
            if min(nt, nc) >= MIN_TARGETS:
                y = w @ data["Y"]
                without, with_r = m33.age_slopes(y, None, X0, age, R)
                boot = np.array([1 - np.divide(*m33.age_slopes(y[i], None, X0[i], age[i], R[i])[::-1])
                                 for i in nw.resamples(data["strata"], NBOOT, rng)])[:, 0]
                row |= dict(slope_without_pct=utils.pct(without), slope_with_pct=utils.pct(with_r[0]),
                            fraction=1 - with_r[0] / without, ci_low=np.percentile(boot, 2.5),
                            ci_high=np.percentile(boot, 97.5))
            mediation.append(row)
    return edges, mediation


def coupling_mediation():
    G["s_focus"], G["u_focus"] = s_u("human", spec("primary", "T1", "human"))
    load = nw.site_map("human")["load"]
    units = [dict(u, load=load) for u in sn.find_units()
             if u["dataset"] in POPULATIONS and u["population"] in POPULATIONS[u["dataset"]]]
    with mp.get_context("fork").Pool(min(WORKERS, len(units))) as pool_:
        res = pool_.map(coupling_unit, units, chunksize=1)
    edges = pd.concat([e for e, _ in res if len(e)], ignore_index=True)
    mediation = pd.DataFrame([r for _, m in res for r in m])
    ctl = edges[edges.family == "control_rbp"]
    coupling = edges[edges.regulator == FOCUS].reset_index(drop=True)

    def ref(r):
        return ctl[ctl.unit == r.unit].b_s

    coupling = coupling.join(versus(coupling, ref, "b_s")).assign(
        control_mean_b_s=[ref(r).mean() for r in coupling.itertuples()],
        control_sd_b_s=[ref(r).std() for r in coupling.itertuples()])
    coupling["z_control"] = (coupling.b_s - coupling.control_mean_b_s) / coupling.control_sd_b_s
    coupling["ci_low"] = coupling.b_s - utils.Z95 * coupling.se_s
    coupling["ci_high"] = coupling.b_s + utils.Z95 * coupling.se_s
    utils.write(ctl, OUT / "igf2bp2_control_edges.tsv")
    utils.write(coupling, OUT / "igf2bp2_coupling.tsv")
    utils.write(mediation, OUT / "igf2bp2_mediation.tsv")
    print(coupling[["unit", "n", "n_genes_s_pos", "b_s", "se_s", "p_s", "b_u", "n_controls", "z_control",
                    "pctile_signed"]].round(4).to_string(), flush=True)
    print(mediation[["population", "dataset", "tier", "comparison", "n_targets", "slope_without_pct", "slope_with_pct",
                     "fraction", "ci_low", "ci_high"]].round(4).to_string(), flush=True)


# YTHDF1-3
def binding_sets():
    """(species, reader, source) -> target set (T1 and bound-unmarked genes of script 37) and whether it is the
    reader's preferred source."""
    table = utils.read(COMP / "yth_binding_sets.tsv", low_memory=False)
    table = table[table.reader.isin(READERS3)]
    return {key: (Genes(set(g.gene_id[g.tier == "T1"]), set(g.gene_id), set(g.gene_id[g.tier == "bound_unmarked"])),
                  bool(g.primary.any())) for key, g in table.groupby(["species", "reader", "source"])}


def yth_job(mode, unit, layer):
    """Age slope of the target scores, coupling of YTHDF3 and the YTHDF module with the bound sets or with the m6A
    load, and mediation of the age slope, in one unit and layer."""
    rng = np.random.default_rng(utils.SEED)
    data = m33.load_data(unit, layer, wildtype=mode == "wildtype")
    species, genes = data["species"], data["genes"]
    arm = nw.ARMS[layer]
    gt = pd.DataFrame(dict(level=genes.level, log_gene_len=genes.log_gene_len,
                           load=genes[nw.LOAD[species]].astype(float)), index=genes.index)
    X, Y, M = data["X"], data["Y"], data["M"]
    age, X0 = X[:, -1], X[:, :-1]
    named = [(r, source, sets, primary, "reader")
             for (sp, r, source), (sets, primary) in G["binding_sets"].items() if sp == species]
    named += [(r, "ENCODE_IDR", sets, True, "control") for r, sets in G["ctrl"][species].items()]
    y1 = []
    for reader, source, sets, primary, kind in named:
        for comparison in COMPARISONS:
            w, nt, nc, _ = score_weights(gt, sets, comparison, own_of(species, reader))
            row = dict(mode=mode, unit=unit, species=species, layer=layer, arm=arm, kind=kind, reader=reader,
                       source=source, primary=primary, comparison=comparison, n_targets=nt, n_comparison=nc, n=len(age))
            y1.append(row | (ols_slope(w @ Y, X, None if M is None else w @ M) if min(nt, nc) >= MIN_TARGETS else {}))
    # items of the coupling model: the YTHDF module and YTHDF3 against the m6A load and against the preferred bound
    # sets, the control RBPs in both
    load = np.log1p(data["load"].reindex(genes.index)).rename("load")
    have = data["regulators"]
    items, meta = [], []
    members = [m for m in READERS3 if m in have]
    on_load = [("YTHDF", members, "module")] if members else []
    on_load += [("YTHDF3", ["YTHDF3"], "reader")] if "YTHDF3" in have else []
    for name, group, family in on_load + [(c, [c], "control_rbp") for c in m33.CONTROLS if c in have]:
        items.append(nw.Item(name, "m6A_load", True, load, None, group, family))
        meta.append(dict(target_type="m6A_load", target_reader="m6A_load", target_source="m6A_load"))
    for reader, source, sets, primary, kind in named:
        if (not primary or len(set(genes.index) & sets.target) < MIN_TARGETS
                or (kind == "control" and reader not in have)):
            continue
        s, u = s_u(species, sets)
        regulators = ([(reader, [reader], "control_rbp")] if kind == "control"
                      else [("YTHDF3", ["YTHDF3"], "reader"), ("YTHDF", members, "module"),
                            (reader, [reader], "reader")])
        seen = set()
        for name, group, family in regulators:
            if name in seen or not group or not all(m in have for m in group):
                continue
            seen.add(name)
            items.append(nw.Item(name, source if kind == "control" else f"{reader}|{source}", True, s, u, group,
                                 family))
            meta.append(dict(target_type="bound_set", target_reader="control" if kind == "control" else reader,
                             target_source=source))
    perm = nw.permutations(data["strata"], rng)
    edges = nw.run_edges(data, nw.Stage1(data), perm, items, nw.stage2_fits(data, items, load))
    edges = pd.concat([edges.reset_index(drop=True), pd.DataFrame(meta)], axis=1).assign(mode=mode)
    # mediation of the age slope of the m6A-load coefficient (m6A load) or of each target score (bound sets)
    R = nw.regulator_matrix(data, items)
    boots = nw.resamples(data["strata"], NBOOT, rng)
    rows_all = np.arange(len(age))
    D = nw.stage2_design(data, load, None, load, np.ones(len(genes), bool))
    parts = [(np.linalg.pinv(D), D[:, :2])]
    load_inverse = np.linalg.pinv(np.column_stack([np.ones(len(genes)), load.to_numpy(), data["covariates"]]))
    y3 = []
    for k, (item, mt) in enumerate(zip(items, meta, strict=True)):
        if item.family == "control_rbp":
            continue
        base = dict(mode=mode, unit=unit, species=species, layer=layer, arm=arm, target_type=mt["target_type"],
                    regulator=item.name,
                    target_reader=mt["target_reader"], target_source=mt["target_source"])
        if mt["target_type"] == "m6A_load":

            def evaluate(rows, r):
                return m33.mediate(Y[:, rows], None if M is None else M[:, rows], X0[rows], age[rows], r, parts,
                                   load_inverse)

            a0, a1 = evaluate(rows_all, R[:, [k]])
            f0 = 1 - a1 / a0
            boot = np.array([1 - np.divide(*evaluate(b, R[b][:, [k]])[::-1]) for b in boots])
            null = np.array([1 - np.divide(*evaluate(rows_all, R[perm[:, p]][:, [k]])[::-1]) for p in range(nw.NPERM)])
            y3.append(base | dict(comparison="m6A_load", n_targets=len(genes), slope_without_pct=utils.pct(a0),
                                  slope_with_pct=utils.pct(a1), fraction=f0, ci_low=np.percentile(boot, 2.5),
                                  ci_high=np.percentile(boot, 97.5), null_mean=null.mean(), null_sd=null.std(ddof=1),
                                  p_perm=nw.twosided(f0 - null.mean(), null - null.mean())))
            continue
        sets = G["binding_sets"][species, mt["target_reader"], mt["target_source"]][0]
        for comparison in COMPARISONS:
            w, nt, nc, _ = score_weights(gt, sets, comparison, own_of(species, mt["target_reader"])
                                         + [data["own"][m] for m in item.members])
            if min(nt, nc) < MIN_TARGETS:
                continue
            y, m = w @ Y, None if M is None else w @ M
            without, with_r = m33.age_slopes(y, m, X0, age, R[:, [k]])
            f0 = 1 - with_r[0] / without

            def fraction(rows, r, y=y, m=m):
                a, b = m33.age_slopes(y[rows], None if m is None else m[rows], X0[rows], age[rows], r)
                return 1 - b / a

            boot = np.array([fraction(b, R[b][:, [k]])[0] for b in boots])
            boot = boot[np.isfinite(boot)]
            null = fraction(rows_all, R[perm, k])
            y3.append(base | dict(comparison=comparison, n_targets=nt, n_comparison=nc,
                                  slope_without_pct=utils.pct(without), slope_with_pct=utils.pct(with_r[0]),
                                  fraction=f0, ci_low=np.percentile(boot, 2.5), ci_high=np.percentile(boot, 97.5),
                                  null_mean=null.mean(), null_sd=null.std(ddof=1),
                                  p_perm=nw.twosided(f0 - null.mean(), null - null.mean())))
    print(mode, unit, layer, "done", len(genes), "genes", flush=True)
    return pd.DataFrame(y1), edges, pd.DataFrame(y3)


def ythdf3_mrna():
    """Age slope of YTHDF3 mRNA in every human bulk unit with an mRNA layer (OLS, HC3, the unit's design); 'falls' =
    negative with BH q < 0.05 over the units."""
    units = utils.read(utils.PREPARED / "units.tsv", index_col="unit")
    gene = nw.symbols["human"]["YTHDF3"]
    rows = []
    for unit in units.index:
        if not (utils.PREPARED / unit / "rna.parquet").exists():
            continue
        _, rna, samples, X = unit_design(unit)
        if gene not in rna.index:
            continue
        y = rna.loc[gene, samples.index].to_numpy(float)
        ok = ~np.isnan(y)
        rows.append(dict(unit=unit, dataset=units.dataset[unit], n=int(ok.sum()))
                    | ols_slope(y[ok], X.to_numpy(float)[ok]))
    table = percent(pd.DataFrame(rows))
    table["q"] = utils.bh(table.p)
    table["falls"] = (table.estimate < 0) & (table.q < utils.Q)
    return table


def yth_summary(y1, edges, y3):
    """Pool the units by species (random effects), compare with the control RBPs and apply BH within each test family of
    one set of units."""
    y1 = percent(y1)
    ok = y1.p.notna() & (y1.kind == "reader")
    y1.loc[ok, "q_unit"] = y1[ok].groupby("arm").p.transform(utils.bh)

    def pooled_controls(ctrl, match):
        """Each control RBP pooled over the units of the tested row where it is available (at least 75% of them)."""
        def ref(r):
            units = r.units.split(",")
            c = ctrl[match(ctrl, r) & ctrl.unit.isin(units)]
            return [utils.random_effects(g.estimate, g.se)["est"] if len(g) > 1 else g.estimate.iloc[0]
                    for _, g in c.groupby("control") if len(g) >= 0.75 * len(units)]
        return ref

    keys = ["species", "arm", "kind", "reader", "source", "primary", "comparison"]
    p1 = percent(pool(y1, keys))
    c1 = y1[y1.kind == "control"].rename(columns={"reader": "control"}).dropna(subset=["estimate", "se"])
    p1 = p1[p1.kind == "reader"].reset_index(drop=True)
    p1 = p1.join(versus(p1, pooled_controls(c1, lambda c, r: (c.species == r.species) & (c.arm == r.arm)
                                            & (c.comparison == r.comparison))))
    primary = p1.primary & p1.p.notna()
    p1.loc[primary, "q_arm"] = p1[primary].groupby("arm").p.transform(utils.bh)  # preferred sets, BH within arm
    p1["q_all"] = utils.bh(p1.p)
    # coupling
    edges = edges.assign(is_control=edges.family == "control_rbp")
    tested = ~edges.is_control
    edges.loc[tested, "q_unit"] = edges[tested].groupby(["target_type", "arm"]).p_s.transform(utils.bh)
    ctrl = edges[edges.is_control]
    edges = edges.join(versus(edges[~edges.is_control],
                              lambda r: ctrl[(ctrl.unit == r.unit) & (ctrl.target_type == r.target_type)].b_s, "b_s"))
    ee = edges.assign(estimate=edges.b_s, se=edges.se_s)
    c2 = ee[ee.is_control].rename(columns={"regulator": "control"}).dropna(subset=["estimate", "se"])
    p2 = pool(ee[~ee.is_control], ["species", "arm", "target_type", "regulator", "target_reader", "target_source"])
    p2 = p2.join(versus(p2, pooled_controls(c2, lambda c, r: (c.species == r.species) & (c.arm == r.arm)
                                            & (c.target_type == r.target_type))))
    p2["q"] = p2.groupby(["target_type", "arm"]).p.transform(utils.bh)
    # mediation
    y3 = y3.assign(q=y3.groupby(["target_type", "arm"]).p_perm.transform(utils.bh),
                   ci_excludes_0=lambda t: (t.ci_low > 0) | (t.ci_high < 0))
    p3 = y3.groupby(["species", "arm", "target_type", "regulator", "target_reader", "target_source", "comparison"]).agg(
        n_units=("unit", "size"), mean_fraction=("fraction", "mean"), median_fraction=("fraction", "median"),
        min_fraction=("fraction", "min"), max_fraction=("fraction", "max"), n_ci_excl0=("ci_excludes_0", "sum"),
        min_q=("q", "min")).reset_index()
    return y1, p1, edges, p2, y3, p3


def predicted_sign(readout, reader):
    """Sign of the target score after YTHDF loss: the translation efficiency or protein of YTHDF1/3 targets falls, the
    mRNA of YTHDF2 targets rises."""
    if readout in ("TE", "protein"):
        return -1 if reader in ("YTHDF1", "YTHDF3") else 0
    return 1 if readout == "RNA" and reader == "YTHDF2" else 0


def yth_loss():
    """Log2 fold change of each reader's m6A-bound targets minus matched non-targets after YTHDF loss or gain, per
    dataset and readout, against the predicted sign (loss lowers the translation efficiency or protein of YTHDF1/3
    targets and raises the mRNA of YTHDF2 targets)."""
    loss = ["KD", "KO", "cKO", "inhibitor", "triple KO"]
    pert = (utils.read(COMP / "perturbation.tsv", low_memory=False)
            .query("reader in @READERS3 or reader == 'YTHDF1/2/3'"))
    level = {"human": level_of("gtex_caudate"), "mouse": level_of("chdi2016_rna")}
    rng = np.random.default_rng(utils.SEED)
    rows = []
    for (dataset, reader, species, direction, readout), d in pert.groupby(
            ["dataset", "reader", "species", "direction", "readout"]):
        d = d.drop_duplicates("gene_id").set_index("gene_id")
        genes = gene_table(species, level[species])
        genes = genes[genes.index.isin(d.index)]
        lfc = d.log2fc.reindex(genes.index).to_numpy(float)
        genes, lfc = genes[np.isfinite(lfc)], lfc[np.isfinite(lfc)]
        for (sp, target_reader, source), (sets, primary) in G["binding_sets"].items():
            if sp != species or (reader != "YTHDF1/2/3" and target_reader != reader):
                continue
            sign = predicted_sign(readout, target_reader)
            for comparison in COMPARISONS:
                row = dict(dataset=dataset, reader=reader, target_reader=target_reader, source=source, primary=primary,
                           species=species, direction=direction, readout=readout, comparison=comparison,
                           n_genes=len(genes), predicted=sign if direction in loss else -sign)
                res = perturbation_effect(lfc, genes, sets, comparison, own_of(species, target_reader), rng)
                if res:
                    ref = [r["log2fc"] for c in G["ctrl"][species].values()
                           if (r := perturbation_effect(lfc, genes, c, comparison, own_of(species, target_reader), rng,
                                                        boot=False))]
                    row |= res | dict(n_controls=len(ref),
                                      pctile_signed=100 * (np.array(ref) < res["log2fc"]).mean() if ref else np.nan,
                                      sign_matches=(bool(np.sign(res["log2fc"]) == row["predicted"])
                                                    if row["predicted"] else np.nan))
                rows.append(row)
    table = pd.DataFrame(rows)
    ok = table.p.notna()
    table.loc[ok, "q"] = utils.bh(table.loc[ok, "p"])
    predicted = ok & (table.predicted != 0)
    table.loc[predicted, "q_predicted"] = utils.bh(table.loc[predicted, "p"])
    utils.write(table, OUT / "yth_perturbation.tsv")
    p = table[predicted & table.primary].assign(
        cls=lambda t: np.where(t.readout == "RNA", "RNA (YTHDF2)", "TE or protein (YTHDF1/3)"))
    p["matches_q05"] = (p.sign_matches == True) & (p.q_predicted < utils.Q)  # noqa: E712
    p["opposite_q05"] = (p.sign_matches == False) & (p.q_predicted < utils.Q)  # noqa: E712
    breakdown = p.groupby(["cls", "comparison"]).agg(
        n=("p", "size"), sign_matches=("sign_matches", "sum"), significant_match=("matches_q05", "sum"),
        significant_opposite=("opposite_q05", "sum")).reset_index()
    utils.write(breakdown, OUT / "yth_perturbation_breakdown.tsv")
    print(table[ok & table.primary][["dataset", "reader", "target_reader", "source", "direction", "readout",
                                     "comparison", "n_targets", "effect_pct", "ci_low", "ci_high", "p", "q",
                                     "predicted", "sign_matches", "pctile_signed"]].round(4).to_string())
    print(breakdown.to_string(), flush=True)


def yth():
    G["binding_sets"] = binding_sets()
    falls = ythdf3_mrna()
    utils.write(falls, OUT / "yth_ythdf3_mrna.tsv")
    units = utils.read(utils.PREPARED / "units.tsv", index_col="unit")
    protein_units = [u for u in units.index if (utils.PREPARED / u / "protein.parquet").exists()]
    jobs = [("main", u, "mrna_adjusted") for u in protein_units if (utils.PREPARED / u / "rna.parquet").exists()]
    jobs += ([("main", u, "protein_only") for u in protein_units]
             + [("main", u, "rna_only") for u in m33.STABILITY_UNITS])
    jobs += [("ythdf3_falls", u, "rna_only") for u in falls.unit[falls.falls]]  # human units where YTHDF3 mRNA falls
    print(len(jobs), "jobs", flush=True)
    with mp.get_context("fork").Pool(WORKERS) as pool_:
        results = pool_.starmap(yth_job, jobs, chunksize=1)
    tables = [[] for _ in range(6)]
    for mode in ("main", "ythdf3_falls"):
        parts = [pd.concat([r[i] for r, j in zip(results, jobs, strict=True) if j[0] == mode], ignore_index=True)
                 for i in range(3)]
        for table, new in zip(tables, yth_summary(*parts), strict=True):
            table.append(new.assign(mode=mode))
    names = ["target_slopes_units", "target_slopes_pooled", "coupling_edges", "coupling_pooled", "mediation_units",
             "mediation_summary"]
    for name, table in zip(names, tables, strict=True):
        utils.write(pd.concat(table, ignore_index=True), OUT / f"yth_{name}.tsv")
    p1 = tables[1][0]
    print(p1[p1.primary][["species", "arm", "reader", "source", "comparison", "n_units", "estimate_pct", "ci_low_pct",
                          "ci_high_pct", "p", "q_arm", "pctile_signed"]].round(4).to_string(), flush=True)
    yth_loss()


def wildtype():
    # the striatal mRNA unit of the Huntington's disease allelic series refitted on its wild-type mice: coupling with
    # the m6A load against the control RBPs of the unit, and the mediation of the age slope
    G["binding_sets"] = binding_sets()
    _, edges, y3 = yth_job("wildtype", m33.STABILITY_UNITS[0], "rna_only")
    edges = edges[edges.target_type == "m6A_load"]
    ctrl = edges[edges.family == "control_rbp"]
    edges = edges[edges.family != "control_rbp"]
    utils.write(edges.join(versus(edges, lambda r: ctrl.b_s, "b_s")), OUT / "yth_wildtype_coupling.tsv")
    utils.write(y3[y3.target_type == "m6A_load"], OUT / "yth_wildtype_mediation.tsv")


if __name__ == "__main__":
    part = sys.argv[1] if len(sys.argv) > 1 else "all"
    OUT.mkdir(parents=True, exist_ok=True)
    G["idr"], G["ctrl"] = control_sets()
    if part in ("all", "igf2bp2"):
        G["sets"] = igf2bp2_sets()
        for step in (striatum_age, proteome_slopes, igf2bp2_loss, nuclei_age, neuron_glia, coupling_mediation):
            step()
    if part in ("all", "yth"):
        yth()
    if part in ("all", "wildtype"):
        wildtype()
