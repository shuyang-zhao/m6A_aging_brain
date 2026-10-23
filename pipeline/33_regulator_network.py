# Regulator-target coupling network of m6A: in every unit, each regulator's level in an individual is added to the
# per-gene age model (stage 1, limma-moderated, one regulator or module at a time), and the moderated t of its
# coefficient is regressed over genes on s (m6A sites within 50 nt of the regulator's binding; human: m6A-seq peaks of
# the brain map), u (binding without m6A), the m6A load and the load-model covariates (stage 2); the coefficient of s
# is the edge, tested against permutations of the regulator across individuals within age strata. Writers and erasers
# have s = m6A load and no u; their edge must also exceed the largest edge of the negative controls METTL5, ZCCHC4 and
# METTL16 (methyltransferases of other RNAs) where these are quantified. Control RBPs (ENCODE eCLIP proteins not
# reported to read or depend on m6A) are fitted the same way (s, u from their IDR peaks; and the load edge) and set
# the null: an edge is m6A-specific only when it also exceeds the 95th percentile of the control RBPs' edges of the
# same kind (binding or load) in the same unit and layer.
# Arms: Output (protein), Translation (protein given the gene's own mRNA, paired cohorts), Stability (mRNA of the
# STABILITY_UNITS, none at present).
# Edges are BH-corrected within the regulator family of each unit and layer, shared across units with mashr
# (network_mash.R); the coupling core is in network.py
# reads scripts 04, 06, 07, 08 (units) and 09, data/raw/encode2020_eclip_idr (ENCODE IDR peaks) and the reader binding
# lists of ref/target_sets (network.py); writes output/33_regulator_network/
import os

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import model
import utils
from network import (
    COVARIATES,
    ERASERS,
    LOAD,
    NEGATIVE,
    READERS,
    WRITERS,
    Stage1,
    binding_links,
    classify,
    features,
    idr_files,
    idr_links,
    permutations,
    regulator_items,
    run_edges,
    site_map,
    stage2_fits,
    strata_of,
    symbols,
    to_mouse,
    with_q,
    z,
)

OUT = utils.OUT / "33_regulator_network"
HERE = Path(__file__).resolve().parent
GSL_LIBRARIES = utils.ROOT / "env" / "lib"  # libgsl and libopenblas that mashr was built against here
WORKERS = 24
# families that are not controls; the m6A regulators of this network are controls of no one
EXCLUDED = {"YTH reader": r"YTHD[CF]\d", "IGF2BP": r"IGF2BP\d", "FMRP family": r"FMR1|FXR[12]",
            "hnRNP C/G/A2B1": r"HNRNPCL?\d?|RBMX|HNRNPA2B1", "PRRC2": r"PRRC2[ABC]", "ELAVL1": r"ELAVL1",
            "G3BP": r"G3BP[12]", "SND1": r"SND1", "EIF3 subunit": r"EIF3.*", "SRSF3/7/10": r"SRSF(3|7|10)",
            "m6A writer or eraser of the network": "|".join(WRITERS + ERASERS)}
STABILITY_UNITS = []  # mouse mRNA units of the Stability arm; the only one (chdi2016_rna) is in utils.EXCLUDED
units = utils.read(utils.PREPARED / "units.tsv", index_col="unit").drop(utils.EXTERNAL.split())  # clock only


def control_rbps():
    """ENCODE eCLIP RBPs with the reason each one is not a control: an excluded family, or protein quantified (at least
    half of the samples) in fewer than half of the proteome units."""
    protein_units = [u for u in units.index if (utils.PREPARED / u / "protein.parquet").exists()]
    rbps = sorted(idr_files.rbp.unique())
    n = pd.Series(0, index=rbps)
    for unit in protein_units:
        protein = pd.read_parquet(utils.PREPARED / unit / "protein.parquet")
        species = "mouse" if protein.index[0].startswith("ENSMUSG") else "human"
        quantified = set(protein.index[protein.notna().mean(axis=1) >= 0.5])
        genes = [symbols[species].get(r if species == "human" else to_mouse.get(r)) for r in rbps]
        n += [g in quantified for g in genes]
    reason = []
    for rbp in rbps:
        excluded = [why for why, pattern in EXCLUDED.items() if re.fullmatch(pattern, rbp)]
        reason.append(excluded[0] if excluded
                      else f"protein quantified in {n[rbp]} of {len(protein_units)} proteome units"
                      if n[rbp] * 2 < len(protein_units) else "")
    cells = idr_files.groupby("rbp").cell.agg(lambda c: ",".join(sorted(set(c))))
    return pd.DataFrame(dict(rbp=rbps, cell_lines=cells.reindex(rbps).to_numpy(), n_units_quantified=n.to_numpy(),
                             n_units=len(protein_units), selected=[r == "" for r in reason], reason=reason))


control_table = control_rbps()
CONTROLS = list(control_table.rbp[control_table.selected])
IDR_READERS = [r for r in READERS if r in set(idr_files.rbp)]  # readers refitted on the same IDR peaks


# data of a unit and layer
def load_data(unit, layer):
    protein, rna, samples, design = utils.load_unit(unit)
    y, adjust = (rna, None) if layer == "rna_only" else (protein, rna if layer == "mrna_adjusted" else None)
    keep = y.index[y.notna().all(axis=1)]
    if adjust is not None:
        keep = keep.intersection(adjust.index[adjust.notna().all(axis=1)])
    species = "mouse" if keep[0].startswith("ENSMUSG") else "human"
    X = model.design(samples, design)
    k = design["hidden_factors"]
    if k:
        pcs = model.protected_pca(rna, X, k)
        X = pd.concat([X.drop(columns="age"), pcs.set_axis([f"PC{i + 1}" for i in range(k)], axis=1), X[["age"]]],
                      axis=1)
    level = ((rna if layer in ("rna_only", "mrna_adjusted") else protein)
             .loc[keep.intersection(features[species].index)].mean(axis=1))
    genes = features[species].join(level.rename("level"), how="inner")
    genes = genes[genes[COVARIATES].notna().all(axis=1)]
    # regulator levels: the protein of the individual, the mRNA where the protein is absent
    regulators, own = {}, {}
    for name in READERS + WRITERS + ERASERS + NEGATIVE + CONTROLS:
        gene = symbols[species].get(name)
        source = protein if layer != "rna_only" and protein is not None and gene in protein.index else rna
        if source is None or gene not in source.index:
            continue
        v = source.loc[gene, samples.index].to_numpy(float)
        v = np.where(np.isnan(v), np.nanmean(v), v)  # a missing level takes the mean of the unit
        regulators[name], own[name] = z(v), gene
    return dict(unit=unit, layer=layer, species=species, samples=samples, X=X.to_numpy(float), genes=genes,
                Y=y.loc[genes.index, samples.index].to_numpy(float),
                M=None if adjust is None else adjust.loc[genes.index, samples.index].to_numpy(float),
                load=genes[LOAD[species]].astype(float),
                covariates=np.column_stack([z(genes[c]) for c in COVARIATES]), regulators=regulators, own=own,
                strata=strata_of(samples))


def run_unit(unit, layer, S, U, tables):
    rng = np.random.default_rng(utils.SEED)
    data = load_data(unit, layer)
    species = data["species"]
    items, load = regulator_items(data, S[species], U[species], tables[species], CONTROLS)
    perm = permutations(data["strata"], rng)
    edges = run_edges(data, Stage1(data), perm, items, stage2_fits(data, items, load))
    print(unit, layer, len(data["genes"]), "genes", len(data["regulators"]), "regulators", flush=True)
    return edges


# sharing across units, classes and controls
def run_mashr(edges):
    """Posterior mean and lfsr of the edges (rows: regulator|source; columns: unit|layer)."""
    use = (edges[~edges.family.isin(["negative_control", "control_rbp"]) & (edges.se_s > 1e-8)]
           .assign(edge=lambda e: e.regulator + "|" + e.source, column=lambda e: e.unit + "|" + e.layer))
    bhat = use.pivot(index="edge", columns="column", values="b_s")
    shat = use.pivot(index="edge", columns="column", values="se_s")
    with tempfile.TemporaryDirectory() as folder:
        bhat.fillna(0).to_csv(f"{folder}/bhat.tsv", sep="\t")
        shat.fillna(100).to_csv(f"{folder}/shat.tsv", sep="\t")
        library_path = f"{GSL_LIBRARIES}:{os.environ.get('LD_LIBRARY_PATH', '')}"
        run = subprocess.run(["Rscript", str(HERE / "network_mash.R"), folder], capture_output=True, text=True,
                             env=os.environ | {"LD_LIBRARY_PATH": library_path})
        if run.returncode:
            raise RuntimeError(run.stderr)
        mash = pd.read_csv(f"{folder}/mash.tsv", sep="\t")
    mash[["unit", "layer"]] = mash.unit.str.split("|", expand=True)
    mash = mash.merge(use[["edge", "unit", "layer", "regulator", "source", "arm", "b_s", "se_s"]],
                      on=["edge", "unit", "layer"])
    mash["dataset"] = mash.unit.map(units.dataset)
    # shared edge: lfsr < 0.05 with the same sign in at least two independent cohorts (datasets)
    hits = (mash[mash.lfsr < utils.Q].assign(sign=lambda m: np.sign(m.posterior_mean))
            .groupby(["edge", "sign"]).dataset.nunique())
    mash["n_cohorts_same_sign"] = [hits.get((e, np.sign(p)), 0) if lfsr < utils.Q else 0
                                   for e, p, lfsr in zip(mash.edge, mash.posterior_mean, mash.lfsr, strict=True)]
    mash["shared_edge"] = mash.n_cohorts_same_sign >= 2
    return mash


def collect(results):
    edges = classify(with_q(pd.concat(results, ignore_index=True)))

    edges["edge"] = edges.regulator + "|" + edges.source
    mash = run_mashr(edges)
    edges = edges.merge(mash[["edge", "unit", "layer", "posterior_mean", "lfsr", "shared_edge"]],
                        on=["edge", "unit", "layer"], how="left")
    columns = ["unit", "species", "layer", "arm", "regulator", "family", "source", "primary", "n", "n_genes",
               "n_genes_s_pos", "b_s", "ci_low", "ci_high", "perm_P", "q", "b_u", "perm_P_u", "control", "b_control",
               "class_vs_comparator", "control_p95", "n_controls", "edge_class", "posterior_mean", "lfsr",
               "shared_edge"]
    short = ["unit", "species", "layer", "arm", "regulator", "n", "n_genes", "b_s", "ci_low", "ci_high", "perm_P", "q"]
    utils.write(edges[~edges.family.isin(["negative_control", "control_rbp"])][columns], OUT / "edges.tsv")
    utils.write(edges[edges.family == "control_rbp"][short + ["kind"]], OUT / "control_edges.tsv")
    utils.write(edges[edges.family == "negative_control"][short], OUT / "negative_controls.tsv")
    utils.write(mash, OUT / "mashr.tsv")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    utils.write(control_table, OUT / "control_rbps.tsv")
    sites = {s: site_map(s) for s in ("mouse", "human")}
    idr = dict(Parallel(n_jobs=WORKERS)(delayed(idr_links)(r, sites["human"]) for r in CONTROLS + IDR_READERS))
    built = {s: binding_links(s, sites[s], CONTROLS, idr) for s in ("mouse", "human")}
    S, U, tables = ({s: b[i] for s, b in built.items()} for i in range(3))
    utils.write(pd.concat(tables.values()), OUT / "sources.tsv")
    primary = tables["mouse"][tables["mouse"].primary & ~tables["mouse"].reader.isin(CONTROLS)]
    # genes with m6A sites in the primary binding of each reader (read by script 34)
    utils.write(pd.concat([S["mouse"][f"{r.reader}|{r.source}"].loc[lambda s: s > 0].rename("s").reset_index()
                           .assign(regulator=r.reader) for r in primary.itertuples()]), OUT / "mouse_targets.tsv.gz")
    protein_units = [u for u in units.index if (utils.PREPARED / u / "protein.parquet").exists()]
    jobs = [(u, "protein_only") for u in protein_units]
    jobs += [(u, "mrna_adjusted") for u in protein_units if (utils.PREPARED / u / "rna.parquet").exists()]
    jobs += [(u, "rna_only") for u in STABILITY_UNITS]
    results = Parallel(n_jobs=WORKERS)(delayed(run_unit)(u, layer, S, U, tables) for u, layer in jobs)
    collect(results)
