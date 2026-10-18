# Coupling core of the m6A regulator network, shared by scripts 33 (bulk units), 37, 38 and single_nucleus.py: s and u
# links of the regulators from binding and m6A maps, stage 1 (per-gene limma-moderated fits of age + one regulator),
# stage 2 (moderated t regressed over genes on s, u, m6A load and covariates), the permutation scheme,
# BH q-values and the control-RBP comparison.
# reads the CLIP peaks, site atlas and gene annotation (scripts 04, 06, 07, 09), data/raw/encode2020_eclip_idr and the
# reader binding lists of ref/target_sets; writes nothing (imported)
import os

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"
import subprocess
import tempfile
from collections import namedtuple
from pathlib import Path

import numpy as np
import pandas as pd

import utils

HERE = Path(__file__).resolve().parent
IDR = utils.DATA / "raw" / "encode2020_eclip_idr"  # ENCODE eCLIP IDR peaks (GRCh38; K562, HepG2, adrenal gland)
NPERM = utils.N_PERMUTATIONS
COVARIATES = utils.LOAD_COVARIATES
WRITERS = ["METTL3", "METTL14", "WTAP", "VIRMA", "ZC3H13", "RBM15", "RBM15B", "CBLL1"]
ERASERS = ["FTO", "ALKBH5"]
READERS = ["YTHDF1", "YTHDF2", "YTHDF3", "YTHDC1", "YTHDC2", "IGF2BP1", "IGF2BP2", "IGF2BP3", "FMR1", "FXR1", "FXR2",
           "HNRNPC", "RBMX", "PRRC2A"]  # HNRNPA2B1: no binding data; ALYREF is an m5C reader
NEGATIVE = ["METTL5", "ZCCHC4", "METTL16"]
MODULES = {
    "FMRP_family": ["FMR1", "FXR1", "FXR2"],
    "YTHDF": ["YTHDF1", "YTHDF2", "YTHDF3"],
    "IGF2BP": ["IGF2BP1", "IGF2BP2", "IGF2BP3"],
    "writer_complex": ["METTL3", "METTL14", "WTAP", "VIRMA"],
    "erasers": ERASERS,
}
ARMS = {"protein_only": "Output", "mrna_adjusted": "Translation", "rna_only": "Stability"}
MIN_BOUND_GENES = 500  # peak sets binding fewer genes rank last
LOAD = {"mouse": "n_sites_m6A", "human": "n_sites_m6A"}  # human: brain peaks (script 09)
Item = namedtuple("Item", "name source primary s u members family")

features = {
    species: utils.read(utils.OUT / "09_transcript_features" / f"{species}_features.tsv.gz").set_index("gene_id")
    for species in ("mouse", "human")
}


def family(name, controls):
    if name in controls:
        return "control_rbp"
    if name in NEGATIVE:
        return "negative_control"
    if name in WRITERS + ["writer_complex"]:
        return "writer"
    return "eraser" if name in ERASERS + ["erasers"] else "module" if name in MODULES else "reader"


def z(x):
    x = np.asarray(x, float)
    sd = x.std()
    return (x - x.mean()) / sd if sd > 0 else x * 0


# s and u: m6A sites within the binding of a reader, and binding without a site
symbols = {s: pd.Series(f.index, index=f.symbol.str.upper()) for s, f in features.items()}
symbols = {s: v[~v.index.duplicated(keep=False)] for s, v in symbols.items()}
orthologs = utils.orthologs()
to_human, to_mouse = (dict(zip(orthologs.mouse, orthologs.human, strict=True)),
                      dict(zip(orthologs.human, orthologs.mouse, strict=True)))
peaks = pd.read_parquet(utils.RESOURCES / "reader_clip" / "peaks.parquet")
exons = {s: pd.read_parquet(utils.RESOURCES / f"gencode_{s}" / "gene_exons.parquet") for s in ("mouse", "human")}
idr_files = utils.read(IDR / "files.tsv")

atlas = utils.read(utils.RESOURCES / "site_atlas" / "sites.parquet").query("modification == 'm6A'")
atlas = atlas[atlas.gene_id.notna() & atlas.tx_offset.notna() & atlas.ref_ok]
SITE_COLUMNS = ["chrom", "start", "end", "strand", "gene_id"]
UNIQUE = ["species", "modification", "chrom", "start", "end", "strand"]


def site_map(species):
    """Sites and m6A load of every gene: the two mouse GLORI maps, or the brain peaks of human (script 09 counts)."""
    sites = atlas.drop_duplicates(UNIQUE)
    sites = sites[sites.species == {"mouse": "Mus musculus", "human": "Homo sapiens"}[species]]
    load = features[species][LOAD[species]].astype(float)
    return dict(sites=sites[SITE_COLUMNS].reset_index(drop=True), load=load)


def near(a, b):
    """Rows of a (start, end; 1-based, closed) with an interval of b within the window."""
    out = np.zeros(len(a), bool)
    groups = dict(list(b.groupby(["chrom", "strand"])))
    for key, rows in a.groupby(["chrom", "strand"]):
        if key not in groups:
            continue
        other = groups[key].sort_values("start")
        starts, ends = other.start.to_numpy(), other.end.to_numpy()
        reach = np.maximum.accumulate(ends + utils.PEAK_WINDOW)
        last = np.searchsorted(starts - utils.PEAK_WINDOW, rows.end.to_numpy(), side="right") - 1
        out[a.index.get_indexer(rows.index)] = (last >= 0) & (reach[np.maximum(last, 0)] >= rows.start.to_numpy())
    return out


def peak_genes(peak, species):
    """(peak, gene) pairs of the peaks (BED, 0-based half-open) that overlap an exon (GTF, 1-based closed)."""
    pairs = []
    for (chrom, strand), rows in peak.groupby(["chrom", "strand"]):
        e = exons[species][(exons[species].chrom == chrom) & (exons[species].strand == strand)].sort_values("start")
        if e.empty:
            continue
        starts, ends, genes = e.start.to_numpy(), e.end.to_numpy(), e.gene_id.to_numpy()
        reach = np.maximum.accumulate(ends)
        a, b, index = rows.start.to_numpy(), rows.end.to_numpy(), rows.index.to_numpy()
        i = np.searchsorted(starts, b, side="right") - 1
        while True:
            ok = i >= 0
            ok[ok] = reach[i[ok]] > a[ok]
            if not ok.any():
                break
            j = i[ok]
            hit = ends[j] > a[ok]
            pairs.append(pd.DataFrame({"peak": index[ok][hit], "gene_id": genes[j[hit]]}))
            i = np.where(ok, i - 1, -1)
    return pd.concat(pairs) if pairs else pd.DataFrame(columns=["peak", "gene_id"])


def links_from_peaks(peak, species, site):
    peak = peak.drop_duplicates(["chrom", "start", "end", "strand"]).reset_index(drop=True)
    closed = peak.assign(start=peak.start + 1)  # BED to the 1-based closed intervals of the site atlas
    sites = site["sites"]
    s = sites[near(sites, closed)].groupby("gene_id").size()
    pairs = peak_genes(peak, species)
    unmarked = pairs[~near(closed, sites)[pairs.peak.to_numpy()]].drop_duplicates(["peak", "gene_id"])
    index = features[species].index
    return (np.log1p(s.reindex(index).fillna(0)), np.log1p(unmarked.groupby("gene_id").size().reindex(index).fillna(0)),
            set(pairs.gene_id))


def links_from_genes(gene_ids, species, site):
    # no positions: binding x m6A load of the listed genes; u marks listed genes without a site
    index = features[species].index
    listed = index.isin(gene_ids)
    load = site["load"].to_numpy(float)
    return (pd.Series(np.log1p(np.where(listed, load, 0)), index),
            pd.Series((listed & (load == 0)).astype(float), index), set(index[listed]))


def gene_ids(symbols_, native, species):
    x = pd.Series(list(symbols_)).dropna().str.upper()
    if native != species:
        x = x.map(to_human if native == "mouse" else to_mouse).dropna()
    return set(symbols[species].reindex(x).dropna())


def idr_peaks(rbp):
    """ENCODE IDR peaks of an RBP, all its experiments pooled."""
    names = ["chrom", "start", "end", "strand"]
    return pd.concat([pd.read_csv(IDR / f"{rbp}_{f.cell}_{f.file}.bed.gz", sep="\t", header=None, usecols=[0, 1, 2, 5],
                                  names=names) for f in idr_files[idr_files.rbp == rbp].itertuples()])


def idr_links(rbp, site):
    """s, u and bound genes from the pooled IDR peaks of an RBP (human; mouse goes by the bound orthologs)."""
    return rbp, links_from_peaks(idr_peaks(rbp), "human", site)


targets = utils.ROOT / "ref" / "target_sets"
tables = {f: pd.read_csv(targets / f, sep="\t") for f in [
    "darnell2011_fmrp.tsv", "maurin2018_fmrp.tsv", "wu2019_prrc2a.tsv", "stage2/shi2026_ythdf1_clip.tsv",
    "ascano2012_fmrp_fxr1_fxr2.tsv", "liu2015_hnrnpc.tsv", "liu2017_rbmx.tsv"]}
maurin = tables["maurin2018_fmrp.tsv"]
LISTS = {  # name: (species of the list, gene symbols)
    "Darnell2011_brain": ("mouse", tables["darnell2011_fmrp.tsv"].gene_symbol),
    "Maurin2018_brain": ("mouse", maurin[maurin.biotype == "mRNA"].gene_symbol),
    "Wu2019_PRRC2A_brain": ("mouse", tables["wu2019_prrc2a.tsv"].query("set.str.contains('NestinKO')").gene_symbol),
    "Shi2026_striatum_CLIP": ("mouse", tables["stage2/shi2026_ythdf1_clip.tsv"].query("arm == 'saline'").gene_symbol),
    "Ascano2012_FXR1": ("human",
                        tables["ascano2012_fmrp_fxr1_fxr2.tsv"].query("set.str.startswith('FXR1_')").gene_symbol),
    "Ascano2012_FXR2": ("human",
                        tables["ascano2012_fmrp_fxr1_fxr2.tsv"].query("set.str.startswith('FXR2_')").gene_symbol),
    "Liu2015_HNRNPC": ("human", tables["liu2015_hnrnpc.tsv"].query("set.str.contains('_HC_')").gene_symbol),
    "Liu2017_RBMX": ("human", tables["liu2017_rbmx.tsv"].gene_symbol),
}
BRAIN = {"FMR1": ["Darnell2011_brain", "Maurin2018_brain"], "PRRC2A": ["Wu2019_PRRC2A_brain"],
         "YTHDF1": ["Shi2026_striatum_CLIP"]}
CELL_LINE = {"FXR1": "Ascano2012_FXR1", "FXR2": "Ascano2012_FXR2", "HNRNPC": "Liu2015_HNRNPC", "RBMX": "Liu2017_RBMX"}


def binding_links(species, site, controls, idr):
    """s and u of every reader and source. Rank of a source: in vivo brain CLIP of the species 0, CLIP peaks of the
    species (ENCODE, PureCLIP) 1, ENCODE IDR peaks 1.2, CLIP gene lists of other cell types 2, lists and peaks mapped by
    orthology 3, 3.5 and 3.6; peak sets of fewer than MIN_BOUND_GENES genes 4. The lowest rank is the primary source.
    The control RBPs follow with their IDR peaks as the only source. idr: rbp -> (s, u, bound genes) of the human peaks;
    site: the m6A map of the species (peaks of the other species sit on its default map)."""
    other = "human" if species == "mouse" else "mouse"
    other_site = site_map(other)

    def idr_in(rbp):
        if species == "human":
            return idr[rbp]
        mapped = gene_ids(features["human"].symbol.reindex(list(idr[rbp][2])).dropna(), "human", "mouse")
        return links_from_genes(mapped, "mouse", site)

    S, U, rows = {}, {}, []
    for reader in READERS:
        sources = [("list", n, 0 if LISTS[n][0] == species else 3) for n in BRAIN.get(reader, [])]
        if reader in set(peaks[peaks.species == species].reader):
            sources.append(("peaks", reader, 1))
        if reader in idr:
            sources.append(("idr", reader, 1.2 if species == "human" else 3.6))
        if reader in CELL_LINE:
            sources.append(("list", CELL_LINE[reader], 2 if LISTS[CELL_LINE[reader]][0] == species else 3))
        if reader in set(peaks[peaks.species == other].reader):
            sources.append(("mapped_peaks", reader, 3.5))
        built = []
        for kind, name, rank in sources:
            if kind == "peaks":
                s, u, bound = links_from_peaks(
                    peaks[(peaks.species == species) & (peaks.reader == name)], species, site)
                rank = rank if len(bound) >= MIN_BOUND_GENES else 4
            elif kind == "mapped_peaks":
                _, _, bound = links_from_peaks(
                    peaks[(peaks.species == other) & (peaks.reader == name)], other, other_site)
                s, u, bound = links_from_genes(
                    gene_ids(features[other].symbol.reindex(list(bound)).dropna(), other, species), species, site)
                name += "_mapped"
            elif kind == "idr":
                s, u, bound = idr_in(name)
                name += "_IDR" if species == "human" else "_IDR_mapped"
            else:
                s, u, bound = links_from_genes(gene_ids(LISTS[name][1], LISTS[name][0], species), species, site)
                name += "" if LISTS[name][0] == species else "_mapped"
            built.append((rank, name, kind, s, u, len(bound)))
        for rank, name, kind, s, u, n_bound in built:
            S[f"{reader}|{name}"], U[f"{reader}|{name}"] = s, u
            rows.append(dict(species=species, reader=reader, source=name, kind=kind, rank=rank, n_bound_genes=n_bound,
                             primary=rank == min(b[0] for b in built), n_genes_s_pos=int((s > 0).sum()),
                             n_genes_u_pos=int((u > 0).sum())))
    for rbp in controls:
        s, u, bound = idr_in(rbp)
        S[f"{rbp}|ENCODE_IDR"], U[f"{rbp}|ENCODE_IDR"] = s, u
        rows.append(dict(species=species, reader=rbp, source="ENCODE_IDR", kind="idr", rank=1, n_bound_genes=len(bound),
                         primary=True, n_genes_s_pos=int((s > 0).sum()), n_genes_u_pos=int((u > 0).sum())))
    table = pd.DataFrame(rows)
    table["primary"] &= ~table.duplicated(["reader", "primary"])  # the first of equal-rank sources
    return pd.DataFrame(S), pd.DataFrame(U), table


# stage 1 (per gene) and stage 2 (across genes)
def strata_of(samples):
    # permutation blocks: age (bins of about ten donors when the ages are continuous), genotype, middle-age flags
    age = (pd.qcut(samples.age, max(2, len(samples) // 10), duplicates="drop") if samples.age.nunique() > 6
           else samples.age)
    parts = [age.astype(str)] + [samples[c].astype(str) for c in samples.columns
                                 if c == "genotype" or c.startswith("middle_")]
    return pd.Series(["|".join(t) for t in zip(*parts, strict=True)]).to_numpy()


def prior_of(variance, df):
    with tempfile.TemporaryDirectory() as folder:
        pd.DataFrame(dict(variance=variance, df=df)).to_csv(f"{folder}/variances.tsv", sep="\t", index=False)
        run = subprocess.run(["Rscript", str(HERE / "network_prior.R"), folder], capture_output=True, text=True)
        if run.returncode:
            raise RuntimeError(run.stderr)
        prior = pd.read_csv(f"{folder}/prior.tsv", sep="\t").iloc[0]
    return float(prior.var_prior), float(prior.df_prior)


class Stage1:
    """Per-gene fits of y ~ covariates + age + r (given the gene's own mRNA in the paired layer) for many r at once."""

    def __init__(self, data):
        Y, M, X = data["Y"], data["M"], data["X"]
        self.M = M
        self.H = np.eye(Y.shape[1]) - X @ np.linalg.pinv(X)
        self.Yr = Y @ self.H
        yy = (self.Yr**2).sum(1)
        if M is not None:
            self.Mr = M @ self.H
            self.mm = np.maximum((self.Mr**2).sum(1), 1e-12)
            self.ym = (self.Yr * self.Mr).sum(1)
            yy = yy - self.ym**2 / self.mm
        self.yy = yy
        self.df = Y.shape[1] - np.linalg.matrix_rank(X) - (M is not None) - 1
        # prior of the model without r
        self.var_prior, self.df_prior = prior_of(yy / (self.df + 1), np.full(len(yy), self.df + 1))

    def fit(self, R):
        """Coefficient and moderated t of every column of R (samples x regulators) in every gene."""
        Rr = self.H @ R
        ry, rr = self.Yr @ Rr, np.broadcast_to((Rr**2).sum(0), (len(self.yy), R.shape[1]))
        if self.M is not None:
            mr = self.Mr @ Rr
            ry = ry - self.ym[:, None] * mr / self.mm[:, None]
            rr = rr - mr**2 / self.mm[:, None]
        rr = np.maximum(rr, 1e-12)
        rss = np.maximum(self.yy[:, None] - ry**2 / rr, 0)
        variance = (self.df_prior * self.var_prior + rss) / (self.df_prior + self.df)
        return ry / rr, ry / rr / np.sqrt(variance / rr)


def stage2_design(data, focal, other, load, keep):
    """[1, focal, other, load, covariates] over the kept genes (no other and load for writers and erasers)."""
    index = data["genes"].index
    columns = [np.ones(keep.sum()), z(focal.reindex(index).fillna(0).to_numpy()[keep])]
    if other is not None:
        columns += [z(other.reindex(index).fillna(0).to_numpy()[keep]), z(load.to_numpy()[keep])]
    return np.column_stack(columns + [data["covariates"][keep]])


def regulator_items(data, S, U, table, controls):
    """Items (name, source, primary, s, u, members, family) of every regulator, module and control of the unit, and the
    log m6A load of its genes; s = load when u is None."""
    load = np.log1p(data["load"].reindex(data["genes"].index)).rename("load")
    have, primary = data["regulators"], table[table.primary].set_index("reader")
    fam = lambda name: family(name, controls)
    items = [Item(r.reader, r.source, r.primary, S[f"{r.reader}|{r.source}"], U[f"{r.reader}|{r.source}"], [r.reader],
                  fam(r.reader)) for r in table.itertuples() if r.reader in have]
    items += [Item(name, "m6A_load", True, load, None, [name], fam(name))
              for name in WRITERS + ERASERS + NEGATIVE + controls if name in have]
    for module, members in MODULES.items():
        present = [m for m in members if m in have]
        if module in ("writer_complex", "erasers"):
            if present:
                items.append(Item(module, "m6A_load", True, load, None, present, fam(module)))
            continue
        present = [m for m in present if m in primary.index]
        if present:
            s = pd.concat([S[f"{m}|{primary.source[m]}"] for m in present], axis=1).mean(axis=1)
            u = pd.concat([U[f"{m}|{primary.source[m]}"] for m in present], axis=1).mean(axis=1)
            items.append(Item(module, "primary_members", True, s, u, present, fam(module)))
    return items, load


def regulator_matrix(data, items):
    """Standardised level of every item (samples x items): the mean of its members' z-scores."""
    return np.column_stack([z(np.column_stack([data["regulators"][m] for m in item.members]).mean(axis=1))
                            for item in items])


def stage2_fits(data, items, load):
    """Kept genes (all but the regulator's own) and pseudo-inverse of the stage-2 design of every item."""
    fits = []
    for item in items:
        keep = ~data["genes"].index.isin([data["own"][m] for m in item.members])
        fits.append((keep, np.linalg.pinv(stage2_design(data, item.s, item.u, load, keep))))
    return fits


def twosided(observed, null):
    return (1 + (np.abs(null) >= abs(observed)).sum()) / (1 + len(null))


def run_edges(data, stage1, perm, items, fits):
    """Edges (coupling of r with age effects) of every regulator, module and control."""
    rows = []
    R = regulator_matrix(data, items)
    for k, (item, (keep, P)) in enumerate(zip(items, fits, strict=True)):
        Rk = np.column_stack([R[:, k], R[perm, k]])
        # rows: 1, s[, u, load], covariates; columns: observed, permutations
        coefficients = P @ stage1.fit(Rk)[1][keep]
        row = dict(unit=data["unit"], species=data["species"], layer=data["layer"], arm=ARMS[data["layer"]],
                   regulator=item.name, family=item.family, source=item.source, primary=item.primary,
                   n=len(data["samples"]), n_genes=int(keep.sum()),
                   n_genes_s_pos=int((item.s.reindex(data["genes"].index)[keep] > 0).sum()))
        for tag, i in (("s", 1), ("u", 2)):
            if i == 2 and item.u is None:
                continue
            observed, null = coefficients[i, 0], coefficients[i, 1:]
            row |= {f"b_{tag}": observed, f"se_{tag}": null.std(ddof=1), f"p_{tag}": twosided(observed, null)}
        if row["n_genes_s_pos"] == 0:  # no gene of the unit carries s: the edge is not defined
            row |= dict(b_s=np.nan, se_s=np.nan, p_s=np.nan)
        rows.append(row)
    return pd.DataFrame(rows)


def permutations(strata, rng):
    """Index matrix (samples x NPERM) that shuffles individuals within their stratum."""
    out = np.empty((len(strata), NPERM), int)
    for g in (np.flatnonzero(strata == s) for s in np.unique(strata)):
        for p in range(NPERM):
            out[g, p] = rng.permutation(g)
    return out


# BH q-values and classes of the edges
def with_q(table):
    """Permutation P, BH q within the regulator family of each unit and layer (alternative sources: own family), CI."""
    table = table.rename(columns={"p_s": "perm_P", "p_u": "perm_P_u"})
    kind = np.where(table.source == "m6A_load", "load", "binding")
    group = np.where(table.primary, table.family, "alt_" + table.family)
    group = pd.Series(group, index=table.index) + np.where(table.family == "control_rbp", "|" + kind, "")
    table["q"] = table.groupby(["unit", "layer", group]).perm_P.transform(utils.bh)
    table["ci_low"], table["ci_high"] = table.b_s - utils.Z95 * table.se_s, table.b_s + utils.Z95 * table.se_s
    return table


def classify(edges):
    """Edge classes. class_vs_comparator: m6A-specific when q passes and |b_s| exceeds the comparator, u (readers) or
    the negative-control regulators (writers, erasers; none in a unit without them: q only). edge_class: also |b_s|
    beyond the 95th percentile of the control RBPs in the unit and layer (reader binding edges against control binding
    edges, load edges against control load edges)."""
    control = (edges[edges.family == "negative_control"].groupby(["unit", "layer"]).b_s.agg(lambda b: b.abs().max())
               .rename("b_negative"))
    edges = edges.join(control, on=["unit", "layer"])
    edges["control"] = np.where(edges.b_u.notna(), "u", np.where(edges.b_negative.notna(), "negative controls", "none"))
    edges["b_control"] = edges.b_u.fillna(edges.b_negative)
    specific = (edges.q < utils.Q) & (edges.b_s.abs() > edges.b_control.abs().fillna(0))
    binding = (edges.perm_P_u < 0.05) & ~(edges.q < utils.Q)
    edges["class_vs_comparator"] = np.where(specific, "m6A-specific", np.where(binding, "Binding-only", "Null"))
    edges["kind"] = np.where(edges.source == "m6A_load", "load", "binding")
    limit = (edges[edges.family == "control_rbp"].groupby(["unit", "layer", "kind"])
             .b_s.agg(control_p95=lambda b: b.abs().quantile(0.95), n_controls="count").reset_index())
    edges = edges.merge(limit, on=["unit", "layer", "kind"], how="left")  # no controls in the unit: nothing is specific
    beyond_controls = edges.b_s.abs() > edges.control_p95
    edges["edge_class"] = np.where((edges.class_vs_comparator == "m6A-specific") & ~beyond_controls, "Null",
                                   edges.class_vs_comparator)
    edges.loc[edges.family.isin(["negative_control", "control_rbp"]), ["class_vs_comparator", "edge_class"]] = ""
    return edges
