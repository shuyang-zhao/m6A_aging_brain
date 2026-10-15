# Sign check of the network edges against public perturbation experiments (descriptive; not used to fit the model).
# A positive edge means the transcripts bound by a regulator with m6A follow its level, so loss of the regulator should
# lower them (gain raise them) relative to other transcripts of the same expression and length; a negative edge
# predicts the opposite. For each experiment: log2 fold change of the regulator's m6A targets (readers: mouse genes with
# m6A in the primary binding source; writers and erasers: genes with a mouse m6A site) minus that of level- and
# length-matched genes outside the target set, with a bootstrap over the matched pairs; the predicted sign comes from
# the median edge of the mouse units.
# reads edges.tsv and mouse_targets.tsv.gz of script 33, features of script 09, the perturbation tables in data/raw
# (GSE140565, GSE143659, GSE171790, GSE100491, GSE113793, GSE154992, GSE325852, GSE207145) and the Koranda 2018 table in
# ref/target_sets; writes output/34_perturbation_signs/sign_check.tsv
import re

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import utils

NETWORK = utils.OUT / "33_regulator_network"
OUT = utils.OUT / "34_perturbation_signs"
RAW = utils.DATA / "raw"
NBOOT = 2000
rng = np.random.default_rng(utils.SEED)
perturbations = []  # (id, regulator, mode, tissue, readout, table with gene, log2fc, level)


def add(name, table, regulator, mode, tissue, readout):
    table = table.dropna(subset=["gene", "log2fc", "level"]).drop_duplicates("gene")
    perturbations.append((name, regulator, mode, tissue, readout, table))


def from_counts(counts, ko, wt):
    cpm = counts[ko + wt].astype(float) / counts[ko + wt].astype(float).sum() * 1e6
    log = np.log2(cpm + 1)
    table = pd.DataFrame({"log2fc": log[ko].mean(axis=1) - log[wt].mean(axis=1),
                          "level": np.log2(cpm.mean(axis=1) + 1)})
    return table[cpm.mean(axis=1) >= 1].rename_axis("gene").reset_index()


def from_deseq(path, sign=1, sep="\t"):
    d = pd.read_csv(path, sep=sep, index_col=0)
    d.index = d.index.astype(str).str.strip('"')
    d = d[d.baseMean >= 10]
    return pd.DataFrame({"gene": d.index, "log2fc": sign * d.log2FoldChange.to_numpy(),
                         "level": np.log2(d.baseMean.to_numpy() + 1)})


# Fmr1 KO cortex: Shu 2019 (RNA, ribosome footprints), Seo 2022 (RNA, footprints, translation efficiency)
add("fmr1ko_cortex_rna_gse140565", from_deseq(RAW / "shu2019_gse140565/GSE140565_RNADEseq_WTxFK.txt.gz", sign=-1),
    "FMR1", "loss", "cortex", "RNA")
add("fmr1ko_cortex_rpf_gse140565", from_deseq(RAW / "shu2019_gse140565/GSE140565_RPFDEseq_FKxWT.txt.gz"),
    "FMR1", "loss", "cortex", "Ribo-seq")
for kind, readout in (("rna", "RNA"), ("ribo", "Ribo-seq"), ("te", "translation efficiency")):
    add(f"fmr1ko_cortex_{kind}_gse143659",
        from_deseq(RAW / f"seo2022_gse143659/GSE143659_{kind}_fxsVwt_all.csv.gz", sep=","),
        "FMR1", "loss", "cortex", readout)
# Fmr1 KO cortex, hippocampus, cerebellum (GSE207145, counts) and RiboTag translatome of Pvalb and Camk2a neurons
# (GSE325852)
d = pd.read_csv(RAW / "gse207145_fmr1_tissues/GSE207145_RNAseq.DE.log2FC_Brain.tsv.gz", sep="\t", index_col=0,
                usecols=lambda c: c == "gene_symbol" or c.startswith("num_"))
for t, tissue in (("CTX", "cortex"), ("HC", "hippocampus"), ("CB", "cerebellum")):
    ko, wt = [f"num_KO{i}_{t}" for i in (1, 2, 3)], [f"num_WT{i}_{t}" for i in (1, 2, 3)]
    add(f"fmr1ko_{tissue}_rna_gse207145", from_counts(d, ko, wt), "FMR1", "loss", tissue, "RNA")
d = pd.read_csv(RAW / "gse325852_fmr1_translatome/GSE325852_ALL_Counts_Pvalb_Camk2.csv.gz", index_col=0)
meta = pd.DataFrame([re.match(r"B2_\d+_([PC])\d+_(WT|KO)_[FM]_(V1|S1)_IP", c).groups() for c in d.columns],
                    index=d.columns, columns=["cell", "genotype", "area"])
for cell, name in (("P", "Pvalb neurons"), ("C", "Camk2a neurons")):
    columns = meta.index[meta.cell == cell]
    cpm = d[columns] / d[columns].sum() * 1e6
    log, m = np.log2(cpm + 1), meta.loc[columns]
    fc = sum(log[m.index[(m.genotype == "KO") & (m.area == a)]].mean(axis=1)
             - log[m.index[(m.genotype == "WT") & (m.area == a)]].mean(axis=1) for a in ("V1", "S1")) / 2
    table = pd.DataFrame({"gene": cpm.index, "log2fc": fc.to_numpy(),
                          "level": np.log2(cpm.mean(axis=1) + 1).to_numpy()})
    add(f"fmr1ko_{cell}_ribotag_gse325852", table[cpm.mean(axis=1).to_numpy() >= 1], "FMR1", "loss",
        f"cortex {name}", "RiboTag")
# Mettl14 cKO striatum (Koranda 2018)
k = pd.read_csv(utils.ROOT / "ref/target_sets/stage2/koranda2018_mettl14.tsv", sep="\t",
                usecols=["comparison", "gene_symbol", "log2fc_exp_over_ctrl", "base_mean_rpkm"], low_memory=False)
for comparison, tissue, tag in (("D1R-Cre_Mettl14cKO_vs_ctrl", "striatum D1 MSN", "d1msn"),
                                ("A2AR-Cre_Mettl14cKO_vs_ctrl", "striatum D2 MSN", "d2msn"),
                                ("AAV-hSyn-Cre_vs_AAV-hSyn-GFP", "striatum AAV-Cre", "aav")):
    x = k[(k.comparison == comparison) & (k.base_mean_rpkm >= 1)]
    add(f"mettl14cko_{tag}_rna_gse114675",
        pd.DataFrame({"gene": x.gene_symbol, "log2fc": x.log2fc_exp_over_ctrl, "level": np.log2(x.base_mean_rpkm + 1)}),
        "METTL14", "loss", tissue, "RNA")
# Mettl3 and Fto cKO in adult excitatory neurons (Engel 2018)
e = pd.read_csv(RAW / "gse113793_mettl3_fto_cortex/RNA_analysis.csv")
e = e[e.Diff_Mettl3_WT_baseMean >= 10]
for regulator in ("METTL3", "FTO"):
    add(f"{regulator.lower()}cko_cortex_rna_gse113793",
        pd.DataFrame({"gene": e.Gene, "log2fc": e[f"Diff_{regulator.capitalize()}_WT_log2FoldChange"],
                      "level": np.log2(e.Diff_Mettl3_WT_baseMean + 1)}),
        regulator, "loss", "cortex and hippocampus neurons", "RNA")
# Prrc2a cKO and Fto transgene (Wu 2019; one sample per group)
for name, wt, ko, regulator, mode, tissue in (
        ("prrc2acko_nestin", "GSM2684967_WT", "GSM2684968_cKO", "PRRC2A", "loss", "brain Nestin-Cre"),
        ("prrc2acko_olig2", "GSM2803005_prrc2a-WT", "GSM2803006_prrc2a-olig2-cre", "PRRC2A", "loss", "brain Olig2-Cre"),
        ("ftotg", "GSM2803007_Fto-WT", "GSM2803008_Fto-Tg", "FTO", "gain", "brain Fto transgene")):
    c = pd.DataFrame({g: pd.read_csv(RAW / f"wu2019_gse100491/{f}-reads-count.txt.gz", sep="\t", index_col=0).iloc[:, 0]
                      for g, f in (("wt", wt), ("ko", ko))})
    add(f"{name}_rna_gse100491", from_counts(c, ["ko"], ["wt"]), regulator, mode, tissue, "RNA")
# Ythdf2 cKO dentate gyrus (Zhuang 2023)
d = pd.read_csv(RAW / "huang2021_gse171790/GSE171790_cKO-RNAseq_gene_count.xls.gz", sep="\t", index_col=0)
add("ythdf2cko_dentate_rna_gse171790", from_counts(d, ["cKO1", "cKO2", "cKO3"], ["Ctrl1", "Ctrl2", "Ctrl3"]),
    "YTHDF2", "loss", "dentate gyrus", "RNA")
# embryonic cortex Mettl3 and Fto cKO, RNA and ribosome footprints (Du 2021)
for kind, name, readout in (("RNA-Seq", "rna", "RNA"), ("Ribo-seq", "ribo", "Ribo-seq")):
    d = pd.read_csv(RAW / f"gse154992_mettl3_ribo/GSE154992_{kind}_read_count.txt.gz", sep="\t", index_col=0)
    d = d[[c for c in d.columns if "Readcount" in c or c.endswith("_count")]]
    for regulator, tag in (("METTL3", "ME"), ("FTO", "FE")):
        add(f"{regulator.lower()}cko_e155_{name}_gse154992",
            from_counts(d, [c for c in d.columns if f"_{tag}" in c], [c for c in d.columns if "_WT" in c]),
            regulator, "loss", "embryonic cortex", readout)

# targets, edges and matched controls
features = utils.read(utils.OUT / "09_transcript_features/mouse_features.tsv.gz",
                      usecols=["gene_id", "symbol", "log_tx_len", "n_sites_m6A"]).dropna(subset=["log_tx_len"])
features = features.drop_duplicates("symbol")
reader_targets = utils.read(NETWORK / "mouse_targets.tsv.gz").groupby("regulator").gene_id.agg(set)
edges = utils.read(NETWORK / "edges.tsv")
edges = edges[edges.primary & (edges.species == "mouse") & (edges.family != "module")]
edges = edges.groupby(["regulator", "layer"]).agg(
    n_units=("b_s", "size"), median_b_s=("b_s", "median"), n_positive=("b_s", lambda b: (b > 0).sum()),
    n_m6A_specific=("edge_class", lambda c: (c == "m6A-specific").sum())).reset_index()


def match(table, marked):
    """One unmarked partner per marked gene, nearest on z-scored level and transcript length, without replacement."""
    m, u = table[marked], table[~marked]
    columns = ["level", "log_tx_len"]
    z = lambda x: (x - table[columns].mean()) / table[columns].std()
    _, nearest = cKDTree(z(u[columns]).to_numpy()).query(z(m[columns]).to_numpy(), k=min(30, len(u)))
    used, pairs = set(), []
    for i in rng.permutation(len(m)):
        for j in nearest[i]:
            if j not in used:
                used.add(j)
                pairs.append((i, j))
                break
    pairs = np.array(pairs)
    return m.log2fc.to_numpy()[pairs[:, 0]], u.log2fc.to_numpy()[pairs[:, 1]]


rows = []
for name, regulator, mode, tissue, readout, table in perturbations:
    by_id = table.gene.str.startswith("ENSMUSG").any()
    table = table.merge(features, left_on="gene", right_on="gene_id" if by_id else "symbol").drop_duplicates("gene_id")
    targets = reader_targets.get(regulator)
    marked = table.gene_id.isin(targets).to_numpy() if targets is not None else (table.n_sites_m6A >= 1).to_numpy()
    if marked.sum() < 50:
        continue
    target, control = match(table, marked)
    boot = (target - control)[rng.integers(0, len(target), (NBOOT, len(target)))].mean(axis=1)
    effect = (target - control).mean()
    p = min(1, 2 * (min((boot <= 0).sum(), (boot >= 0).sum()) + 1) / (NBOOT + 1))
    for edge in edges[edges.regulator == regulator].itertuples():
        predicted = -np.sign(edge.median_b_s) if mode == "loss" else np.sign(edge.median_b_s)
        rows.append(dict(perturbation=name, regulator=regulator, mode=mode, tissue=tissue, readout=readout,
                         n_targets=len(target), effect=effect, ci_low=np.quantile(boot, 0.025),
                         ci_high=np.quantile(boot, 0.975), boot_P=p, layer=edge.layer, n_units=edge.n_units,
                         median_b_s=edge.median_b_s, n_positive=edge.n_positive, n_m6A_specific=edge.n_m6A_specific,
                         predicted_sign=predicted, observed_sign=np.sign(effect),
                         agrees=float(predicted == np.sign(effect))))
result = pd.DataFrame(rows)
result["q"] = result.groupby("layer").boot_P.transform(utils.bh)
utils.write(result, OUT / "sign_check.tsv")
print(result.groupby(["regulator", "layer"]).agg(n=("agrees", "size"), agree=("agrees", "sum"),
                                                 significant=("q", lambda q: (q < utils.Q).sum())).to_string())
print(f"{int(result.agrees.sum())} of {len(result)} comparisons agree; "
      f"among q < 0.05: {int(result[result.q < utils.Q].agrees.sum())} of {(result.q < utils.Q).sum()}")
