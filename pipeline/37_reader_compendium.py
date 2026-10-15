# Compendium of IGF2BP2 and YTHDF1-3: open binding and perturbation datasets, harmonised for the bridge tests of
# script 38.
# Binding: a source is a peak set (CLIP; replicate libraries combined by overlap, lifted to GRCh38/mm10 where needed) or
# a gene list (RIP-seq, or a published list of CLIP targets). A gene is bound by a source when a peak overlaps its
# exons, and carries m6A when a site of the primary m6A map (mouse: the single-base maps; human: the brain m6A-seq
# peaks) lies within 50 nt of a peak (gene lists: when the gene has any site). IGF2BP2 tiers: T1 bound with m6A, T2
# the T1 genes whose m6A peak falls after METTL14 knockdown (HepG2 m6A-seq of Huang 2018), T3 the T1 genes that fall
# after IGF2BP2 loss. YTHDF target sets (read by script 38): the preferred source of each paralog and species (in vivo
# brain CLIP, else the consensus of at least two neural sources or the single neural source, else the consensus of all
# its sources, else the single source), the consensus of all sources and every source alone, each as T1 and
# bound-unmarked genes.
# Perturbation: log2 fold change of every open knockdown, knockout, inhibitor or overexpression dataset, with the BH
# padj where counts or deposited statistics give one, by readout (RNA, translation efficiency from Ribo-seq, protein).
# reads data/raw (script 00; the ENCODE IDR peaks, clipflow and clip_standard peaks come with the data release), the
# site maps (04), the CLIP peaks (07), gene features (09) and ref/target_sets; writes output/37_reader_compendium/
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from pyliftover import LiftOver

import network as nw
import utils

OUT = utils.OUT / "37_reader_compendium"
RAW = utils.DATA / "raw"
HERE = Path(__file__).resolve().parent
NAMES = ["chrom", "start", "end", "strand"]
CHROMS = {f"chr{c}" for c in [*range(1, 23), "X", "Y"]}
chains = {"hg19": "ucsc_liftover/hg19ToHg38.over.chain.gz", "mm9": "preitner2016_gse83822/mm9ToMm10.over.chain.gz"}
chains = {build: LiftOver(str(RAW / path)) for build, path in chains.items()}
human, mouse = nw.symbols["human"], nw.symbols["mouse"]


def gene_id(x, species):
    """Ensembl ids (version stripped) kept when known, other entries mapped from their symbol."""
    x = pd.Series(x).astype(str)
    ensembl = x.str.split(".").str[0]
    return np.where(ensembl.isin(nw.features[species].index), ensembl, x.str.upper().map(nw.symbols[species]))


def de(mode, matrix, ref, test, rpf=None):
    """log2 fold change and BH padj of test against ref with perturbation_de.R; te modes: matrix is the RNA, rpf the
    Ribo-seq."""
    groups = pd.DataFrame({"sample": ref + test, "group": ["ref"] * len(ref) + ["test"] * len(test)})
    matrix = matrix[ref + test]
    if rpf is not None:
        matrix = pd.concat([matrix.add_prefix("rna_"), rpf[ref + test].add_prefix("rpf_")], axis=1).dropna()
        groups = pd.DataFrame({"sample": matrix.columns, "assay": [c[:3] for c in matrix.columns],
                               "group": ["ref" if c[4:] in ref else "test" for c in matrix.columns]})
    with tempfile.TemporaryDirectory() as folder:
        matrix.to_csv(f"{folder}/matrix.tsv", sep="\t")
        groups.to_csv(f"{folder}/groups.tsv", sep="\t", index=False)
        run = subprocess.run(["Rscript", str(HERE / "perturbation_de.R"), mode, f"{folder}/matrix.tsv",
                              f"{folder}/groups.tsv", f"{folder}/out.tsv"], capture_output=True, text=True)
        if run.returncode:
            raise RuntimeError(run.stderr)
        return pd.read_csv(f"{folder}/out.tsv", sep="\t").rename(columns={"gene": "gene_id"})


# perturbation datasets
rows = []


def add(dataset, reader, species, cell, direction, readout, table, accession, note):
    table = (table.assign(gene_id=gene_id(table.gene_id, species)).dropna(subset=["gene_id", "log2fc"])
             .drop_duplicates("gene_id"))
    rows.append(table[np.isfinite(table.log2fc)].assign(
        dataset=dataset, reader=reader, accession=accession, species=species, cell_type=cell, direction=direction,
        readout=readout, note=note)[["dataset", "reader", "accession", "species", "cell_type", "direction", "readout",
                                     "gene_id", "log2fc", "padj", "note"]])


def rna(dataset, reader, species, cell, direction, matrix, ref, test, accession, note, mode="edger"):
    add(dataset, reader, species, cell, direction, "RNA", de(mode, matrix, ref, test), accession, note)


def deposited(dataset, reader, species, cell, direction, readout, gene, log2fc, padj, accession, note):
    add(dataset, reader, species, cell, direction, readout,
        pd.DataFrame({"gene_id": gene, "log2fc": log2fc, "padj": padj}), accession, note)


def sample_counts(folder, files, **kw):
    return pd.concat({name: pd.read_csv(RAW / folder / f, sep="\t", header=None, index_col=0, **kw)[1]
                      for name, f in files.items()}, axis=1)


def ribo_rpkm(name):
    """Ribo-seq RPKM of a Lasman 2020 library, summed over the isoforms of each gene."""
    d = pd.read_csv(next((RAW / "lasman2020_gse148039").glob(f"*_{name}.ribo.txt.gz")), sep="\t", usecols=[0, 2],
                    header=0, names=["gene", "rpkm"])
    return pd.to_numeric(d.rpkm, errors="coerce").groupby(d.gene).sum()


def perturbations():
    # IGF2BP2: ENCODE shRNA RNA-seq (DESeq2 of the ENCODE pipeline, 2 vs 2)
    for name, cell, path, accession in [
            ("encode_hepg2_shRNA", "HepG2", "encode_regulator_rnaseq/ENCSR478FJK/ENCFF302MAE.tsv", "ENCSR478FJK"),
            ("encode_k562_shRNA", "K562", "encode_igf2bp2_kd_k562/ENCFF407VNW.tsv", "ENCSR952RRH")]:
        d = pd.read_csv(RAW / path, sep="\t")
        deposited(name, "IGF2BP2", "human", cell, "KD", "RNA", d.index, d.log2FoldChange, d.padj, accession,
                  "ENCODE DESeq2, 2 vs 2")
    # Degrauwe 2016 (GSE73845): glioma stem-like cells IMP2s knockdown and IMP2a overexpression, Imp2 knockout mouse
    # brain
    for tag, direction, cell, species in [("IMP2s_KD", "KD", "MGG4+MGG8 glioma stem cells", "human"),
                                          ("IMP2a_OE", "OE", "MGG4+MGG8 glioma stem cells", "human"),
                                          ("brain_IMP2_KO", "KO", "mouse brain (E15 + P21)", "mouse")]:
        d = pd.read_csv(RAW / f"degrauwe2016_gse73845/GSE73845_RNA_{tag}.txt.gz", sep="\t")
        lfc = next(c for c in d.columns if c.endswith("_LFC"))
        deposited(f"degrauwe2016_{tag}", "IGF2BP2", species, cell, direction, "RNA", d.symbol, d[lfc],
                  d[lfc.replace("LFC", "FDR")], "GSE73845", "deposited")
    # Huang 2018 (GSE90684): HepG2 shIGF2BP2, average RPKM only (no replicates, no padj)
    d = pd.read_csv(RAW / "huang2018_gse90684/GSE90684_expression_matrix_of_genes_rpkm_average_result.txt.gz", sep="\t")
    d = d[d.shNs_exp >= 1]
    deposited("huang2018_hepg2_shIGF2BP2", "IGF2BP2", "human", "HepG2", "KD", "RNA", d.GeneName,
              np.log2((d.shIGF2BP2_exp + 0.1) / (d.shNs_exp + 0.1)), np.nan, "GSE90684",
              "average RPKM, RPKM >= 1 in control; no padj")
    # GSE146726: bladder cancer lines, siIGF2BP2 against scramble, 3 vs 3 deposited normalised values (limma-trend on
    # log2(x + 1))
    for cell in ("J82", "UMUC3"):
        d = pd.read_csv(RAW / f"gse146726_igf2bp2_kd/GSE146726_{cell}_DESeq2_IGF2BP2.txt.gz", sep="\t")
        d.index = d.index.astype(str).str.split(".").str[0]
        rna(f"gse146726_{cell}_siIGF2BP2", "IGF2BP2", "human", cell, "KD", np.log2(d + 1),
            [f"scramble_{i}" for i in "ABC"], [f"IGF2BP2_{i}" for i in "ABC"], "GSE146726",
            "limma-trend on log2(x + 1) of the deposited normalised values", "limma")
    # GSE292558: bone-marrow macrophages, Igf2bp2 knockout (deposited DESeq2); GSE119801: adult muscle Imp2 knockout, 3
    # vs 3 counts
    d = pd.read_csv(RAW / "gse292558_igf2bp2_ko_bmm/GSE292558_DESeq2_IGF2BP2_KO_vs_WT.txt.gz", sep="\t")
    deposited("gse292558_bmm_ko", "IGF2BP2", "mouse", "bone-marrow macrophages", "KO", "RNA", d.ensembl_gene_id,
              d.log2FoldChange, d.padj, "GSE292558", "Igf2bp2 mRNA itself not reduced")
    wt = {f"WT{i}": f"GSM33840{27 + i}_WT.total.{i}.count.txt.gz" for i in (1, 2, 3)}
    ko = {f"KO{i}": f"GSM33840{30 + i}_KO.Total.{i}.count.txt.gz" for i in (1, 2, 3)}
    rna("regue2019_muscle_ko", "IGF2BP2", "mouse", "adult skeletal muscle", "KO",
        sample_counts("regue2019_gse119801", wt | ko), list(wt), list(ko), "GSE119801", "edgeR QL, 3 vs 3")

    # YTHDF, neural first. Yu 2021 (GSE153688): cerebellar granule cells shYthdf2 RNA-seq; shYthdf1 proteome (regulated
    # proteins of Table S3 only)
    m = pd.read_csv(RAW / "yu2021_gse153688/GSE153688_Raw_gene_counts_matrix_KD.txt.gz", sep="\t", index_col=0)
    rna("yu2021_cgc_shYthdf2", "YTHDF2", "mouse", "cerebellar granule cells (P8)", "KD", m,
        ["KDC_1", "KDC_2", "KDC_4"], ["KD2_3_1", "KD2_3_2", "KD2_3_4"], "GSE153688", "edgeR QL, 3 vs 3")
    p = pd.read_excel(RAW / "yu2021_advsci_supp_pmc8596126/ADVS-8-2101329-s002.xlsx")
    deposited("yu2021_cgc_shYthdf1_proteome", "YTHDF1", "mouse", "cerebellar granule cells (P8)", "KD", "protein",
              p["Gene name"], np.log2(p["shYthdf1/shCtrlRatio"]), np.nan, "Yu 2021 Table S3",
              "regulated proteins only (P < 0.05); no full table, no padj")
    d = pd.read_csv(RAW / "huang2021_gse171790/GSE171790_cKO-RNAseq_gene_count.xls.gz", sep="\t", index_col=0)
    rna("zhuang2023_dg_Ythdf2_cKO", "YTHDF2", "mouse", "hippocampal dentate gyrus granule cells (P40)", "cKO", d,
        ["Ctrl1", "Ctrl2", "Ctrl3"], ["cKO1", "cKO2", "cKO3"], "GSE171790", "edgeR QL, 3 vs 3")
    d = pd.read_excel(RAW / "gse306850_ythdf1/GSE306850_gene_count_matrix.xlsx").set_index("gene_id")
    rna("gse306850_acc_Ythdf1_KO", "YTHDF1", "mouse", "anterior cingulate cortex", "KO", d,
        [f"WT_ACC_{i}" for i in (1, 2, 3)], [f"KO_ACC_{i}" for i in (1, 2, 3)], "GSE306850", "edgeR QL, 3 vs 3")
    d = pd.read_excel(RAW / "shi2018_gse106607/GSE106607_C_Hippocampal_RNA_seq_WT_KO_processed.xlsx").set_index("rpkm")
    d.index = gene_id(d.index, "mouse")
    d = d[pd.notna(d.index) & ~d.index.duplicated()]
    rna("shi2018_hippocampus_Ythdf1_KO", "YTHDF1", "mouse", "adult hippocampus", "KO", np.log2(d + 1),
        [f"WT_{i}" for i in (1, 2, 3)], [f"KO_{i}" for i in (1, 2, 3)], "GSE106607",
        "limma-trend on log2(RPKM + 1) of the deposited RPKM, 3 vs 3", "limma")
    # Li 2018 (GSE104867): Ythdf2 cKO neural progenitors, input RNA-seq
    ids = {"WT1": "00", "WT2": "01", "WT3": "02", "KO1": "03", "KO2": "04", "KO3": "05"}
    m = sample_counts("li2018_gse104867", {k: f"GSM28091{v}_input_{k}.raw.counts.txt.gz" for k, v in ids.items()},
                      skiprows=1)
    rna("li2018_nspc_Ythdf2_cKO", "YTHDF2", "mouse", "neural stem/progenitor cells (E14.5)", "cKO",
        m[~m.index.str.startswith("__")], ["WT1", "WT2", "WT3"], ["KO1", "KO2", "KO3"], "GSE104867",
        "edgeR QL, 3 vs 3")
    m = (pd.read_csv(RAW / "gse266096_ythdf2_nsc/GSE266096_Raw_count_table_RNA-seq.txt.gz", sep=r"\s+", index_col=0)
         .groupby(level=0).sum())
    rna("gse266096_nsc_Ythdf2_cKO", "YTHDF2", "mouse", "hippocampal neural stem cells (P0.5)", "cKO", m,
        [f"Decay0h_WT{i}" for i in (1, 2, 3)], [f"Decay0h_cKO{i}" for i in (1, 2, 3)], "GSE266096",
        "0 h of the actinomycin D series; Ythdf2 mRNA is higher in the cKO (exon-deleted transcript)")
    # Zou 2023 (GSE214882): primary neurons, YTHDF1 inhibitor SAC against DMSO, Ribo-seq and input RNA, 2 vs 2, wild
    # type and Fmr1 knockout
    x = (pd.read_excel(RAW / "zou2023_gse214882/GSE214882_posterior_mean_estimate.xlsx", sheet_name="Neuron_RSEM_PMC")
         .set_index("symbol").drop(columns="ens_id"))
    x = x[~x.index.duplicated()].fillna(0).round()
    for genotype, label in (("wt", "wild-type"), ("ko", "Fmr1-KO")):
        ref, test = [f"{genotype}_dmso_r{i}" for i in (1, 2)], [f"{genotype}_sac_r{i}" for i in (1, 2)]
        counts, rpf = (pd.DataFrame({k: x[f"{assay}_{k}"] for k in ref + test}) for assay in ("input", "rpf"))
        for readout, table in (("TE", de("te", counts, ref, test, rpf)), ("RNA", de("edger", counts, ref, test))):
            add(f"zou2023_neuron_{genotype}_YTHDF1inh", "YTHDF1", "mouse", f"primary neurons ({label})", "inhibitor",
                readout, table, "GSE214882",
                "SAC against DMSO; edgeR QL, 2 vs 2; "
                + ("assay x treatment interaction" if readout == "TE" else "input RNA"))

    # Ribo-seq in cell lines. Lasman 2020 (GSE148536, GSE147849): mESC single and triple Ythdf knockouts; the deposited
    # Ribo-seq values are RPKM, so TE is the assay x genotype interaction of log2 values (RNA: DESeq-normalised counts;
    # Ribo-seq: RPKM summed over isoforms)
    counts = pd.read_csv(RAW / "lasman2020_gse147849/GSE147849_DESeq_mESCs.txt.gz", sep="\t", index_col=0)
    counts = (counts[[c for c in counts.columns if c.endswith(".normalized")]]
              .rename(columns=lambda c: c.replace(".normalized", "")))
    library = {"WT_A": "WT_A_S1", "WT_B": "WT_B_S2", "Ythdf1_KO_A": "Ythdf1_KO_A_S3", "Ythdf1_KO_B": "Ythdf1_KO_B_S4",
               "Ythdf2_KO_A": "Ythdf2_KO_A_S5", "Ythdf2_KO_B": "Ythdf2_KO_B_S6", "Ythdf3_KO_A": "Ythdf3_KO_A_S7",
               "Ythdf3_KO_B": "Ythdf3_KO_B_S8", "Ythdf321_KO_A": "Ythdf3_2_1_KO_A_S9",
               "Ythdf321_KO_B": "Ythdf3_2_1_KO_B_S10"}
    rpf = pd.concat({k: ribo_rpkm(k) for k in library}, axis=1)
    counts = counts[list(library.values())].set_axis(list(library), axis=1)
    common = counts.index.intersection(rpf.index)
    counts, rpf = counts.loc[common], rpf.loc[common]
    for genotype, reader, direction in (("Ythdf1_KO", "YTHDF1", "KO"), ("Ythdf2_KO", "YTHDF2", "KO"),
                                        ("Ythdf3_KO", "YTHDF3", "KO"), ("Ythdf321_KO", "YTHDF1/2/3", "triple KO")):
        ref, test = ["WT_A", "WT_B"], [f"{genotype}_A", f"{genotype}_B"]
        keep = (counts[ref + test].mean(axis=1) >= 10) & (rpf[ref + test].mean(axis=1) >= 1)
        for readout, table in (("TE", de("te_limma", np.log2(counts[keep] + 1), ref, test, np.log2(rpf[keep] + 1))),
                               ("RNA", de("limma", np.log2(counts + 1), ref, test))):
            add(f"lasman2020_mESC_{genotype}", reader, "mouse", "mouse embryonic stem cells", direction, readout, table,
                "GSE148536/GSE147849",
                "limma-trend, 2 vs 2" + ("; assay x genotype interaction on log2(RPKM + 1) and "
                                         "log2(normalised counts + 1)" if readout == "TE" else ""))
    # HeLa Ribo-seq with deposited log2 ratios (no padj): Wang 2015 siYTHDF1 (mean of 2 replicates), Shi 2017 siYTHDF3
    # (2 replicates), Wang 2014 siYTHDF2 (source data)
    d = pd.read_excel(RAW / "wang2015_gse63591/GSE63591_C-Y1-Ribosome_profiling.xlsx", header=None, skiprows=2)
    for readout, column in (("TE", 17), ("RNA", 16)):
        deposited("wang2015_HeLa_siYTHDF1", "YTHDF1", "human", "HeLa", "KD", readout, d[0],
                  pd.to_numeric(d[column], errors="coerce"), np.nan, "GSE63591",
                  "deposited mean of 2 replicates, no padj")
    d = pd.read_excel(RAW / "shi2017_gse86214/GSE86214_siYTHDF3_RiboProfile.xlsx", header=None, skiprows=2)
    for readout, columns in (("TE", (14, 15)), ("RNA", (9, 10))):
        deposited("shi2017_HeLa_siYTHDF3", "YTHDF3", "human", "HeLa", "KD", readout, d[0],
                  d[list(columns)].apply(pd.to_numeric, errors="coerce").mean(axis=1, skipna=False), np.nan,
                  "GSE86214", "mean of the 2 replicate log2 ratios, no padj")
    d = pd.read_excel(utils.ROOT / "ref/target_sets/raw/wang2014/41586_2014_BFnature12730_MOESM262_ESM.xlsx",
                      header=None, skiprows=2)
    for readout, column in (("TE", 9), ("RNA", 8)):
        deposited("wang2014_HeLa_siYTHDF2", "YTHDF2", "human", "HeLa", "KD", readout, d[1],
                  pd.to_numeric(d[column], errors="coerce"), np.nan, "GSE49339",
                  "deposited log2(siYTHDF2/siControl), no padj")

    # RNA only. Zou 2023 (GSE201070): HeLa siYTHDF1-3, deposited log2FC (YTHDF3 from the replicate means), no padj
    d = pd.read_csv(RAW / "zou2023_gse201070/GSE201070_HeLa_DF_Knockdown.csv.gz")
    for k, reader in (("DF1", "YTHDF1"), ("DF2", "YTHDF2")):
        deposited(f"zou2023_HeLa_si{k}", reader, "human", "HeLa", "KD", "RNA", d.ensembl_gene_id, d[f"log2FC_{k}KD"],
                  np.nan, "GSE201070", "deposited log2FC, no padj")
    ratio = np.log2((d[["siDF3_rep1", "siDF3_rep2"]].mean(axis=1) + 1e-4)
                    / (d[["siCtrl_rep1", "siCtrl_rep2"]].mean(axis=1) + 1e-4))
    deposited("zou2023_HeLa_siDF3", "YTHDF3", "human", "HeLa", "KD", "RNA", d.ensembl_gene_id, ratio, np.nan,
              "GSE201070", "log2 of the replicate-mean deposited values (pseudocount 1e-4), no padj")
    # Heck 2020: iPSC YTHDF2 knockdown, 3 vs 4
    m = pd.read_csv(RAW / "heck2020_gse133898/GSE133898_Control_YTH2KD_HTseq_raw_counts.txt.gz", sep="\t", index_col=0)
    rna("heck2020_iPSC_siYTHDF2", "YTHDF2", "human", "iPSC (neural-differentiation model)", "KD", m,
        [f"Control_{i}" for i in (1, 2, 3)], [f"YTH2KD_{i}" for i in (1, 2, 3, 4)], "GSE133898", "edgeR QL, 3 vs 4")
    # SNU638 gastric cancer cells, YTHDF1 overexpression (both arms IFN-gamma treated)
    m = pd.read_csv(RAW / "gse270417_ythdf_gastric/GSE270417_raw_counts.txt.gz", sep="\t")
    m = (m.groupby("Gene_Symbol")[[c for c in m.columns if c.endswith("Read_Count")]].sum()
         .rename(columns=lambda c: c.replace("_Read_Count", "")))
    rna("gse270417_SNU638_YTHDF1_OE", "YTHDF1", "human", "SNU638 gastric cancer cells", "OE", m,
        ["Control_1", "Control_2"], ["YTHDF1_1", "YTHDF1_2"], "GSE270417",
        "edgeR QL, 2 vs 2; both arms IFN-gamma treated")


perturbations()
perturbation = pd.concat(rows, ignore_index=True)
utils.write(perturbation, OUT / "perturbation.tsv")
print(perturbation.groupby(["dataset", "reader", "readout"]).size().to_string(), flush=True)


# binding sources
def lift(peaks, build):
    """Lift 0-based half-open intervals; unmapped or split intervals are dropped."""
    out = []
    for chrom, start, end, strand in zip(peaks.chrom, peaks.start, peaks.end, peaks.strand, strict=True):
        a, b = (chains[build].convert_coordinate(chrom, int(start)),
                chains[build].convert_coordinate(chrom, int(end) - 1))
        if a and b and a[0][0] == b[0][0]:
            low, high = sorted([a[0][1], b[0][1]])
            out.append((a[0][0], low, high + 1, strand))
    return pd.DataFrame(out, columns=NAMES)


def combine(libraries):
    """Peaks of a dataset: one library as it is; several: the peaks of a library that overlap a peak of another library
    (same strand)."""
    libraries = [x[x.chrom.isin(CHROMS)].drop_duplicates(NAMES).reset_index(drop=True) for x in libraries]
    if len(libraries) == 1:
        return libraries[0]
    kept = []
    for i, a in enumerate(libraries):
        others, hit = pd.concat(libraries[:i] + libraries[i + 1:]), np.zeros(len(a), bool)
        for key, rows in a.groupby(["chrom", "strand"]):
            o = others[(others.chrom == key[0]) & (others.strand == key[1])].sort_values("start")
            if o.empty:
                continue
            reach = np.maximum.accumulate(o.end.to_numpy())
            last = np.searchsorted(o.start.to_numpy(), rows.end.to_numpy(), side="left") - 1
            hit[rows.index.to_numpy()] = (last >= 0) & (reach[np.maximum(last, 0)] > rows.start.to_numpy())
        kept.append(a[hit])
    return pd.concat(kept).drop_duplicates(NAMES).reset_index(drop=True)


def bed(path, columns, **kw):
    return pd.read_csv(path, sep="\t", header=None, usecols=list(columns), names=list(columns.values()), **kw)


def library_peaks(folder, reader, species, libraries):
    """Binding regions of the IP libraries (clip_standard: one file per library; clipflow_peaks: one per reader)."""
    if folder == "clip_standard":
        return [bed(RAW / f"clip_standard/{species}_{reader}_{lib}.bed.gz",
                    {0: "chrom", 1: "start", 2: "end", 5: "strand"}) for lib in libraries]
    b = bed(RAW / f"clipflow_peaks/{species}_{reader}.bed.gz",
            {0: "chrom", 1: "start", 2: "end", 3: "lib", 5: "strand"})
    return [b[b.lib == lib][NAMES] for lib in libraries]


def rip(counts, ip, control):
    """Genes enriched in the IP over the control (edgeR QL): log2FC > 1 and FDR < 0.05."""
    table = de("edger", counts, control, ip)
    return table[(table.log2fc > 1) & (table.padj < 0.05)].gene_id.unique()


sources = []


def source(reader, name, species, cell, tissue_class, rank, method, accession, data, note):
    """data: peaks (DataFrame) or the ids of the bound genes. rank: 1 in vivo brain CLIP, 2 brain or neural CLIP, 3
    brain or neural RIP, 4 cell-line source."""
    sources.append(dict(reader=reader, source=name, species=species, cell=cell, tissue_class=tissue_class, rank=rank,
                        method=method, accession=accession, level="peak" if isinstance(data, pd.DataFrame) else "gene",
                        data=data, note=note))


def mouse_ids(x, symbols=False):
    x = pd.Series(x).astype(str)
    return x.str.upper().map(mouse).dropna().to_numpy() if symbols else x[x.isin(nw.features["mouse"].index)].to_numpy()


def binding_sources():
    # IGF2BP2, human: Conway 2016 H9 eCLIP (input-normalised peaks, p <= 0.001 and fold >= 8 in both replicates)
    reps = []
    for n, rep in ((43, "Rep1"), (46, "Rep2")):
        d = bed(RAW / f"conway2016_gse78507/GSM20717{n}_IGF2BP2_H9ES_{rep}_eCLIP.InputNormalizedPeaks.bed.gz",
                {0: "chrom", 1: "start", 2: "end", 3: "logp", 4: "log2fc", 5: "strand"})
        reps.append(lift(d[(d.logp >= 3) & (d.log2fc >= 3)][NAMES], "hg19"))
    source("IGF2BP2", "conway2016_H9_eCLIP", "human", "H9 hESC", "neural-adjacent (pluripotent)", np.nan, "eCLIP",
           "GSE78507", combine(reps), "peaks in both replicates, hg19 lifted to hg38")
    parts = []  # Degrauwe 2016: MGG8 PAR-CLIP of the IMP2a and IMP2s isoforms, clusters in both replicates
    for k in ("IMP2a", "IMP2s"):
        parts.append(combine([lift(bed(RAW / f"degrauwe2016_gse73845/GSE73845_PARCLIP_{k}_r{i}.bed.gz",
                                       {0: "chrom", 1: "start", 2: "end", 5: "strand"}, skiprows=1), "hg19")
                              for i in (1, 2)]))
    source("IGF2BP2", "degrauwe2016_MGG8_PARCLIP", "human", "MGG8 glioma stem cells", "neural (glioma)", np.nan,
           "PAR-CLIP", "GSE73845", pd.concat(parts).drop_duplicates(NAMES).reset_index(drop=True),
           "clusters in both replicates, hg19 lifted to hg38")
    h = [bed(RAW / f"han2022_gse192792/GSM57654{n}_IMP2_IP_{i}.peaks.bed.annot.txt.gz",
             {0: "chrom", 1: "start", 2: "end", 3: "logp", 5: "strand"}, skiprows=1) for n, i in ((15, 1), (16, 2))]
    source("IGF2BP2", "han2022_A172_eCLIP", "human", "A172 glioma", "neural (glioma)", np.nan, "eCLIP", "GSE192792",
           combine([x[x.logp >= 3][NAMES] for x in h]), "-log10 P >= 3 in both replicates")
    source("IGF2BP2", "encode2020_K562_eCLIP_IDR", "human", "K562", "cell line", np.nan, "eCLIP IDR", "ENCSR062NNB",
           combine([nw.idr_peaks("IGF2BP2")]), "ENCODE IDR peaks")
    source("IGF2BP2", "encode2020_K562_eCLIP_PureCLIP", "human", "K562", "cell line", np.nan, "eCLIP PureCLIP",
           "ENCSR062NNB", combine([nw.peaks[(nw.peaks.reader == "IGF2BP2") & (nw.peaks.species == "human")][NAMES]]),
           "PureCLIP sites of script 07")
    d = pd.read_excel(RAW / "huang2018_gse90639/GSE90639_IGF2BPs_RIP_bingding_gene_and_exp_info.xlsx",
                      sheet_name="IGF2BP2_bingding_gene")
    source("IGF2BP2", "huang2018_HEK293T_RIPseq", "human", "HEK293T", "cell line", np.nan, "RIP-seq", "GSE90639",
           d.geneName.astype(str).str.upper().map(human).dropna().unique(), "deposited binding-gene list")
    # IGF2BP2, mouse: Preitner 2016 E14 brain CLIP clusters (mm9), MEF RIP-seq (GSE133242, IP against IgG, edgeR QL),
    # adipocyte RIP-seq list (GSE192705)
    d = (pd.read_csv(RAW / "preitner2016_gse83822/GSE83822_IMP2_peaks_RefSeq_mRNA.txt.gz", sep="\t").iloc[:, :4]
         .set_axis(NAMES, axis=1))
    source("IGF2BP2", "preitner2016_E14brain_CLIP", "mouse", "E14 embryonic brain", "brain", np.nan, "CLIP", "GSE83822",
           combine([lift(d.drop_duplicates(), "mm9")]), "deposited RefSeq-mRNA clusters, mm9 lifted to mm10")
    samples = {**{f"IgG{i}": f"GSM{3903890 + i}_IP.IgG.{i}.txt.gz" for i in (1, 2, 3)},
               **{f"IP{i}": f"GSM{3903896 + i}_IP.IMP2.{i}.txt.gz" for i in (1, 2, 3)}}
    source("IGF2BP2", "gse133242_MEF_RIPseq", "mouse", "mouse embryonic fibroblasts", "fibroblast", np.nan, "RIP-seq",
           "GSE133242", mouse_ids(rip(sample_counts("gse133242_imp2_mef_rip", samples), [f"IP{i}" for i in (1, 2, 3)],
                                      [f"IgG{i}" for i in (1, 2, 3)])), "IP against IgG, 3 vs 3")
    d = pd.read_csv(RAW / "gse192705_igf2bp2_rip_adipocyte/GSE192705_DESeq.input_noAb_IgG.Igfbp2.over10_either_over2"
                          "-fold_enrichment.underp05_Adjpval.txt.gz", sep="\t")
    source("IGF2BP2", "gse192705_adipocyte_RIPseq", "mouse", "cultured adipocytes", "adipocyte", np.nan, "RIP-seq",
           "GSE192705", mouse_ids(d[(d["log2(Igfbp2/CTRL)"] > 1) & (d.padj < 0.05)].id),
           "deposited list, >2-fold, padj < 0.05")

    # YTHDF1, mouse: Shi 2026 striatum CLIP (rank 1, repo list), Shi 2018 hippocampus CLIP, RIP-seq of neural tissues
    # and cells (Niu 2022, Flamand 2022, Yu 2021)
    t = (pd.read_csv(utils.ROOT / "ref/target_sets/stage2/shi2026_ythdf1_clip.tsv", sep="\t")
         .query("arm == 'saline'"))
    source("YTHDF1", "shi2026_striatum_CLIP", "mouse", "adult striatum (saline arm)", "brain in vivo", 1,
           "CLIP-seq (CLIPper)", "GSE326611", mouse_ids(t.gene_id), "genes with a peak in each of 3 replicates")
    d = pd.read_excel(RAW / "shi2018_gse106607/GSE106607_A_Ythdf1_CLIP_seq_peak_list.xlsx",
                      sheet_name="High-confidence Y1-CLIP peaks")
    source("YTHDF1", "shi2018_hippocampus_CLIP", "mouse", "adult hippocampus", "brain in vivo", 2,
           "CLIP-seq (mutation-based)", "GSE106607",
           combine([pd.DataFrame({"chrom": "chr" + d.chr.astype(str), "start": d.start - 1, "end": d.end,
                                  "strand": d.Strand})]), "deposited high-confidence peaks")
    niu = pd.read_csv(RAW / "niu2022_gse199962/GSE199962_gene_count.txt.gz", sep="\t", index_col=0)
    # the cortex IP is close to its input: no enriched gene
    for tissue, label in (("Cortex", "embryonic cortex (E12.5)"), ("Retina", "neonatal retina (P0)")):
        for paralog in ("Ythdf1", "Ythdf2"):
            source(paralog.upper(), f"niu2022_{tissue.lower()}_{paralog}_RIP", "mouse", label,
                   "neural primary tissue", 3, "RIP-seq", "GSE199962",
                   mouse_ids(rip(niu, [f"{tissue}_{paralog}_Rep{i}" for i in (1, 2)],
                                 [f"{tissue}_Input_Rep{i}" for i in (1, 2)])), "IP against input, 2 vs 2")
    for i in (1, 2, 3):
        d = pd.read_csv(RAW / f"flamand2022_gse194207/GSE194207_DF{i}_RIP_DESeq.txt.gz", sep="\t")
        source(f"YTHDF{i}", f"flamand2022_hippocampal_neurons_YTHDF{i}_RIP", "mouse", "cultured hippocampal neurons",
               "neural primary cells", 3, "RIP-seq", "GSE194207",
               mouse_ids(d[(d.log2FoldChange > 1) & (d.padj < 0.05)].ensembl_gene_id),
               "IP against input, 3 vs 3, deposited DESeq2")
    for sheet, reader in (("s004", "YTHDF1"), ("s008", "YTHDF2")):
        d = pd.read_excel(RAW / f"yu2021_advsci_supp_pmc8596126/ADVS-8-2101329-{sheet}.xlsx")
        source(reader, f"yu2021_cerebellar_granule_{reader}_RIP", "mouse", "cerebellar granule cells (P8)",
               "neural primary tissue", 3, "RIP-seq", "GSE153689", mouse_ids(d[d.iloc[:, 4] > 1].Gene_ID),
               "deposited RIP targets, IP/input log2FC > 1 and P < 0.05")
    m = (pd.read_csv(RAW / "gse266096_ythdf2_nsc/GSE266096_Raw_count_table_YTHDF2-RIP-seq.txt.gz", sep=r"\s+",
                     index_col=0).groupby(level=0).sum())
    source("YTHDF2", "gse266096_hippocampal_NSC_YTHDF2_RIP", "mouse", "hippocampal neural stem cells (P0.5)",
           "neural primary cells", 3, "RIP-seq", "GSE266096",
           mouse_ids(rip(m, [f"WT{i}_YTHDF2" for i in (1, 2, 3)], [f"WT{i}_IgG" for i in (1, 2, 3)]), symbols=True),
           "IP against IgG, 3 vs 3")
    # YTHDF2, human neural: Dixit 2020 iCLIP of glioblastoma stem cells (lines 387, 4121) and neural progenitors (HNP1),
    # clipflow truncation peaks
    for lab, libraries, cell, kind in (("387", ["GSM4785369", "GSM4785370"], "glioblastoma stem cells (line 387)",
                                        "neural cell line (glioma)"),
                                       ("4121", ["GSM4785372", "GSM4785373"], "glioblastoma stem cells (line 4121)",
                                        "neural cell line (glioma)"),
                                       ("HNP1", ["GSM4785375", "GSM4785376"], "human neural progenitors (HNP1)",
                                        "neural cell line")):
        source("YTHDF2", f"dixit2020_{lab}_iCLIP", "human", cell, kind, 2, "iCLIP", "GSE158742",
               combine(library_peaks("clipflow", "YTHDF2", "human", libraries)), "peaks in both libraries")
    # cell lines
    source("YTHDF1", "lasman2020_mESC_eCLIP", "mouse", "mouse ESC", "cell line", 4, "eCLIP (PureCLIP)", "GSE148039",
           combine(library_peaks("clip_standard", "YTHDF1", "mouse", ["GSM4590952"])), "single library")
    d = bed(RAW / "zhang2024_gse264498/GSE264498_YTHDF1_peaks_WT.bed.gz",
            {0: "chrom", 1: "start", 2: "end", 5: "strand"})
    source("YTHDF1", "zhang2024_HEK293T_CLIP", "human", "HEK293T", "cell line", 4, "CLIP-seq", "GSE264498",
           combine([d]), "deposited WT peaks (replicates merged)")
    w = pd.read_excel(RAW / "wang2015_gse63591/GSE63591_A-Y1-PARCLIP.xlsx", sheet_name=["rep1", "rep2"])
    reps = []
    for k in ("rep1", "rep2"):  # two listed clusters per gene
        parts = [pd.DataFrame({"chrom": w[k][c].astype(str).str.lower(), "start": w[k][s] - 1, "end": w[k][e],
                               "strand": w[k][st]})
                 for c, s, e, st in (("Chromosome", "ClusterStart", "ClusterEnd", "Strand"),
                                     ("Chromosome.1", "ClusterStart.1", "ClusterEnd.1", "Strand.1"))]
        reps.append(lift(pd.concat(parts).dropna(), "hg19"))
    source("YTHDF1", "wang2015_HeLa_PARCLIP", "human", "HeLa", "cell line", 4, "PAR-CLIP", "GSE63591", combine(reps),
           "clusters in both replicates, hg19 lifted to hg38")
    t = (pd.read_csv(utils.ROOT / "ref/target_sets/wang2014_ythdf2.tsv", sep="\t")
         .query("set == 'YTHDF2_Wang2014_HeLa_PARCLIP'"))
    source("YTHDF2", "wang2014_HeLa_PARCLIP", "human", "HeLa", "cell line", 4, "PAR-CLIP", "GSE49339",
           t.gene_symbol_current.str.upper().map(human).dropna().unique(),
           "published targets bound in all 3 replicates")
    source("YTHDF2", "encore2026_K562_eCLIP", "human", "K562", "cell line", 4, "eCLIP (PureCLIP)", "GSE315406",
           combine(library_peaks("clip_standard", "YTHDF2", "human", ["GSM9426917", "GSM9426919"])),
           "peaks in both replicates")
    source("YTHDF2", "lasman2020_mESC_eCLIP", "mouse", "mouse ESC", "cell line", 4, "eCLIP (PureCLIP)", "GSE148039",
           combine(library_peaks("clipflow", "YTHDF2", "mouse", ["GSM4590953"])), "single library")
    source("YTHDF2", "ito2022_CD4Tcell_iCLIP", "mouse", "CD4 T cells", "primary immune cells", 4, "iCLIP (PureCLIP)",
           "GSE188846", combine(library_peaks("clip_standard", "YTHDF2", "mouse", ["GSM5690950", "GSM5690951"])),
           "peaks in both libraries")
    source("YTHDF2", "gse230846_HEK293_iCLIP", "human", "HEK293", "cell line", 4, "iCLIP", "GSE230846",
           combine(library_peaks("clipflow", "YTHDF2", "human", ["GSM7245061", "GSM7245062"])),
           "peaks in both libraries")
    source("YTHDF3", "shi2017_HeLa_PARCLIP", "human", "HeLa", "cell line", 4, "PAR-CLIP", "GSE86214",
           combine(library_peaks("clipflow", "YTHDF3", "human", ["GSM2424844", "GSM2424845", "GSM2424846"])),
           "peaks in at least 2 of 3 libraries")
    source("YTHDF3", "lasman2020_mESC_eCLIP", "mouse", "mouse ESC", "cell line", 4, "eCLIP (PureCLIP)", "GSE148039",
           combine(library_peaks("clip_standard", "YTHDF3", "mouse", ["GSM4590954"])), "single library")


binding_sources()

# bound genes with m6A within 50 nt
maps = {("human", "brain"): nw.site_map("human"), ("mouse", "all"): nw.site_map("mouse")}
bound = []
for s in sources:
    species, features = s["species"], nw.features[s["species"]]
    if s["level"] == "peak":
        peaks = s["data"].reset_index(drop=True)
        n_peaks = nw.peak_genes(peaks, species).drop_duplicates(["peak", "gene_id"]).groupby("gene_id").size()
        table = pd.DataFrame({"gene_id": n_peaks.index, "n_peaks": n_peaks.to_numpy()})
        sites = {k[1]: np.expm1(nw.links_from_peaks(peaks, species, m)[0]).round()
                 for k, m in maps.items() if k[0] == species}
    else:
        table = pd.DataFrame({"gene_id": [g for g in pd.unique(s["data"]) if g in features.index], "n_peaks": np.nan})
        sites = {k[1]: m["load"] for k, m in maps.items() if k[0] == species}
    for name, count in sites.items():
        table[f"m6a_{name}"] = count.reindex(table.gene_id).fillna(0).astype(int).to_numpy()
    bound.append(table.assign(**{k: v for k, v in s.items() if k != "data"}, symbol=table.gene_id.map(features.symbol)))
    print(s["reader"], s["source"], len(table), "bound genes", flush=True)
bound = pd.concat(bound, ignore_index=True)
bound["m6a_primary"] = bound[["m6a_brain", "m6a_all"]].max(axis=1).fillna(0).astype(int)  # human: brain map; mouse: all
bound["human_gene_id"] = np.where(bound.species == "human", bound.gene_id,
                                  bound.symbol.str.upper().map(nw.to_human).map(human))
bound["T1"] = bound.m6a_primary > 0
bound = bound.drop(columns=["m6a_brain", "m6a_all"])
utils.write(bound, OUT / "bound_genes.tsv")

# IGF2BP2 tiers: T1 per source, T2 (m6A falls after METTL14 knockdown) and T3 (falls after IGF2BP2 loss)
t1 = bound[bound.T1 & (bound.reader == "IGF2BP2")].assign(tier="T1", binding_source=lambda d: d.source)
h = pd.read_excel(RAW / "huang2018_gse90642/GSE90642_shM14_vs_shCont_diff_gene_and_diff_m6A.xlsx")
h = h[h.Type.str.startswith("hypo")].drop_duplicates("geneName")
h = (h.assign(human_gene_id=h.geneName.str.upper().map(human)).dropna(subset=["human_gene_id"])
     .set_index("human_gene_id"))
t2 = t1[t1.human_gene_id.isin(h.index)].assign(tier="T2", source="huang2018_HepG2_shMETTL14_m6Aseq")
t2["m6a_log2fc_mettl14kd"], t2["expr_log2fc_mettl14kd"] = (h[c].reindex(t2.human_gene_id).to_numpy()
                                                          for c in ("m6A_log2Fc", "exp_log2Fc"))
t3 = []
loss = perturbation[(perturbation.reader == "IGF2BP2") & (perturbation.padj < 0.05) & (perturbation.log2fc < 0)
                    & (perturbation.direction != "OE")]
for dataset, d in loss.groupby("dataset"):
    d = d.set_index("gene_id")
    t = t1[(t1.species == d.species.iloc[0]) & t1.gene_id.isin(d.index)].assign(tier="T3", source=dataset)
    t3.append(t.assign(pert_log2fc=d.log2fc.reindex(t.gene_id).to_numpy(),
                       pert_padj=d.padj.reindex(t.gene_id).to_numpy()))
tiers = pd.concat([t1, t2, *t3], ignore_index=True)
tiers = tiers[["reader", "tier", "source", "binding_source", "species", "gene_id", "symbol", "human_gene_id",
               "m6a_primary", "m6a_log2fc_mettl14kd", "expr_log2fc_mettl14kd", "pert_log2fc", "pert_padj"]]
utils.write(tiers, OUT / "tiers_igf2bp2.tsv")
print(tiers.groupby(["tier", "species"]).gene_id.nunique().to_string(), flush=True)

# YTHDF sets: T1 (bound with m6A) and bound-unmarked genes per source, consensus and preferred set
sets = []


def add_set(frame, species, reader, name, primary):
    for tier, ids in (("T1", frame.gene_id[frame.T1]), ("bound_unmarked", frame.gene_id[~frame.T1])):
        sets.append(pd.DataFrame(dict(species=species, reader=reader, source=name, primary=primary, tier=tier,
                                      gene_id=ids.unique())))


def consensus(frame):
    """Genes bound in at least 2 sources: T1 when marked in at least 2, unmarked when marked in none."""
    bound_in, marked_in = frame.groupby("gene_id").source.nunique(), frame[frame.T1].groupby("gene_id").source.nunique()
    marked_in = marked_in.reindex(bound_in.index).fillna(0)
    keep = (marked_in >= 2) | ((bound_in >= 2) & (marked_in == 0))
    return pd.DataFrame(dict(gene_id=bound_in.index[keep], T1=(marked_in >= 2)[keep].to_numpy()))


for (species, reader), g in bound[bound.reader.isin(["YTHDF1", "YTHDF2", "YTHDF3"])].groupby(["species", "reader"]):
    vivo, neural = g[g.tissue_class == "brain in vivo"], g[g.tissue_class.str.startswith("neural")]
    if len(vivo):
        preferred = vivo.sort_values("rank").source.iloc[0]
    elif neural.source.nunique():
        preferred = "consensus_neural" if neural.source.nunique() > 1 else neural.source.iloc[0]
    else:
        preferred = "consensus_all" if g.source.nunique() > 1 else g.source.iloc[0]
    for name, rows in g.groupby("source"):
        add_set(rows, species, reader, name, preferred == name)
    for name, rows in (("consensus_all", g), ("consensus_neural", neural)):
        if rows.source.nunique() > 1:
            add_set(consensus(rows), species, reader, name, preferred == name)
mapped = (bound[(bound.source == "shi2026_striatum_CLIP") & bound.human_gene_id.notna()]
          .assign(gene_id=lambda d: d.human_gene_id).drop_duplicates("gene_id"))
add_set(mapped, "human", "YTHDF1", "shi2026_striatum_CLIP_mapped", False)  # the in vivo striatal CLIP by orthology
sets = pd.concat(sets, ignore_index=True)
utils.write(sets, OUT / "yth_binding_sets.tsv")
print(sets[sets.primary].groupby(["species", "reader", "source", "tier"]).size().to_string())
