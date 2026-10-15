# Per-gene features of the canonical transcripts: length, exon structure, GC, K/R content, cell-type origin, half-life,
# protein lifetime, GO programme, site counts per modification (human m6A: also the m6A-seq peaks of the brain map,
# n_peaks_m6A_brain), DRACH motifs outside m6A sites; mouse also the
# hippocampal GLORI sites with summed stoichiometry, RiboNN translation efficiency and Orthrus components
# reads the gene models and site atlas (04) and the per-gene resources (06);
# writes output/09_transcript_features/{mouse,human}_features.tsv.gz
import re

import numpy as np
import pandas as pd

import utils

OUT = utils.OUT / "09_transcript_features"
SPECIES = {"mouse": "Mus musculus", "human": "Homo sapiens"}
N_ORTHRUS_PCS = 64
KR_WINDOW = 10  # residues
KR_CODONS = {"AAA", "AAG", "CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}
DE_CODONS = {"GAT", "GAC", "GAA", "GAG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}
DRACH = re.compile(r"(?=[AGT][AG]AC[ACT])")  # overlapping motifs; the methylated A is the third base
GLORI = ["m6a_hippocampus_cyto", "m6a_hippocampus_syn"]
BRAIN_TE = "brain Cerebellum Forebrain hippocampal Striatal_cells Primary_cortical_neurons dentate_gyrus".split()


def gc_fraction(sequences):
    return np.array([(text.count("G") + text.count("C")) / len(text) if text else np.nan for text in sequences])


def protein_features(cds):
    # K/R and D/E codons, stop codon removed
    codons = [cds[i : i + 3] for i in range(0, len(cds) - 2, 3)]
    while codons and codons[-1] in STOP_CODONS:
        codons.pop()
    if not codons:
        return np.nan, np.nan, np.nan
    kr = np.array([codon in KR_CODONS for codon in codons], float)
    window = np.convolve(kr, np.ones(KR_WINDOW), "valid") if len(kr) >= KR_WINDOW else kr.sum()
    return np.max(window), kr.mean(), (kr.sum() - sum(codon in DE_CODONS for codon in codons)) / len(codons)


origin = utils.read(utils.RESOURCES / "cell_type_origin" / "shares.parquet")

# young-brain protein lifetime (5 months), mean of cortex and cerebellum, by mouse symbol
lifetimes = utils.read(utils.RESOURCES / "kluever2022_sciadv" / "kluever_lifetimes.parquet")
uniprot = utils.read(utils.PREPARED / "ref" / "uniprot_gene_mouse.tsv").set_index("uniprot").symbol
lifetimes["symbol"] = lifetimes.accession.str.split("-").str[0].map(uniprot).str.upper()
lifetime = lifetimes.groupby(["symbol", "tissue"]).log2_t_young.mean().groupby("symbol").mean()
human_to_mouse = utils.orthologs().set_index("human").mouse

half_life = utils.read(utils.RESOURCES / "saluki_halflife" / "half_life.parquet")
programmes = utils.read(utils.RESOURCES / "go_release" / "programmes.parquet")
sites = utils.read(utils.RESOURCES / "site_atlas" / "sites.parquet")
# sites on the canonical transcript, reference base matched, counted once across maps; the human brain m6A-seq peaks
# are counted apart (the m6A load of human genes)
sites = sites[sites.gene_id.notna() & sites.tx_offset.notna() & sites.ref_ok]
# GLORI sites of mouse hippocampus with the mean stoichiometry of the two maps
glori = sites[sites["map"].isin(GLORI)].groupby(["gene_id", "tx_offset"]).stoich.mean().groupby("gene_id")
sites = sites.drop_duplicates(["species", "modification", "chrom", "start", "end", "strand"])
peaks = sites[sites["map"] == utils.BRAIN_M6A]
brain_m6a = peaks.groupby("gene_id").size()
sites = sites[sites["map"] != utils.BRAIN_M6A]

for species, latin in SPECIES.items():
    folder = utils.RESOURCES / f"gencode_{species}"
    tx = utils.read(folder / "transcripts.parquet")
    exons = utils.read(folder / "exons.parquet")
    seq = utils.read(folder / "sequences.parquet").set_index("tx_id").seq.reindex(tx.tx_id).to_numpy()
    assert tx.gene_id.is_unique and (tx.cds_len > 0).all()
    assert (tx.tx_len.to_numpy() == [len(text) for text in seq]).all()
    # gene-map symbol, GENCODE name where missing
    symbol = tx.gene_id.map(utils.ensg_to_symbol(species)).fillna(tx.gene_name)
    features = pd.DataFrame({"gene_id": tx.gene_id, "symbol": symbol, "chrom": tx.chrom})
    # length and exon structure
    features["log_tx_len"] = np.log(tx.tx_len)
    features["log1p_utr5_len"] = np.log1p(tx.utr5_len)
    features["log_cds_len"] = np.log(tx.cds_len)
    features["log1p_utr3_len"] = np.log1p(tx.utr3_len)
    sizes = exons.sort_values(["tx_id", "exon_no"]).groupby("tx_id")["len"]
    first, last = tx.tx_id.map(sizes.first()), tx.tx_id.map(sizes.last())
    internal = ((tx.tx_id.map(sizes.sum()) - first - last) / (tx.n_exons - 2)).where(tx.n_exons > 2, 0)
    features["log_n_exons"] = np.log(tx.n_exons)
    features["log1p_internal_exon_mean_len"] = np.log1p(internal)
    features["log_first_exon_len"], features["log_last_exon_len"] = np.log(first), np.log(last)
    features["junction_density_per_kb"] = (tx.n_exons - 1) / (tx.tx_len / 1000)
    # genomic span of the gene and intronic length of the canonical transcript
    features["log_gene_len"] = np.log(tx.gene_end - tx.gene_start + 1)
    features["log1p_intron_len"] = np.log1p(tx.tx_end - tx.tx_start + 1 - tx.tx_len)
    # GC of the transcript and its regions
    utr5 = [text[:cut] for text, cut in zip(seq, tx.cds_start, strict=True)]
    cds = [text[cut:stop] for text, cut, stop in zip(seq, tx.cds_start, tx.cds_end, strict=True)]
    utr3 = [text[stop:] for text, stop in zip(seq, tx.cds_end, strict=True)]
    for name, region in {"tx": seq, "utr5": utr5, "cds": cds, "utr3": utr3}.items():
        features[f"gc_{name}"] = gc_fraction(region)
    # coding sequence: K/R content and net charge
    protein = np.array([protein_features(text) for text in cds])
    features["kr_max10"], features["kr_frac"], features["net_charge_per_res"] = protein.T
    # cell-type origin, half-life, protein lifetime (human via mouse ortholog)
    features = features.join(origin[origin.species == species].filter(like="ct_frac_"), on="gene_id")
    half = half_life[half_life.species == latin].set_index("gene_id").hl_pc1
    features["mrna_half_life_pc1"] = features.gene_id.map(half)
    key = features.symbol.str.upper()
    features["log2_protein_lifetime"] = (key if species == "mouse" else key.map(human_to_mouse)).map(lifetime)
    # functional programme: GO myelin, lysosome, mitochondrion and synapse
    for name in ["myelin", "lysosome", "mitochondrion", "synapse"]:
        members = programmes.symbol[(programmes.species == latin) & (programmes.programme == name)]
        features[name] = features.symbol.isin(members).astype(int)
    # sites per gene, one column per modification
    per_gene = sites[sites.species == latin].groupby(["modification", "gene_id"]).size()
    for mod in utils.MODIFICATIONS:
        features[f"n_sites_{mod}"] = features.gene_id.map(per_gene[mod]).fillna(0).astype(int)
    # DRACH motifs whose A is not an m6A site (mouse: the single-base sites; human: outside the brain m6A-seq peaks)
    if species == "mouse":
        marked = sites[(sites.species == latin) & (sites.modification == "m6A")].groupby("tx_id").tx_offset
        marked = {tx_id: set(offsets.astype(int)) for tx_id, offsets in marked}
    else:
        radius = ((peaks.end - peaks.start) / 2).to_numpy()
        intervals = pd.DataFrame(dict(tx_id=peaks.tx_id, low=peaks.tx_offset - radius, high=peaks.tx_offset + radius))
        intervals = {tx_id: group[["low", "high"]].to_numpy() for tx_id, group in intervals.groupby("tx_id")}
    drach, unmarked = [], []
    for tx_id, text in zip(tx.tx_id, seq, strict=True):
        a = np.array([m.start() + 2 for m in DRACH.finditer(text)])
        if species == "mouse":
            free = ~np.isin(a, list(marked.get(tx_id, ())))
        else:
            spans = intervals.get(tx_id, np.empty((0, 2)))
            free = ~((a[:, None] >= spans[:, 0]) & (a[:, None] <= spans[:, 1])).any(axis=1)
        drach.append(len(a))
        unmarked.append(int(free.sum()))
    features["n_drach"], features["n_drach_unmarked"] = drach, unmarked
    if species == "human":
        features["n_peaks_m6A_brain"] = features.gene_id.map(brain_m6a).fillna(0).astype(int)
    else:
        features["n_sites_m6A_glori"] = features.gene_id.map(glori.size()).fillna(0).astype(int)
        features["stoich_m6A_glori"] = features.gene_id.map(glori.sum()).fillna(0)
        # RiboNN translation efficiency, mean and brain tissues
        te = utils.read(utils.RESOURCES / "ribonn_mouse" / "predictions.parquet").set_index("tx_id")
        features["ribonn_te_mean"] = tx.tx_id.map(te.mean_predicted_TE).to_numpy()
        brain = te[[f"predicted_TE_{tissue}" for tissue in BRAIN_TE]].mean(axis=1)
        features["ribonn_te_brain"] = tx.tx_id.map(brain).to_numpy()
        # leading principal components of the standardised Orthrus embeddings
        orthrus = utils.read(utils.RESOURCES / "orthrus_mouse" / "embeddings.parquet")
        assert set(orthrus.gene_id) <= set(features.gene_id)
        scaled = orthrus.filter(like="emb")
        scaled = ((scaled - scaled.mean()) / (scaled.std(ddof=0) + 1e-8)).to_numpy()
        scores = scaled @ np.linalg.svd(scaled, full_matrices=False)[2][:N_ORTHRUS_PCS].T
        columns = [f"orthrus_pc{i + 1}" for i in range(N_ORTHRUS_PCS)]
        features = features.join(pd.DataFrame(scores, index=orthrus.gene_id, columns=columns), on="gene_id")
    utils.write(features, OUT / f"{species}_features.tsv.gz")
    print(f"{species}: {len(features)} genes, {features.shape[1] - 3} features")
