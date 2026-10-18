# Targets of the m6A readers: transcripts with an m6A site (mouse, GLORI) or brain m6A-seq peak (human) within 50 nt of
# the reader's CLIP peaks, interval to interval (target), against the other m6A-marked transcripts of the species
# (comparison)
# reads the site atlas, the reader CLIP peaks and the regulator table; writes output/40_target_sets/target_genes.tsv.gz
import pandas as pd

import utils

OUT = utils.OUT / "40_target_sets"
SPECIES = {"Homo sapiens": "human", "Mus musculus": "mouse"}

columns = ["modification", "species", "gene_id", "tx_offset", "ref_ok", "chrom", "start", "end", "strand"]
atlas = utils.read(utils.RESOURCES / "site_atlas" / "sites.parquet", columns=columns)
sites = atlas[atlas.gene_id.notna() & atlas.tx_offset.notna() & atlas.ref_ok]
sites = sites[sites.modification == "m6A"].assign(species=lambda d: d.species.map(SPECIES))
regulators = utils.read(utils.ROOT / "ref/gene_sets/rna_mod_regulators.tsv")
readers = regulators.query("role == 'reader' and modification == 'm6A'")
peaks = utils.read(utils.RESOURCES / "reader_clip" / "peaks.parquet")
members = []
for (reader, species), reader_peaks in peaks[peaks.reader.isin(readers.human_symbol)].groupby(["reader", "species"]):
    marked = sites[sites.species == species]
    target = set(marked.gene_id[utils.near_peaks(marked, reader_peaks, utils.PEAK_WINDOW)])
    print(reader, species, f"{len(target)} target genes")
    for group, genes in (("target", target), ("comparison", set(marked.gene_id) - target)):
        key = dict(modification="m6A", regulator=reader, species=species, window=utils.PEAK_WINDOW)
        members.append(pd.DataFrame({"gene_id": sorted(genes), "group": group}).assign(**key))
targets = pd.concat(members)[["modification", "regulator", "species", "window", "gene_id", "group"]]
utils.write(targets, OUT / "target_genes.tsv.gz")
