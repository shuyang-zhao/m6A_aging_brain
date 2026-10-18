# CLIP peaks of the m6A readers, one set per reader and species
# reads data/raw/clip_standard (PureCLIP binding regions of nf-core/clipseq, made by extra/clip_peaks.sh) and the
# regulator table ref/gene_sets/rna_mod_regulators.tsv
# writes output/prepared/resources/reader_clip/peaks.parquet and provenance.tsv (the provenance of every library)
import pandas as pd

import utils

PEAKS = utils.DATA / "raw/clip_standard"
COLUMNS = ["chrom", "start", "end", "name", "score", "strand"]

# m6A readers of the regulator table (the CLIP folder also holds m5C readers)
readers = set(
    utils.read(utils.ROOT / "ref/gene_sets/rna_mod_regulators.tsv")
    .query("role == 'reader' and modification == 'm6A'")
    .human_symbol
)

# bed file of each IP library, named {species}_{reader}_{library}.bed.gz
parts, libraries = [], []
for path in sorted(PEAKS.glob("*.bed.gz")):
    species, reader, library = path.name.removesuffix(".bed.gz").split("_")
    if reader not in readers:
        continue
    bed = pd.read_csv(path, sep="\t", header=None, names=COLUMNS, usecols=["chrom", "start", "end", "strand"])
    parts.append(bed.assign(reader=reader, species=species))
    libraries.append((species, reader, library))

# peaks of the libraries of a reader pooled; the libraries are those of the per-reader provenance tables
provenance = pd.concat([utils.read(path) for path in sorted(PEAKS.glob("*.provenance.tsv"))], ignore_index=True)
provenance = provenance[provenance.reader.isin(readers)]
assert sorted(libraries) == sorted(zip(provenance.species, provenance.reader, provenance.library, strict=True))
peaks = pd.concat(parts, ignore_index=True)[["reader", "species", "chrom", "start", "end", "strand"]]
peaks = peaks.drop_duplicates().sort_values(["species", "reader", "chrom", "start", "end", "strand"], ignore_index=True)
n_peaks = peaks.groupby(["species", "reader"]).size()
utils.write(peaks, utils.RESOURCES / "reader_clip/peaks.parquet")
utils.write(provenance, utils.RESOURCES / "reader_clip/provenance.tsv")
print(n_peaks.to_string())
print(f"{len(peaks)} peaks of {len(n_peaks)} reader-species pairs")
