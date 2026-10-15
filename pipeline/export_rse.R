# Counts, row ranges and sample table of a SummarizedExperiment .rda as plain files
#
# Usage: Rscript export_rse.R <file.rda> <outdir>
# Writes <outdir>/counts.int32       integer counts, one sample after another
#        <outdir>/rows.tsv.gz        row ranges and row data
#        <outdir>/coldata.tsv.gz     sample table
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(data.table)
})

args <- commandArgs(TRUE)
outdir <- args[2]
env <- new.env()
object <- get(load(args[1], envir = env)[1], env)

# integer counts, written column by column
counts <- assay(object, "counts")
storage.mode(counts) <- "integer"
writeBin(as.vector(counts), file.path(outdir, "counts.int32"), size = 4)

# list columns become the mean of numbers or the distinct values
flat <- function(column) {
  if (!(is(column, "List") || is.list(column))) return(as.vector(column))
  vapply(as.list(column), function(values) {
    values <- unlist(values)
    if (is.numeric(values)) as.character(mean(values)) else paste(unique(values), collapse = ";")
  }, character(1))
}

ranges <- rowRanges(object)
rows <- data.frame(
  row_id = names(ranges), seqnames = as.character(seqnames(ranges)), start = start(ranges),
  end = end(ranges), strand = as.character(strand(ranges)), lapply(mcols(ranges), flat), check.names = FALSE
)
fwrite(rows, file.path(outdir, "rows.tsv.gz"), sep = "\t")
coldata <- data.frame(sample = colnames(object), lapply(colData(object), flat), check.names = FALSE)
fwrite(coldata, file.path(outdir, "coldata.tsv.gz"), sep = "\t")
