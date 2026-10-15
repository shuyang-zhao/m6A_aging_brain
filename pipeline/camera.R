# Competitive gene set test of a gene-level statistic (limma cameraPR, default inter-gene correlation).
#
# Usage: Rscript camera.R <workdir>
# Reads  <workdir>/stat.tsv   gene_id, t
#        <workdir>/sets.tsv   set, gene_id
# Writes <workdir>/camera.tsv set, n_genes, direction, p
suppressPackageStartupMessages(library(limma))

workdir <- commandArgs(TRUE)[1]
stat <- read.delim(file.path(workdir, "stat.tsv"), colClasses = c(gene_id = "character"))
sets <- read.delim(file.path(workdir, "sets.tsv"), colClasses = c(set = "character", gene_id = "character"))

t <- setNames(stat$t, stat$gene_id)
index <- lapply(split(sets$gene_id, sets$set), function(genes) which(names(t) %in% genes))
res <- cameraPR(t, index, sort = FALSE)

out <- data.frame(set = rownames(res), n_genes = res$NGenes, direction = res$Direction, p = res$PValue)
write.table(out, file.path(workdir, "camera.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
