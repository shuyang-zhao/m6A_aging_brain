# Two-group comparison of a gene x sample matrix (script 37).
#
# Usage: Rscript perturbation_de.R <mode> <matrix.tsv> <groups.tsv> <out.tsv>
# matrix.tsv  gene x sample; counts (edger, te) or log2 values (limma, te_limma)
# groups.tsv  sample, group (ref or test) and, for the te modes, assay (rna or rpf)
# mode        edger: edgeR quasi-likelihood; limma: limma-trend; te, te_limma: the assay x group interaction
#             (translation efficiency change) with the same two engines
# out.tsv     gene, log2fc, padj
suppressPackageStartupMessages({library(edgeR); library(limma)})

a <- commandArgs(TRUE)
m <- as.matrix(read.delim(a[2], row.names = 1, check.names = FALSE))
g <- read.delim(a[3])
g <- g[match(colnames(m), g$sample), ]
group <- factor(g$group, levels = c("ref", "test"))
te <- startsWith(a[1], "te")
if (te) {
  design <- model.matrix(~ factor(g$assay, levels = c("rna", "rpf")) * group)
  coef <- ncol(design)  # assay x group
} else {
  design <- model.matrix(~group)
  coef <- 2
}

if (a[1] %in% c("edger", "te")) {
  y <- DGEList(m)
  keep <- if (te) filterByExpr(y, design = design, min.count = 10, min.total.count = 30) else filterByExpr(y, design)
  y <- estimateDisp(calcNormFactors(y[keep, , keep.lib.sizes = FALSE]), design)
  tt <- topTags(glmQLFTest(glmQLFit(y, design), coef = coef), n = Inf)$table
  res <- data.frame(gene = rownames(tt), log2fc = tt$logFC, padj = tt$FDR)
} else {
  if (!te) m <- m[rowMeans(m) > 0.5, ]
  tt <- topTable(eBayes(lmFit(m, design), trend = TRUE), coef = coef, number = Inf)
  res <- data.frame(gene = rownames(tt), log2fc = tt$logFC, padj = tt$adj.P.Val)
}
write.table(res, a[4], sep = "\t", quote = FALSE, row.names = FALSE)
