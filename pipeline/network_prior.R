# Empirical-Bayes prior of the per-gene variances of a unit (limma squeezeVar, robust); one prior df for all genes, the
# median of the robust per-gene df.prior (moderate.R keeps them per gene)
#
# Usage: Rscript network_prior.R <workdir>
# Reads  <workdir>/variances.tsv  variance, df
# Writes <workdir>/prior.tsv      var_prior, df_prior
suppressPackageStartupMessages(library(limma))

workdir <- commandArgs(TRUE)[1]
d <- read.delim(file.path(workdir, "variances.tsv"))
v <- squeezeVar(d$variance, df = d$df, robust = TRUE)
write.table(data.frame(var_prior = v$var.prior, df_prior = median(v$df.prior)), file.path(workdir, "prior.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
