# Sharing of regulator-target edges across units: mashr on the edge estimates with their permutation SDs.
#
# Usage: Rscript network_mash.R <workdir>
# Reads  <workdir>/bhat.tsv, shat.tsv   edges (rows) x units (columns); absent edges: Bhat 0, Shat 100
# Writes <workdir>/mash.tsv             edge, unit, posterior_mean, posterior_sd, lfsr
# Covariances: canonical plus data-driven (PCA of the strong edges)
suppressPackageStartupMessages(library(mashr))

workdir <- commandArgs(TRUE)[1]
read_matrix <- function(f) as.matrix(read.delim(file.path(workdir, f), row.names = 1, check.names = FALSE))
bhat <- read_matrix("bhat.tsv")
shat <- read_matrix("shat.tsv")

data <- mash_set_data(bhat, shat)
covs <- cov_canonical(data)
strong <- get_significant_results(mash_1by1(data), thresh = 0.05)
if (length(strong) >= 5) {
  strong_data <- mash_set_data(bhat[strong, , drop = FALSE], shat[strong, , drop = FALSE])
  pca <- cov_pca(strong_data, min(5, ncol(bhat), length(strong)))
  covs <- c(covs, pca)
}
fit <- mash(data, covs, verbose = FALSE)

long <- function(m, name)
  setNames(data.frame(rownames(m)[row(m)], colnames(m)[col(m)], as.vector(m)), c("edge", "unit", name))
out <- Reduce(function(a, b) merge(a, b, by = c("edge", "unit")),
              list(long(get_pm(fit), "posterior_mean"), long(get_psd(fit), "posterior_sd"),
                   long(get_lfsr(fit), "lfsr")))
write.table(out, file.path(workdir, "mash.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
