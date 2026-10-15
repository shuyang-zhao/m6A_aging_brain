# Heterogeneity of the effect of a continuous gene property on the age effect (grf causal forest).
#
# Usage: Rscript causal_forest.R <workdir> <threads>
# Reads  <workdir>/data.tsv     y, w, x_* covariates (missing allowed), g_* gene classes (1, 0 or missing)
# Writes <workdir>/effects.tsv  class, quantity (all, in, out, difference), estimate, se
# all, in, out: average partial effect of w (AIPW) over all genes and inside or outside the class; difference: best
# linear projection of the gene-level effects on the class indicator (HC3)
suppressPackageStartupMessages(library(grf))

args <- commandArgs(TRUE)
workdir <- args[1]
data <- read.delim(file.path(workdir, "data.tsv"), check.names = FALSE)
X <- as.matrix(data[grep("^x_", names(data))])
# grf default of 2000 trees, fixed seed
forest <- causal_forest(X, data$y, data$w, num.trees = 2000, seed = 1, num.threads = as.integer(args[2]))

row <- function(class, quantity, value)
  data.frame(class = class, quantity = quantity, estimate = value[1], se = value[2])
rows <- list(row("all", "all", average_treatment_effect(forest)))
for (column in grep("^g_", names(data), value = TRUE)) {
  class <- sub("^g_", "", column)
  member <- data[[column]]
  known <- !is.na(member)
  blp <- best_linear_projection(forest, A = data.frame(member = member[known]), subset = known)
  rows <- c(rows, list(
    row(class, "in", average_treatment_effect(forest, subset = known & member == 1)),
    row(class, "out", average_treatment_effect(forest, subset = known & member == 0)),
    row(class, "difference", blp["member", c("Estimate", "Std. Error")])
  ))
}
write.table(do.call(rbind, rows), file.path(workdir, "effects.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
