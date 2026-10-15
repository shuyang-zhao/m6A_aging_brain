# Empirical-Bayes moderation of the per-gene joint-model variances (limma squeezeVar, robust).
#
# Usage: Rscript moderate.R <workdir>
# Reads  <workdir>/fit.tsv        feature_id, estimate, variance, df, unscaled
# Writes <workdir>/moderated.tsv  feature_id, se, t, df_total, p, ci_low, ci_high
#
# lmFit assumes one design for all features, but every gene has its own mRNA column and its own
# missing values. The OLS fits are therefore computed per gene and only the variances are
# moderated: the prior is estimated from the per-gene variances and their residual df.
suppressPackageStartupMessages(library(limma))

workdir <- commandArgs(TRUE)[1]
d <- read.delim(file.path(workdir, "fit.tsv"), check.names = FALSE,
                colClasses = c(feature_id = "character"))

v <- squeezeVar(d$variance, df = d$df, robust = TRUE)
d$se <- sqrt(v$var.post * d$unscaled)
d$t <- d$estimate / d$se
d$df_total <- pmin(d$df + v$df.prior, sum(d$df))
d$p <- 2 * pt(-abs(d$t), d$df_total)
half_width <- qt(0.975, d$df_total) * d$se
d$ci_low <- d$estimate - half_width
d$ci_high <- d$estimate + half_width

write.table(d[, c("feature_id", "se", "t", "df_total", "p", "ci_low", "ci_high")],
            file.path(workdir, "moderated.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
