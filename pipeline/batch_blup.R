# Remove a random batch effect (BLUP) from protein, gene by gene, protecting age, own mRNA and
# the fixed covariates:  protein ~ age + rna + covariates + (1 | batch),  REML, lme4.
#
# Usage: Rscript batch_blup.R <workdir> <cores>
# Reads  <workdir>/protein.tsv, rna.tsv   genes x samples (first column feature_id)
#        <workdir>/design.tsv             Intercept, covariates and age, one row per sample
#        <workdir>/samples.tsv            sample_id, batch
# Writes <workdir>/corrected.tsv          batch-corrected protein, same layout as protein.tsv
#        <workdir>/batch_audit.tsv        feature_id, status
# status is "ok", "insufficient" (too few observations; the row is returned as NA) or the
# message of an lme4 error.
suppressPackageStartupMessages({
  library(lme4)
  library(parallel)
})

args <- commandArgs(TRUE)
workdir <- args[1]
cores <- as.integer(args[2])
P <- as.matrix(read.delim(file.path(workdir, "protein.tsv"), row.names = 1, check.names = FALSE))
R <- as.matrix(read.delim(file.path(workdir, "rna.tsv"), row.names = 1, check.names = FALSE))
D <- read.delim(file.path(workdir, "design.tsv"), check.names = FALSE)
S <- read.delim(file.path(workdir, "samples.tsv"), check.names = FALSE,
                colClasses = c(sample_id = "character"))
stopifnot(identical(colnames(P), S$sample_id), identical(dimnames(P), dimnames(R)))

C <- D[, !(colnames(D) %in% c("Intercept", "age")), drop = FALSE]
names(C) <- paste0("V", seq_len(ncol(C)))
base <- cbind(data.frame(age = D$age, batch = factor(S$batch)), C)
f <- as.formula(paste("protein ~", paste(c("age", "rna", names(C), "(1|batch)"), collapse = " + ")))

one <- function(i) {
  d <- cbind(base, protein = as.numeric(P[i, ]), rna = as.numeric(R[i, ]))
  ok <- complete.cases(d)
  if (sum(ok) < ncol(C) + 5) return(list(y = rep(NA_real_, ncol(P)), status = "insufficient"))
  fit <- tryCatch(
    suppressMessages(suppressWarnings(
      lmer(f, data = d[ok, ], REML = TRUE, control = lmerControl(check.conv.singular = "ignore"))
    )),
    error = function(e) e
  )
  if (inherits(fit, "error")) return(list(y = rep(NA_real_, ncol(P)), status = conditionMessage(fit)))
  blup <- ranef(fit)$batch[, 1]
  names(blup) <- rownames(ranef(fit)$batch)
  y <- P[i, ]
  y[ok] <- y[ok] - blup[as.character(d$batch[ok])]
  list(y = y, status = "ok")
}

res <- mclapply(seq_len(nrow(P)), one, mc.cores = cores)
out <- do.call(rbind, lapply(res, `[[`, "y"))
dimnames(out) <- dimnames(P)
write.table(data.frame(feature_id = rownames(out), out, check.names = FALSE),
            file.path(workdir, "corrected.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
audit <- data.frame(feature_id = rownames(P), status = vapply(res, `[[`, "", "status"))
write.table(audit, file.path(workdir, "batch_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
