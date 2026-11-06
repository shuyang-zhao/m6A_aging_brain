# Loads the panel library; run from the repository root
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
for (f in c("theme_pub.R", "palettes.R", "grammar.R", "recipes.R", "compose.R")) {
  sys.source(file.path("figures/pub", f), envir = globalenv())
}
