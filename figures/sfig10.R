# Additional file 1 Fig. S10: what drives the clock in disease, the disease coefficient split exactly over groups of
# genes (the partial scores of the groups sum to the score)
# a, every study by cell source and by m6A sites; b, every region of Merrihew 2023 and Blumenreich 2024; c, cell source
# crossed with m6A sites and with the sign of the weight, AD combined and HD
# reads source_data/sfig10a-c (script 53); writes output/figures/sfig10.png and .pdf
source("figures/style.R")

SOURCES <- c(neuron = "Neuron", astrocyte = "Astrocyte", oligodendrocyte = "Oligodendrocyte", microglia = "Microglia",
             ubiquitous = "Ubiquitous", "no cell-type data" = "No cell-type data")
SITES <- c("0" = "0", "1" = "1", "2-3" = "2–3", "≥4" = "≥4")
STUDIES <- c("Ping 2020", "Merrihew 2023", "AD combined", "Blumenreich 2024", "Umoh 2018", "Langfelder 2016")
CONTRIBUTION <- "Contribution to the disease coefficient (SD of the clock score; HD per 100 CAG repeats)"
EXPECTED <- "Expected from the share of the weights"
studies <- read_panel("sfig10a")
total <- studies[partition == "total", setNames(signed(est), study)]
study_strip <- function(study, disease, by = "") {
  s <- paste0(fifelse(study == "AD combined", study, paste0(study, ", ", disease)), by, "\ntotal ", total[study])
  factor(s, unique(s[order(match(study, STUDIES))]))
}

# a, b: contribution of each group (bar and 95% CI; red, P < 0.05) and the contribution expected from its share of the
# weights (tick, total coefficient x weight share)
groups <- function(d) {
  d <- d[partition != "total"]
  d[, `:=`(row = factor(c(SOURCES, SITES)[group], rev(c(SOURCES, SITES))),
           partition = factor(partition, c("cell source", "m6A sites"), c("Cell source", "m6A sites")))]
}
contribution <- function(d, x, scales = "free") {
  ggplot(d, aes(est, row)) +
    ZERO_X +
    geom_col(aes(fill = p < 0.05), width = 0.66, colour = EDGE, linewidth = 0.5) +
    geom_linerange(aes(xmin = ci_low, xmax = ci_high), linewidth = 0.6, colour = EDGE) +
    geom_errorbar(aes(xmin = expected, xmax = expected, linetype = EXPECTED), width = 0.86, linewidth = 1,
                  colour = INK, key_glyph = "vline") +
    facet_grid(partition ~ strip, scales = scales, space = "free_y") +
    scale_fill_manual(values = c("TRUE" = RED, "FALSE" = "white"), breaks = c("TRUE", "FALSE"),
                      labels = c("TRUE" = "P < 0.05", "FALSE" = "P ≥ 0.05"), name = NULL) +
    scale_linetype_manual(values = setNames("solid", EXPECTED), name = NULL) +
    scale_x_continuous(labels = MINUS_LABELS, breaks = scales::breaks_pretty(3), expand = expansion(mult = 0.1)) +
    guides(fill = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
    labs(x = x, y = NULL) +
    theme_rows() +
    BOXED +
    theme(strip.text.y = element_text(angle = 270), panel.spacing.x = unit(4, "mm"))
}

a <- groups(studies)[, strip := study_strip(study, disease)]
a <- contribution(a, CONTRIBUTION) + theme(legend.position = "bottom")

# b: per region on one axis per study, the study named in the axis title
REGIONS <- c("Superior and middle temporal gyri (SMTG)" = "SMTG", "Inferior parietal lobule (IPL)" = "IPL",
             "Middle temporal gyrus" = "Middle temporal\ngyrus")
b <- read_panel("sfig10b")
region_total <- b[partition == "total", setNames(signed(est), region)]
b <- groups(b)
b[, strip := paste0(fifelse(region %in% names(REGIONS), REGIONS[region], region), "\ntotal ", region_total[region])]
b[, strip := factor(strip, unique(strip))]
by_study <- function(s, first) {
  d <- b[study == s]
  contribution(d, sprintf("Contribution in %s, %s (SD of the clock score)", s, d$disease[1]), "free_y") +
    theme(legend.position = "none", panel.spacing.x = unit(7, "mm"),
          axis.text.y = if (first) element_text() else element_blank(),
          strip.text.y = if (first) element_blank() else element_text(angle = 270))
}
b <- (by_study("Merrihew 2023", TRUE) + by_study("Blumenreich 2024", FALSE) + plot_layout(widths = c(4, 5))) + FLUSH

# c: the crosses, cell source by m6A sites and by the sign of the weight (genes predicted to decline or rise with age),
# named in the strips; bar, CI and tick as in a
cross <- read_panel("sfig10c")
cross[, row := factor(SOURCES[cell_source], rev(SOURCES))]
crossed <- function(d, colours, name, y_labels = TRUE) {
  dodge <- position_dodge(0.84, preserve = "single", reverse = TRUE)
  ggplot(d, aes(est, row, fill = subgroup)) +
    ZERO_X +
    geom_col(width = 0.84, colour = EDGE, linewidth = 0.4, position = dodge) +
    geom_linerange(aes(xmin = ci_low, xmax = ci_high), linewidth = 0.5, colour = EDGE, position = dodge) +
    geom_errorbar(aes(xmin = expected, xmax = expected), width = 0.94, linewidth = 0.9, colour = INK,
                  position = dodge) +
    facet_grid(. ~ strip, scales = "free_x") +
    scale_fill_manual(values = colours, name = name) +
    scale_x_continuous(labels = MINUS_LABELS, breaks = scales::breaks_pretty(3), expand = expansion(mult = 0.08)) +
    labs(x = CONTRIBUTION, y = NULL) +
    theme_rows() +
    BOXED +
    theme(legend.position = "bottom", panel.spacing.x = unit(4, "mm"),
          axis.text.y = if (y_labels) element_text() else element_blank())
}
sites <- cross[partition == "cell source x m6A sites"][, subgroup := factor(SITES[subgroup], SITES)]
sites[, strip := study_strip(study, disease, ", by sites")]
direction <- cross[partition == "cell source x direction"][, strip := study_strip(study, disease, ", by direction")]
direction[, subgroup := factor(subgroup, c("declining", "rising"), c("Declining (weight < 0)", "Rising (weight > 0)"))]
cc <- (crossed(sites, setNames(c(GREY, FOCAL_LIGHT, FOCAL, FOCAL_DARK), SITES), "m6A sites") +
         crossed(direction, setNames(c(FOCAL, RED), levels(direction$subgroup)), NULL, FALSE) +
         plot_layout(axis_titles = "collect")) + FLUSH

fig <- panel_grid(list(a = a, b = b, c = cc), ncol = 1, heights = c(100, 100, 110))
save_fig(fig, "sfig10", 310)
