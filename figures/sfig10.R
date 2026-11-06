# Fig. S10: what drives the clock in disease, the age acceleration split exactly over groups of genes (the partial
# predictions of the groups sum to the prediction)
# a, every study by cell source and by m6A sites; b, every region of Merrihew 2023 and Blumenreich 2024; c, cell source
# crossed with m6A sites and with the sign of the weight, AD combined and HD
# reads source_data/sfig10a-c (script 53); writes output/figures/sfig10
source("figures/style.R")

SOURCES <- c(neuron = "Neuron", astrocyte = "Astrocyte", oligodendrocyte = "Oligodendrocyte", microglia = "Microglia",
             ubiquitous = "Ubiquitous", "no cell-type data" = "No cell-type data")
SITE_NAMES <- c("0" = "0", "1" = "1", "2-3" = "2–3", "≥4" = "≥4")
GROUPS <- c("all genes" = "All genes", SOURCES, SITE_NAMES)
STUDIES <- c("Ping 2020", "Merrihew 2023", "AD combined", "Blumenreich 2024", "Umoh 2018", "Langfelder 2016")
BAR <- 0.7
P <- function(p) patchwork::wrap_elements(full = p) # panels keep their own widths, so the grid columns stay equal

# a, b: contribution to the acceleration of every group (colour and number; bold, P < 0.05) next to the contribution
# expected from its share of the weights (in parentheses); the total of all genes on top
partition_rows <- function(d) {
  d[, `:=`(row = factor(GROUPS[group], rev(GROUPS)),
           partition = factor(partition, c("total", "cell source", "m6A sites"), c(" ", "Cell source", "m6A sites")))]
}
long <- function(d, by) {
  keep <- c(by, "partition", "group", "digits")
  rbind(d[, c(.SD, list(kind = "est", value = est, p = p)), .SDcols = keep],
        d[partition != "total", c(.SD, list(kind = "exp", value = expected, p = NA_real_)), .SDcols = keep])
}
heat_table <- function(d, columns, legend = FALSE) {
  d[, text := fifelse(kind == "est", mapply(number, value, digits), paste0("(", mapply(number, value, digits), ")"))]
  ggplot(d, aes(kind, row, fill = value)) +
    geom_tile(colour = "white", linewidth = lw(0.6)) +
    geom_text(aes(label = text, fontface = fifelse(!is.na(p) & p < 0.05, "bold", "plain")), size = pt_size(6)) +
    ggh4x::facet_nested(rows = vars(partition), cols = columns, scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_signed(6, "Contribution to the acceleration (years; HD, months per 100 CAG repeats)", light = TRUE) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_tile() +
    theme(axis.text.y = element_text(size = 6), axis.text.x = element_blank(), strip.placement = "outside",
          strip.text = element_text(size = 6, lineheight = 0.95), strip.text.y.left = element_text(angle = 90),
          panel.spacing.x = unit(1.6, "mm"), panel.spacing.y = unit(1.2, "mm"),
          legend.position = if (legend) "bottom" else "none", legend.margin = margin(t = 1, unit = "mm"))
}

# a: the six studies
a <- read_panel("sfig10a")
HEADERS <- c("Ping 2020" = "Ping 2020\nAD", "Merrihew 2023" = "Merrihew 2023\nAD", "AD combined" = "AD combined\n",
             "Blumenreich 2024" = "Blumenreich 2024\nPD", "Umoh 2018" = "Umoh 2018\nALS", "Langfelder 2016" = "Langfelder 2016\nHD, per 100 CAG")
a[, digits := fifelse(study == "Langfelder 2016", 2, 1)]
a <- partition_rows(long(a, "study"))
a[, study := factor(HEADERS[study], HEADERS[STUDIES])]
a <- heat_table(a, vars(study))

# b: the regions of Merrihew 2023 and Blumenreich 2024
REGIONS <- c("Superior and middle temporal gyri (SMTG)" = "SMTG", "Inferior parietal lobule (IPL)" = "IPL",
             Hippocampus = "Hippocampus", Caudate = "Caudate", Striatum = "Striatum", "Substantia nigra" = "Substantia\nnigra",
             "Cingulate gyrus" = "Cingulate\ngyrus", "Middle temporal gyrus" = "Middle\ntemporal gyrus",
             "Occipital cortex" = "Occipital\ncortex")
b <- read_panel("sfig10b")
b[, `:=`(digits = 1, study = paste0(study, ", ", disease), region = unname(REGIONS[region]))]
b <- partition_rows(long(b, c("study", "region")))
b[, `:=`(study = factor(study, unique(study)), region = factor(region, unique(region)))]
b <- heat_table(b, vars(study, region), legend = TRUE)

# c: cell source crossed with m6A sites (stacked bars, the number of sites in the ramp of the cohort set) and with the
# sign of the weight (bar in the colour of the cell source, 95% CI, tick: expected from the weight share; genes predicted to
# decline with age to the left, to rise to the right)
cross <- read_panel("sfig10c")
cross[, row := factor(SOURCES[cell_source], rev(SOURCES))]
SIDES <- c("m6A sites", "Declining\n(weight < 0)", "Rising\n(weight > 0)")
crossed <- function(s, hue, limits, breaks, x_label) {
  sites <- cross[study == s & partition == "cell source x m6A sites"]
  sites[, `:=`(side = factor(SIDES[1], SIDES), subgroup = factor(SITE_NAMES[subgroup], SITE_NAMES))]
  weight <- cross[study == s & partition == "cell source x direction"]
  weight[, `:=`(side = factor(fifelse(subgroup == "declining", SIDES[2], SIDES[3]), SIDES), fill = CELL[SOURCES[cell_source]])]
  ggplot(mapping = aes(y = droplevels(row))) +
    geom_vline(xintercept = 0, linewidth = lw(0.5)) +
    geom_col(data = sites, aes(est, fill = subgroup), width = BAR, position = position_stack(reverse = TRUE)) +
    geom_col(data = weight, aes(est, fill = I(fill)), width = BAR) +
    geom_errorbar(data = weight, aes(xmin = ci_low, xmax = ci_high), width = 0, linewidth = lw(0.5), colour = EDGE) +
    geom_errorbar(data = weight, aes(xmin = expected, xmax = expected), width = BAR, linewidth = lw(0.9)) +
    ggh4x::facet_nested(. ~ study + side, scales = "free_x", independent = "x") +
    ggh4x::facetted_pos_scales(x = list(scale_x_continuous(labels = MINUS_LABELS, breaks = scales::breaks_pretty(3)),
                                        scale_x_reverse(limits = rev(limits), labels = MINUS_LABELS, breaks = breaks),
                                        scale_x_continuous(limits = limits, labels = MINUS_LABELS, breaks = breaks))) +
    ggh4x::force_panelsizes(cols = unit(c(30, 19, 19), "mm")) +
    scale_fill_manual(values = sites_ramp(hue), name = NULL, breaks = SITE_NAMES) +
    labs(x = x_label, y = NULL) +
    theme_box() +
    theme(strip.text = element_text(size = 6, lineheight = 0.95), axis.ticks.y = element_blank(),
          panel.spacing.x = unit(1.2, "mm"), legend.position = "top", legend.justification = "left",
          legend.key.size = unit(2.4, "mm"), legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0.8, "mm"))
}
c_ <- P(crossed("AD combined", HUMAN, c(-5.5, 7.2), c(0, 5), "Contribution (years)")) |
  P(crossed("Langfelder 2016", MOUSE, c(-0.8, 1.7), c(0, 1), "Contribution (months per 100 CAG repeats)"))

fig <- fig_compose(P(a), P(b), P(c_), design = "A\nB\nC", heights = c(52, 56, 50))
save_fig(fig, "sfig10", 165)
