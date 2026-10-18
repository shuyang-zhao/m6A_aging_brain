# Additional file 1 Fig. S7: YTHDF3 and the YTHDF readers
# reads source_data/sfig7a, sfig7b, sfig7b_controls and sfig7c; writes output/figures/sfig7.png and .pdf
source("figures/style.R")

stars <- function(q) as.character(cut(q, c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", "")))

# a: YTHDF3 age effect; mouse protein per cohort (late versus early age group) and combined (random effects over
# cohorts), human mRNA combined per context (per decade); bars, estimate; whiskers, 95% CI; * q < 0.05
HUMAN <- c(striatum = "Striatum", cortex = "Cortex", hippocampus = "Hippocampus", neurons = "Neurons")
d <- read_panel("sfig7a")
d[, `:=`(combined = label == "Combined" | panel == "Human mRNA",
         label = fifelse(panel == "Human mRNA", HUMAN[label], label))]
d[, `:=`(panel = factor(panel, c("Mouse protein", "Human mRNA")), label = factor(label, unique(label)),
         inside = fcase(!combined, GREY, panel == "Mouse protein", FOCAL, default = RED))]
bars <- function(dd, xlab) {
  ggplot(dd, aes(effect, label, fill = I(inside))) +
    ZERO_X +
    geom_col(width = 0.68, colour = EDGE, linewidth = 0.5) +
    geom_linerange(aes(xmin = ci_low, xmax = ci_high), linewidth = 0.6, colour = EDGE) +
    geom_text(data = dd[q < 0.05], aes(x = pmin(ci_low, 0), label = stars(q)), hjust = 1.2, vjust = 0.75, size = PT,
              colour = INK) +
    facet_wrap(~panel) +
    scale_x_continuous(labels = MINUS_LABELS, expand = expansion(mult = c(0.12, 0.05))) +
    scale_y_discrete(limits = rev) +
    labs(x = xlab, y = NULL) +
    theme_rows() +
    BOXED
}
a <- (bars(d[panel == "Mouse protein"], "YTHDF3 age effect (%, late versus early)") /
        bars(d[panel == "Human mRNA"], "YTHDF3 age effect (% per decade)") +
        plot_layout(heights = d[, .N, by = panel]$N + 1)) + FLUSH

# b: pooled load edge of the YTHDF module and YTHDF3 against the control RBPs (density of one pooled edge per RBP;
# random effects over units); numbers in the strip, percentile of the pooled edge among the controls
ctrl <- read_panel("sfig7b_controls")
ctrl[, panel := layer_rows(species, layer)]
hit <- read_panel("sfig7b")
hit[, `:=`(panel = layer_rows(species, layer), who = c(YTHDF = "Module", YTHDF3 = "YTHDF3")[regulator])]
hit[, label := paste0(who, " ", ordinal(pctile_signed))]
strips <- hit[order(who), .(strip = paste0(panel[1], "\n", paste(label, collapse = ", "))), by = panel][order(panel)]
strip_of <- function(p) factor(strips$strip[match(p, strips$panel)], strips$strip)
ctrl[, panel := strip_of(panel)]
hit[, panel := strip_of(panel)]
b <- ggplot(ctrl, aes(edge)) +
  ZERO_X +
  geom_density(fill = LIGHT, colour = GREY, linewidth = 0.5, adjust = 0.8) +
  geom_vline(data = hit, aes(xintercept = estimate, colour = who), linewidth = 1.2) +
  facet_wrap(~panel, nrow = 2, scales = "free") +
  scale_colour_manual(values = c(Module = CONTRAST, YTHDF3 = FOCAL),
                      labels = c(Module = "YTHDF module", YTHDF3 = "YTHDF3"), name = NULL) +
  scale_x_continuous(labels = MINUS_LABELS, n.breaks = 3, expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "YTHDF load edge (b)", y = "Density of control RBPs") +
  theme_fig() +
  BOXED +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "inside",
        legend.position.inside = c(0.85, 0.2), legend.key.width = unit(6, "mm"), panel.spacing = unit(5, "mm"),
        strip.text = element_text(lineheight = 0.95))

# c: change of the target score after loss of a YTHDF reader (m6A-marked targets of the preferred binding source, against
# bound transcripts without m6A or unbound m6A genes); filled, the direction predicted from the reader's known function
# (loss lowers the translation efficiency of YTHDF1/3 targets, raises the mRNA of YTHDF2 targets); open, the opposite
# direction; the triple knockout has a row for each reader whose targets are scored
x <- read_panel("sfig7c")
x[, panel := factor(fifelse(readout == "RNA", "YTHDF2 loss: target mRNA", "YTHDF1/3 loss: translation efficiency"),
                    c("YTHDF1/3 loss: translation efficiency", "YTHDF2 loss: target mRNA"))]
LABELS <- c(
  gse266096_nsc_Ythdf2_cKO = "Zhang 2025, neural stem cells, Ythdf2 cKO",
  heck2020_iPSC_siYTHDF2 = "Heck 2020, iPSC-derived cells, siYTHDF2",
  lasman2020_mESC_Ythdf1_KO = "Lasman 2020, mESC, Ythdf1 KO",
  lasman2020_mESC_Ythdf2_KO = "Lasman 2020, mESC, Ythdf2 KO",
  lasman2020_mESC_Ythdf3_KO = "Lasman 2020, mESC, Ythdf3 KO",
  lasman2020_mESC_Ythdf321_KO = "Lasman 2020, mESC, Ythdf1/2/3 triple KO",
  li2018_nspc_Ythdf2_cKO = "Li 2018, neural progenitors, Ythdf2 cKO",
  shi2017_HeLa_siYTHDF3 = "Shi 2017, HeLa, siYTHDF3",
  wang2014_HeLa_siYTHDF2 = "Wang 2014, HeLa, siYTHDF2",
  wang2015_HeLa_siYTHDF1 = "Wang 2015, HeLa, siYTHDF1",
  yu2021_cgc_shYthdf2 = "Yu 2021, cerebellar granule cells, shYthdf2",
  zhuang2023_dg_Ythdf2_cKO = "Zhuang 2023, dentate gyrus, Ythdf2 cKO",
  zou2023_HeLa_siDF2 = "Zou 2023, HeLa, siYTHDF2",
  zou2023_neuron_ko_YTHDF1inh = "Zou 2023, neurons, YTHDF1 inhibitor (KO background)",
  zou2023_neuron_wt_YTHDF1inh = "Zou 2023, neurons, YTHDF1 inhibitor (wild type)"
)
SERIES <- c("Bound without m6A" = BLUE, "Unbound m6A genes" = CONTRAST)
x[, `:=`(label = LABELS[dataset],
         comparison = factor(c(bound_unmarked = "Bound without m6A", unbound_m6A = "Unbound m6A genes")[comparison],
                             names(SERIES)))]
x[reader != target_reader, label := paste0(label, " (", target_reader, " targets)")]
x[, label := factor(label, unique(label[order(panel, -effect_pct)]))]
x[, `:=`(inside = fifelse(sign_matches, SERIES[as.character(comparison)], "white"),
         edge = fifelse(sign_matches, EDGE, GREY))]
dodge <- position_dodge(0.6)
cc <- ggplot(x, aes(effect_pct, label, shape = comparison, fill = I(inside), colour = I(edge))) +
  ZERO_X +
  geom_linerange(aes(xmin = ci_low_pct, xmax = ci_high_pct), orientation = "y", linewidth = 0.6, colour = CI,
                 position = dodge) +
  geom_point(size = 2.8, stroke = 0.5, position = dodge) +
  facet_wrap(~panel, scales = "free") +
  scale_shape_manual(values = c("Bound without m6A" = 21, "Unbound m6A genes" = 24), name = NULL,
                     guide = guide_legend(override.aes = list(fill = SERIES, colour = EDGE))) +
  scale_x_continuous(labels = MINUS_LABELS) +
  scale_y_discrete(limits = rev, labels = scales::label_wrap(36)) +
  labs(x = "Change of target score after YTHDF loss (%)", y = NULL) +
  theme_rows() +
  BOXED +
  theme(legend.position = "bottom")

fig <- panel_grid(list(a = a, b = b, c = cc), design = "AB\nCC", widths = c(70, 113), heights = c(75, 70))
save_fig(fig, "sfig7", 290)
