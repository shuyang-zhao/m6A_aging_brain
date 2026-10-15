# Additional file 1 Fig. S4: YTHDF3, a member of the reading axis
# reads source_data/sfig4a, sfig4b, sfig4b_controls, sfig4c and sfig4d; writes output/figures/sfig4.png and .pdf
source("figures/style.R")

ordinal <- function(x) {
  paste0(round(x), ifelse(round(x) %in% 11:13, "th", c("th", "st", "nd", "rd", rep("th", 6))[round(x) %% 10 + 1]))
}
stars <- function(q) as.character(cut(q, c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", "")))

# a: YTHDF3 age effect; mouse protein per cohort (late versus early age group) and combined (random effects over
# cohorts), human mRNA combined per context (per decade); bars, estimate; whiskers, 95% CI; * q < 0.05
HUMAN <- c(striatum = "Striatum", cortex = "Cortex", hippocampus = "Hippocampus", cerebellum = "Cerebellum",
           neurons = "Neurons")
d <- read_panel("sfig4a")
d[, `:=`(combined = label == "Combined" | panel == "Human mRNA",
         label = fifelse(panel == "Human mRNA", HUMAN[label], label))]
d[, `:=`(panel = factor(panel, c("Mouse protein", "Human mRNA")), label = factor(label, unique(label)),
         key = fcase(!combined, "Cohort", panel == "Mouse protein", "Protein", label == "Neurons", "Single-nucleus RNA",
                     default = "Bulk RNA"))]
d[, inside := fifelse(key == "Cohort", GREY, LAYER[key])]
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
        plot_layout(heights = c(8, 5))) + FLUSH

# b: pooled load edge of the YTHDF module and YTHDF3 against the control RBPs (density of one pooled edge per RBP;
# random effects over units)
PLEV <- c("Mouse Output", "Mouse Translation", "Mouse Stability", "Human Output", "Human Translation",
          "Human Stability")
panel_of <- function(d) factor(paste(c(mouse = "Mouse", human = "Human")[d$species], d$arm), PLEV)
ctrl <- read_panel("sfig4b_controls")
ctrl[, panel := panel_of(.SD)]
hit <- read_panel("sfig4b")
hit[, `:=`(panel = panel_of(.SD), who = c(YTHDF = "Module", YTHDF3 = "YTHDF3")[regulator])]
hit <- hit[!is.na(panel) & panel %in% ctrl$panel]
hit[, label := paste0(who, " ", ordinal(pctile_signed))]
names_p <- hit[order(who), .(strip = paste0(panel[1], "\n", paste(label, collapse = ", "))), by = panel]
strip_of <- function(p) factor(names_p$strip[match(p, names_p$panel)], names_p$strip[order(names_p$panel)])
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
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "top",
        legend.margin = margin(0, 0, -2, 0), legend.key.width = unit(6, "mm"),
        panel.spacing = unit(5, "mm"), strip.text = element_text(lineheight = 0.95))

# c: fraction of the age slope of the m6A load that the YTHDF terms take over; one dot per unit with its
# bootstrap 95% CI where at most two units are available (the intervals of single units in many-unit arms exceed
# 100%; a triangle marks an interval cut at the axis); bar, median over units
m <- read_panel("sfig4c")
m[, `:=`(set = factor(c(YTHDF = "YTHDF module", YTHDF3 = "YTHDF3")[regulator], c("YTHDF module", "YTHDF3")),
         panel = factor(paste("Mouse", arm), c("Mouse Output", "Mouse Translation", "Mouse Stability")),
         f = fraction_pct, lo = ci_low_pct, hi = ci_high_pct)]
m[, n_units := .N, by = .(panel, set)]
m[, x := as.numeric(set) + (seq_len(.N) - (.N + 1) / 2) / .N * 0.7, by = .(panel, set)]
med <- m[, .(f = median(f)), by = .(panel, set)]
LIM <- c(-30, 120)
ci <- m[n_units <= 2]
cc <- ggplot(m, aes(x, f)) +
  ZERO_Y +
  geom_segment(data = med, aes(x = as.numeric(set) - 0.4, xend = as.numeric(set) + 0.4, y = f, yend = f),
               linewidth = 1.4, inherit.aes = FALSE, colour = INK) +
  geom_linerange(data = ci, aes(ymin = pmax(lo, LIM[1]), ymax = pmin(hi, LIM[2])), linewidth = 0.6, colour = CI) +
  geom_point(data = ci[hi > LIM[2]], aes(y = LIM[2], fill = set), shape = 24, size = 2, stroke = 0.5, colour = EDGE,
             show.legend = FALSE) +
  geom_point(data = ci[lo < LIM[1]], aes(y = LIM[1], fill = set), shape = 25, size = 2, stroke = 0.5, colour = EDGE,
             show.legend = FALSE) +
  geom_point(aes(fill = set), shape = 21, size = 2.4, stroke = 0.5, colour = EDGE) +
  facet_grid(. ~ panel, scales = "free_x", space = "free_x", labeller = label_wrap_gen(8)) +
  scale_fill_manual(values = c("YTHDF module" = CONTRAST, YTHDF3 = FOCAL), name = NULL) +
  scale_x_continuous(breaks = NULL, expand = expansion(add = 0.5)) +
  coord_cartesian(ylim = LIM) +
  scale_y_continuous(labels = MINUS_LABELS) +
  labs(x = NULL, y = "m6A-linked age slope mediated (%)") +
  theme_fig() +
  BOXED +
  theme(legend.position = "bottom", axis.ticks.x = element_blank(), strip.text = element_text(lineheight = 0.95),
        panel.spacing = unit(4, "mm"))

# d: change of the target score after loss of a YTHDF reader (m6A-marked targets of the preferred binding source,
# against bound transcripts without m6A or unbound m6A genes); filled, the direction predicted from the reader's known
# function (loss lowers the translation efficiency of YTHDF1/3 targets, raises the mRNA of YTHDF2 targets); open, the
# opposite direction; the triple knockout has a row for each reader whose targets are scored
x <- read_panel("sfig4d")
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
x[, `:=`(label = LABELS[dataset],
         comparison = factor(c(bound_unmarked = "Bound without m6A", unbound_m6A = "Unbound m6A genes")[comparison],
                             c("Bound without m6A", "Unbound m6A genes")))]
x[reader != target_reader, label := paste0(label, "\n(", target_reader, " targets)")]
x[, label := factor(label, unique(label[order(panel, -effect_pct)]))]
x[, label := factor(label, rev(levels(label)))]
SERIES <- c("Bound without m6A" = BLUE, "Unbound m6A genes" = CONTRAST)
x[, `:=`(inside = fifelse(sign_matches, SERIES[as.character(comparison)], "white"),
         edge = fifelse(sign_matches, EDGE, GREY))]
dd <- ggplot(x, aes(effect_pct, label, shape = comparison, fill = I(inside), colour = I(edge))) +
  ZERO_X +
  geom_errorbar(aes(xmin = ci_low_pct, xmax = ci_high_pct), orientation = "y", width = 0, linewidth = 0.6,
                colour = CI, position = position_dodge(0.6)) +
  geom_point(size = 2.8, stroke = 0.5, position = position_dodge(0.6)) +
  facet_grid(panel ~ ., scales = "free", space = "free_y", labeller = label_wrap_gen(16)) +
  scale_shape_manual(values = c("Bound without m6A" = 21, "Unbound m6A genes" = 24), name = NULL,
                     guide = guide_legend(override.aes = list(fill = SERIES, colour = EDGE))) +
  scale_x_continuous(labels = MINUS_LABELS) +
  labs(x = "Change of target score after YTHDF loss (%)", y = NULL) +
  theme_rows() +
  BOXED +
  theme(legend.position = "bottom")

fig <- panel_grid(list(a = a, b = b, c = cc, d = dd), design = "AB\nCD", widths = c(76, 107), heights = c(74, 92))
save_fig(fig, "sfig4", 332)
