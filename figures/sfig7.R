# Fig. S7: YTHDF3 and the YTHDF readers
# reads source_data/sfig7a, sfig7b, sfig7b_controls and sfig7c; writes output/figures/sfig7
source("figures/style.R")

# a: YTHDF3 age effect (95% CI): mouse protein per cohort (late versus early age group) and combined (random effects
# over cohorts), human mRNA per context (per decade); bars sorted by effect, combined estimates in the colour of the
# species
CONTEXTS <- c(striatum = "Striatum", cortex = "Cortex", hippocampus = "Hippocampus", cerebellum = "Cerebellum",
              neurons = "Neurons")
d <- read_panel("sfig7a")
d[, `:=`(combined = label == "Combined" | panel == "Human mRNA", label = fifelse(panel == "Human mRNA", CONTEXTS[label], label))]
d[, fill := fcase(!combined, GREY, panel == "Mouse protein", MOUSE, default = HUMAN)]
d[, panel := factor(panel, c("Mouse protein", "Human mRNA"), c("Mouse protein,\nlate versus early", "Human mRNA,\nper decade"))]
d <- d[order(panel, label == "Combined", effect)]
d[, row := factor(paste(panel, label), rev(paste(panel, label)))]
a <- ggplot(d, aes(effect, row)) +
  geom_vline(xintercept = 0, linewidth = lw(0.5)) +
  geom_col(aes(fill = I(fill)), width = 0.68) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), linewidth = lw(0.5), colour = EDGE) +
  geom_text(data = d[!is.na(q)], aes(x = pmin(ci_low, 0), label = fmt_p(q, prefix = "italic(q)")), parse = TRUE,
            hjust = 1.08, size = pt_size(6)) +
  ggh4x::facet_grid2(panel ~ ., scales = "free", space = "free_y", independent = "x") +
  scale_x_continuous(labels = MINUS_LABELS, expand = expansion(mult = c(0.4, 0.04))) +
  scale_y_discrete(labels = setNames(d$label, d$row)) +
  labs(x = "YTHDF3 age effect (%)", y = NULL) +
  theme_box() +
  theme(strip.text.y = element_text(angle = -90, size = 6))

# b: load coupling of the YTHDF module and of YTHDF3 (ticks) against the density of the same coupling of the control
# proteins, in units of the 95th percentile of their absolute coupling in each unit; right, percentile of the module
# and of YTHDF3 among the controls
ctrl <- read_panel("sfig7b_controls")
hit <- read_panel("sfig7b")
ctrl[, row := layer_rows(species, layer)]
hit[, row := factor(paste(SPECIES_NAME[species], layer), levels(ctrl$row))]
ctrl[, unit := quantile(abs(edge), 0.95), by = row]
hit[ctrl[, .(unit = unit[1]), by = row], unit := i.unit, on = "row"]
PEAK <- 1.2
curves <- ctrl[, {
  k <- density(edge / unit[1], adjust = 0.8)
  .(x = k$x, h = PEAK * k$y / max(k$y))
}, by = row]
curves[, base := nlevels(row) - as.integer(row)]
hit[, `:=`(base = nlevels(row) - as.integer(row), x = estimate / unit, tick = c(YTHDF = 1.1, YTHDF3 = 0.6)[regulator],
           colour = SPECIES[SPECIES_NAME[species]])]
hit[regulator == "YTHDF3", colour := shade(colour, 0.4)]
setorder(hit, regulator)
percentile <- hit[, .(text = paste(ordinal(pctile_signed[regulator == "YTHDF"]), ordinal(pctile_signed[regulator == "YTHDF3"]),
                                   sep = ", "), base = base[1]), by = row]
ridge <- function(r) {
  geom_ribbon(data = curves[row == r], aes(x, ymin = base, ymax = base + h), fill = LIGHT, colour = GREY,
              linewidth = lw(0.5), outline.type = "upper")
}
XR <- range(curves$x)
KEY <- nlevels(ctrl$row) - 1 + PEAK + 0.3 # key row above the top ridge, ticks at half their length
b <- ggplot() +
  ref_line(xintercept = 0) +
  lapply(levels(curves$row), ridge) +
  geom_segment(data = hit, aes(x = x, xend = x, y = base, yend = base + tick, colour = I(colour)), linewidth = lw(1.4)) +
  geom_text(data = percentile, aes(x = XR[2] + 0.04 * diff(XR), y = base + 0.45, label = text), hjust = 0, size = pt_size(6)) +
  annotate("segment", x = XR[1] + c(0.05, 1.3), xend = XR[1] + c(0.05, 1.3), y = KEY, yend = KEY + c(0.55, 0.3),
           linewidth = lw(1.4)) +
  note(XR[1] + c(0.2, 1.45), KEY + 0.25, c("Module", "YTHDF3"), hjust = 0) +
  note(XR[2] + 0.04 * diff(XR), KEY + c(0.45, 0.05), c("Percentile:", "module, YTHDF3"), hjust = 0) +
  scale_x_continuous(labels = MINUS_LABELS, expand = c(0, 0)) +
  scale_y_continuous(breaks = nlevels(ctrl$row) - seq_len(nlevels(ctrl$row)), labels = levels(ctrl$row),
                     expand = expansion(add = c(0.1, 0))) +
  coord_cartesian(clip = "off", xlim = XR, ylim = c(0, KEY + 0.7)) +
  labs(x = "YTHDF load coupling / control 95th percentile", y = NULL) +
  theme_pub() +
  theme(plot.margin = margin(1.5, 17, 1.5, 1.5, "mm"))

# c: change of the target score after loss of a YTHDF reader (%; m6A-marked targets of the preferred binding source
# against bound transcripts without m6A or unbound m6A genes); bold, the change agrees with the reader's known function
# (loss lowers the translation efficiency of YTHDF1/3 targets, raises the mRNA of YTHDF2 targets)
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
BLOCKS <- c("YTHDF1/3 loss:\ntranslation efficiency", "YTHDF2 loss:\ntarget mRNA")
x <- read_panel("sfig7c")
x[, `:=`(panel = factor(fifelse(readout == "RNA", BLOCKS[2], BLOCKS[1]), BLOCKS),
         comparison = factor(c(bound_unmarked = "Bound, no m6A", unbound_m6A = "Unbound, m6A")[comparison],
                             c("Bound, no m6A", "Unbound, m6A")),
         label = LABELS[dataset])]
x[reader != target_reader, label := paste0(label, " (", target_reader, " targets)")]
x[, `:=`(label = factor(label, unique(label[order(panel, -effect_pct)])), text = number(effect_pct, 1))]
block <- function(name) {
  d <- x[panel == name][, label := droplevels(label)]
  ggplot(d, aes(comparison, label, fill = effect_pct)) +
    geom_tile(colour = "white", linewidth = lw(0.6)) +
    geom_text(aes(label = text, fontface = fifelse(sign_matches, "bold", "plain")), size = pt_size(6)) +
    facet_wrap(~panel) +
    scale_fill_signed(20, "Change of the target score (%)", light = TRUE) +
    scale_x_discrete(position = "top", expand = c(0, 0), drop = FALSE) +
    scale_y_discrete(limits = rev, labels = scales::label_wrap(30), expand = c(0, 0)) +
    theme_tile() +
    theme(strip.placement = "outside", strip.text = element_text(size = 6),
          axis.text.x.top = element_text(vjust = 0, margin = margin(t = 1.5, b = 0.5, unit = "mm")))
}
rows <- x[, uniqueN(label), keyby = panel]$V1
cc <- block(BLOCKS[1]) / block(BLOCKS[2]) + plot_layout(heights = rows, guides = "collect") &
  theme(legend.position = "bottom")

fig <- fig_compose(a, b, cc, design = "AC\nBC", widths = c(88, 95), heights = c(1, 1))
save_fig(fig, "sfig7", 142)
