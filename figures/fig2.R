# Figure 2: m6A load marks mouse proteins that decline with age
# reads source_data/fig2a-e and fig2c_load; writes output/figures/fig2.png and .pdf
source("figures/style.R")
suppressPackageStartupMessages(library(ggbeeswarm))

LAYERS <- c(Protein = "Protein", "mRNA-adjusted protein" = "Protein\ngiven mRNA", mRNA = "mRNA")

# a: signature z of the m6A-marked genes in every mouse unit (cameraPR against unmodified genes), by layer (colour of
# the layer; neurons, single-nucleus RNA); filled, P < 0.05; open, not significant; bar, median over units
d <- read_panel("fig2a")
d[, `:=`(key = fcase(context == "mouse_neurons", "Single-nucleus RNA", layer == "Protein", "Protein",
                     layer == "mRNA-adjusted protein", "Protein given mRNA", default = "Bulk RNA"),
         layer = factor(LAYERS[layer], LAYERS))]
d[, c("edge", "inside") := p_colours(p, key)]
a <- ggplot(d, aes(layer, z)) +
  ZERO_Y +
  stat_summary(aes(colour = layer), fun = median, geom = "crossbar", width = 0.6, linewidth = 0.7) +
  geom_quasirandom(aes(colour = I(edge), fill = I(inside)), shape = 21, size = 3.6, stroke = 0.5, width = 0.28) +
  scale_colour_manual(values = unname(LAYER[c("Protein", "Protein given mRNA", "Bulk RNA")]), guide = "none") +
  scale_y_continuous(labels = MINUS_LABELS) +
  labs(x = NULL, y = "m6A signature (cameraPR z)") +
  theme_fig() +
  theme(axis.ticks.x = element_blank())

# b: age effect of mouse genes by number of m6A sites: protein, mean over units with a 95% interval across units
# (number, the load coefficient over cohorts as cited in the text); protein given mRNA, one line per unit with its 95%
# CI
b <- sites_panel("fig2b", "Mouse", axes = labs(x = "m6A sites per gene", y = "Age effect (%)"))

# c: age effect of every protein by the number of m6A sites of its gene in the largest mouse cohort (Wang 2025, 40
# brains); violins and boxes over genes, the central 98% of effects shown; number, the load coefficient of the unit
genes <- read_panel("fig2c")
genes[, sites := cut(n_sites_m6A, c(-Inf, 0, 1, 3, Inf), c("0", "1", "2–3", "≥4"))]
counts <- genes[, .N, by = sites]
unit_coef <- read_panel("fig2c_load")
cc <- ggplot(genes, aes(sites, effect_pct)) +
  ZERO_Y +
  geom_violin(aes(fill = sites), colour = NA, scale = "width", width = 0.85, adjust = 1.2) +
  geom_boxplot(width = 0.16, outliers = FALSE, fill = "white", colour = EDGE, linewidth = 0.6) +
  geom_text(data = counts, aes(y = -Inf, label = format(N, big.mark = ",")), vjust = -0.5, size = PT, colour = INK) +
  annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1.2, size = PT, colour = INK, lineheight = 0.95,
           label = paste0(signed(unit_coef$effect), "% per log unit\n", fmt_p(unit_coef$p))) +
  scale_fill_manual(values = c("#E3E3E3", "#D7E9EE", FOCAL_LIGHT, "#5EA8BC"), guide = "none") +
  scale_y_continuous(labels = MINUS_LABELS) +
  coord_cartesian(ylim = genes[, quantile(effect_pct, c(0.01, 0.99))]) +
  labs(x = "m6A sites per gene", y = "Age effect, Wang 2025 (%)") +
  theme_fig() +
  theme(axis.ticks.x = element_blank())

# d: age effect of the m6A load with the covariates added in turn (OLS with HC3 standard errors, cohorts combined by
# random effects; 95% CI); filled, P < 0.05; number, the coefficient of the main model as a percentage of the
# coefficient without covariates
STEPS <- c(none = "No covariates", lengths = "+ Four lengths", expression = "+ Expression",
           "neuron share" = "+ Neuron share", programmes = "+ Programmes (main model)",
           "other loads" = "+ Ψ, m5C, m7G loads")
d <- read_panel("fig2d")
d[, `:=`(panel = factor(panel, CONTEXTS), step = factor(STEPS[step], rev(STEPS)))]
d[, key := KEYS[as.character(panel)]]
d[, c("edge", "inside") := p_colours(p, key)]
main <- d[step == STEPS[["programmes"]]]
dd <- ggplot(d, aes(effect, step)) +
  ZERO_X +
  geom_path(aes(group = panel, colour = I(LAYER_LIGHT[key])), linewidth = 0.9) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = CI, linewidth = 0.6) +
  geom_point(aes(colour = I(edge), fill = I(inside)), shape = 21, size = 3.8, stroke = 0.5) +
  geom_text(data = main, aes(x = Inf, label = paste0(round(retained_pct), "%")), hjust = 1.15, size = PT,
            colour = INK) +
  facet_wrap(~panel, nrow = 1, scales = "free_x", labeller = label_wrap_gen(14)) +
  scale_x_continuous(labels = MINUS_LABELS, breaks = \(l) Filter(\(x) x <= 0, scales::breaks_pretty(3)(l)),
                     expand = expansion(mult = c(0.06, 0.42))) +
  labs(x = "Age effect of the m6A load per log unit (%)", y = NULL) +
  theme_rows() +
  BOXED +
  theme(panel.spacing = unit(5, "mm"), strip.text = element_text(lineheight = 0.95))

# e: average partial effect of the m6A load on the protein age effect among the long-lived proteins and among the other
# proteins with a measured lifetime (causal forests, mouse protein; random effects over cohorts; 95% CI); q, of the
# difference
h <- read_panel("fig2e")
inout <- h[quantity %in% c("in", "out")]
inout[, row := factor(quantity, c("out", "in"), c("Other proteins\nwith a lifetime", "Long-lived\nproteins"))]
e <- ggplot(inout, aes(effect, row, colour = quantity, fill = quantity)) +
  ZERO_X +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = CI, linewidth = 0.6) +
  geom_point(shape = 21, size = 4.4, stroke = 0.5) +
  annotate("text", x = -Inf, y = 2.45, hjust = 0, size = PT, colour = INK,
           label = fmt_p(h[quantity == "difference", q], "q")) +
  scale_colour_manual(values = c("in" = EDGE, "out" = GREY), guide = "none") +
  scale_fill_manual(values = c("in" = FOCAL, "out" = "white"), guide = "none") +
  scale_x_continuous(labels = MINUS_LABELS) +
  labs(x = "Age effect of the m6A load per log unit (%)", y = NULL) +
  theme_rows()

fig <- panel_grid(list(a = a, b = b, c = cc, d = dd, e = e), design = "ABC\nDDE", widths = c(52, 76, 55),
                  heights = c(86, 76))
save_fig(fig, "fig2", 324)
