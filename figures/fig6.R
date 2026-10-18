# Figure 6: the reading layer
# reads source_data/fig6a, fig6b, fig6b_controls, fig6c_fmr1 and fig6c_ythdf (script 54 and the edges of script 33);
# writes output/figures/fig6.png and .pdf
source("figures/style.R")

species_strip <- labeller(species = function(x) SPECIES_NAME[x])
RIGHT_STRIP <- theme(panel.spacing.y = unit(4, "mm"), strip.text.y = element_text(angle = -90))

# a: Spearman rho of the individual m6A load index with age in every cohort (point; filled where the Fisher-z 95%
# interval excludes 0) and the random-effects estimate over the cohorts of the species (diamond; I2)
a <- read_panel("fig6a")
a[, `:=`(species = factor(species, names(SPECIES_COL)), pooled = kind == "pooled")]
setorder(a, species, pooled, cohort)
a[, label := fifelse(pooled, "Pooled", sprintf("%s, n = %d", cohort, n))]
a[, label := factor(label, rev(c(unique(label[!pooled]), "Pooled")))]
a[, `:=`(sig = pooled | ci_low > 0 | ci_high < 0)]
a[, `:=`(fill = fifelse(sig, SPECIES_COL[as.character(species)], "white"), edge = fifelse(sig, EDGE, GREY))]
pooled_text <- a[pooled == TRUE, .(species, label, x = fifelse(species == "mouse", -0.3, 0.15),
                                   text = MINUS(sprintf("ρ = %.2f, I² = %.0f%%", rho, i2)))]
a <- ggplot(a, aes(rho, label)) +
  ZERO_X +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = CI, linewidth = 0.6) +
  geom_point(data = ~ .x[pooled == FALSE], aes(fill = I(fill), colour = I(edge)), shape = 21, size = 3.6, stroke = 0.5) +
  geom_point(data = ~ .x[pooled == TRUE], aes(fill = I(fill)), colour = EDGE, shape = 23, size = 5.2, stroke = 0.5) +
  geom_label(data = pooled_text, aes(x, label, label = text), hjust = 0, size = PT, colour = INK, fill = "white",
             border.colour = NA, label.padding = unit(0.6, "mm"), inherit.aes = FALSE) +
  facet_grid(species ~ ., scales = "free_y", space = "free_y", labeller = species_strip) +
  scale_x_continuous(labels = MINUS_LABELS, breaks = seq(-1, 1, 0.5), limits = c(-1, 1.25)) +
  labs(x = "Spearman ρ of the m6A load index with age", y = NULL) +
  theme_rows() +
  BOXED +
  RIGHT_STRIP

# b: coefficient of the reader level on the m6A load index (SD of the index per SD of the reader; random effects over
# the units of a species; 95% CI; filled, q < 0.05) against the same coefficient of the control RBPs (violin, dots);
# on the right, the percentile of the reader among the controls
ROWS_B <- c("FMR1", "FMRP family", "YTHDF module", "YTHDF3", "Control RBPs")
rd <- read_panel("fig6b")
ct <- read_panel("fig6b_controls")
rd[, `:=`(species = factor(species, names(SPECIES_COL)), row = factor(regulator, rev(ROWS_B)), sig = q < 0.05)]
rd[, `:=`(fill = fifelse(sig, SPECIES_COL[as.character(species)], "white"), edge = fifelse(sig, EDGE, GREY))]
ct[, `:=`(species = factor(species, names(SPECIES_COL)), row = factor("Control RBPs", rev(ROWS_B)))]
txt <- rbind(rd[, .(species, row, text = paste(ordinal(control_pct), "percentile"))],
             rd[, .(row = factor("Control RBPs", rev(ROWS_B)), text = paste0("n = ", n_controls[1])), by = species])
b <- ggplot(rd, aes(est, row)) +
  ZERO_X +
  geom_violin(data = ct, aes(x = est), orientation = "y", fill = LIGHT, colour = GREY, linewidth = 0.5,
              scale = "width", width = 0.8) +
  geom_point(data = ct, aes(x = est), colour = GREY, size = 1.4, alpha = 0.9,
             position = position_jitter(height = 0.18, width = 0, seed = 1)) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = CI, linewidth = 0.6) +
  geom_point(aes(fill = I(fill), colour = I(edge)), shape = 21, size = 3.8, stroke = 0.5) +
  geom_text(data = txt, aes(x = Inf, y = row, label = text), hjust = 1.04, size = PT, colour = INK, inherit.aes = FALSE) +
  facet_grid(species ~ ., labeller = species_strip) +
  scale_x_continuous(labels = MINUS_LABELS, limits = c(-0.045, 0.085), breaks = seq(-0.04, 0.04, 0.02)) +
  labs(x = "Change of the m6A load index per SD of the reader level (SD)", y = NULL) +
  theme_rows() +
  BOXED +
  RIGHT_STRIP

# c, left: edge of FMR1 to its m6A-marked targets (b_s; circles, filled where |b_s| exceeds the 95th percentile of the
# control RBPs of the unit) and to its bound targets without m6A (b_u; grey dots) in every unit, by layer, for the brain
# CLIP lists of Darnell 2011 and Maurin 2018 (human units through orthologs)
d <- read_panel("fig6c_fmr1")
d[, `:=`(row = layer_rows(species, layer), above = abs(b_s) > control_p95, key = SPECIES_COL[species])]
d[, fill := fifelse(above, key, "white")]
jitter <- position_jitter(height = 0.22, width = 0, seed = 2)
c1 <- ggplot(d, aes(y = row)) +
  ZERO_X +
  geom_point(aes(x = b_u), colour = GREY, size = 1.4, alpha = 0.8, position = jitter) +
  geom_point(aes(x = b_s, fill = I(fill), colour = I(fifelse(above, EDGE, GREY))), shape = 21, size = 3.2, stroke = 0.5,
             position = jitter) +
  facet_wrap(~list, ncol = 1) +
  scale_x_continuous(labels = MINUS_LABELS, breaks = c(-0.2, 0, 0.2)) +
  scale_y_discrete(limits = rev) +
  labs(x = "FMR1 edge (b); grey, targets without m6A", y = NULL) +
  theme_rows() +
  BOXED

# c, right: percentile of the load edge of the YTHDF module and of YTHDF3 among the load edges of the control RBPs
# (random effects over units); dashed line, 95th percentile
p <- read_panel("fig6c_ythdf")
p[, `:=`(set = factor(c(YTHDF = "YTHDF module", YTHDF3 = "YTHDF3")[regulator], c("YTHDF module", "YTHDF3")),
         row = layer_rows(species, layer), key = SPECIES_COL[species])]
c2 <- ggplot(p, aes(pctile_signed, row)) +
  geom_vline(xintercept = 95, linewidth = 0.4, colour = INK, linetype = "22") +
  geom_line(aes(group = row), colour = EDGE, linewidth = 0.6) +
  geom_point(aes(fill = I(key), shape = set), colour = EDGE, size = 3.8, stroke = 0.5) +
  scale_shape_manual(values = c("YTHDF module" = 21, YTHDF3 = 22), name = NULL,
                     guide = guide_legend(override.aes = list(fill = GREY))) +
  scale_x_continuous(breaks = c(0, 50, 95), expand = expansion(add = c(2, 0))) +
  scale_y_discrete(limits = rev) +
  coord_cartesian(xlim = c(0, 100), clip = "off") +
  labs(x = "Percentile among\ncontrol RBPs", y = NULL) +
  theme_rows() +
  theme(legend.position = "top", legend.justification = "right", legend.margin = margin(0, 0, -2, 0))
cc <- (c1 | c2) + plot_layout(widths = c(1, 1.15)) + FLUSH

fig <- panel_grid(list(a = a, b = b, c = cc), design = "AB\nAC", widths = c(150, 216), heights = c(95, 140))
save_fig(fig, "fig6", 235)
