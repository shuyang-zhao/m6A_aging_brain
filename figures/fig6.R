# Figure 6: the proteome clock in neurodegeneration
# reads source_data/fig6a-e (script 53); writes output/figures/fig6.png and .pdf
source("figures/style.R")

tests <- read_panel("fig6c")
CLOCK_Y <- "Clock score"

# a: AD, donor scores of cases (red) and controls (grey) in two studies (box, median and quartiles); above, the
# coefficient of AD on the score (SD; adjusted for age, sex and batch) with its P
a <- read_panel("fig6a")
note <- tests[disease == "AD" & study %in% a$study, .(strip = factor(study, unique(a$study)),
                                                      text = paste0(signed(est), " SD\n", fmt_p(p)))]
a[, `:=`(group = factor(group, c("Control", "AD")), strip = factor(study, unique(study)))]
a <- ggplot(a, aes(group, score)) +
  geom_boxplot(aes(colour = group), fill = NA, width = 0.5, linewidth = 0.6, outliers = FALSE, show.legend = FALSE) +
  geom_point(aes(fill = group), colour = EDGE, shape = 21, size = 2.8, stroke = 0.4, alpha = 0.9,
             position = position_jitter(width = 0.14, height = 0, seed = 1)) +
  geom_text(data = note, aes(x = -Inf, y = Inf, label = text), hjust = -0.08, vjust = 1.15, size = PT, lineheight = 0.95,
            colour = INK, inherit.aes = FALSE) +
  facet_wrap(~strip) +
  scale_fill_manual(values = c(Control = GREY, AD = RED), guide = "none") +
  scale_colour_manual(values = c(Control = EDGE, AD = EDGE)) +
  scale_y_continuous(labels = MINUS_LABELS, expand = expansion(mult = c(0.05, 0.3))) +
  labs(x = NULL, y = CLOCK_Y) +
  theme_fig() +
  BOXED

# b: HD allelic series at 6 months, score of each mouse against the CAG length of its genotype (red, expanded alleles;
# grey, 7 and 20 repeats, below the 36 of the HD range); line, least squares
b <- read_panel("fig6b")
b[, expanded := cag >= 36]
b <- ggplot(b, aes(cag, score)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = EDGE, linewidth = 0.8) +
  geom_point(aes(fill = expanded), colour = EDGE, shape = 21, size = 2.8, stroke = 0.4, alpha = 0.9,
             position = position_jitter(width = 4, height = 0, seed = 1)) +
  scale_fill_manual(values = c("TRUE" = RED, "FALSE" = GREY), guide = "none") +
  scale_x_continuous(breaks = c(0, 100, 200), expand = expansion(add = 12)) +
  scale_y_continuous(labels = MINUS_LABELS) +
  labs(x = sprintf("CAG repeats of the genotype (%d-month-old mice)", b$age[1]), y = CLOCK_Y) +
  theme_fig() +
  BOXED

# c and e: coefficient of the disease (SD; point and 95% CI; filled, P < 0.05) for each study, AD studies combined
# (diamond), HD per 100 CAG repeats; P on the right; c on the clock score, e on the individual m6A index
STUDY_ROWS <- c("Ping 2020", "Merrihew 2023", "AD combined", "Blumenreich 2024", "Umoh 2018", "Langfelder 2016")
forest <- function(d, key, x) {
  d <- copy(d)
  d[, label := fifelse(study == "AD combined", "Combined", study)]
  d[, `:=`(row = factor(label, rev(label[match(STUDY_ROWS, study)])), sig = p < 0.05,
           disease = factor(disease, c("AD", "PD", "ALS", "HD"), c("AD", "PD", "ALS", "HD, per\n100 CAG")))]
  ggplot(d, aes(est, row)) +
    ZERO_X +
    geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = CI, linewidth = 0.6) +
    geom_point(data = ~ .x[study != "AD combined"], aes(fill = sig), colour = EDGE, shape = 21, size = 3.8, stroke = 0.5) +
    geom_point(data = ~ .x[study == "AD combined"], aes(fill = sig), colour = EDGE, shape = 23, size = 5.2, stroke = 0.5) +
    geom_text(aes(x = Inf, label = fmt_p(p)), hjust = -0.12, size = PT, colour = INK) +
    facet_grid(disease ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = c("TRUE" = key, "FALSE" = "white"), guide = "none") +
    scale_x_continuous(labels = MINUS_LABELS, breaks = scales::breaks_pretty(4)) +
    coord_cartesian(clip = "off") +
    labs(x = x, y = NULL) +
    theme_rows() +
    BOXED +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0), panel.spacing.y = unit(3, "mm"),
          plot.margin = margin(10.4, 28, 3, 2, "mm"))
}
cc <- forest(tests, RED, "Coefficient on the clock score (SD)")

# d: the clock is a weighted sum over genes, so the disease coefficient splits exactly into the coefficients of the
# partial scores of groups of genes; contribution of each group (bar, 95% CI) to the coefficient of AD (studies
# combined) and of HD (per 100 CAG repeats), by cell source and by m6A sites; tick, the contribution expected from the
# group's share of the clock weights (total x weight share); right, share of the total
GROUPS <- c(neuron = "Neuron", astrocyte = "Astrocyte", oligodendrocyte = "Oligodendrocyte", microglia = "Microglia",
            ubiquitous = "Ubiquitous", "no cell-type data" = "No cell-type data",
            "0" = "0 sites", "1" = "1 site", "2-3" = "2–3 sites", "≥4" = "≥4 sites")
PARTS <- c("cell source" = "Cell source", "m6A sites" = "m6A sites")
dd <- read_panel("fig6d")
heads <- dd[partition == "total", setNames(fifelse(disease == "HD", sprintf("HD: %s SD per 100 CAG", signed(est)),
                                                   sprintf("AD combined: %s SD", signed(est))), study)]
dd <- dd[partition != "total"]
dd[, `:=`(row = factor(GROUPS[group], rev(GROUPS)), part = factor(PARTS[partition], PARTS),
          head = factor(heads[study], heads[c("AD combined", "Langfelder 2016")]))]
d <- ggplot(dd, aes(est, row)) +
  geom_col(aes(fill = part), width = 0.66, colour = EDGE, linewidth = 0.5, show.legend = FALSE) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = EDGE, linewidth = 0.6) +
  geom_errorbar(aes(xmin = expected, xmax = expected, linetype = "Expected from the share of the clock weights"),
                width = 0.8, colour = INK, linewidth = 1, key_glyph = "vline") +
  geom_text(aes(x = Inf, label = sprintf("%.0f%%", 100 * share_of_total)), hjust = 1.15, size = PT, colour = INK) +
  facet_grid(part ~ head, scales = "free", space = "free_y") +
  scale_fill_manual(values = c("Cell source" = BLUE, "m6A sites" = FOCAL)) +
  scale_linetype_manual(values = "solid", name = NULL) +
  scale_x_continuous(labels = MINUS_LABELS, breaks = scales::breaks_pretty(4), expand = expansion(mult = c(0.03, 0.3))) +
  labs(x = "Contribution to the disease coefficient (SD of the clock score)", y = NULL) +
  theme_rows() +
  BOXED +
  theme(panel.spacing.y = unit(4, "mm"), strip.text.y = element_text(angle = -90, margin = margin(0, 0, 0, 2, "mm")),
        legend.position = "top", legend.justification = "right", legend.margin = margin(0, 0, -2, 0))

e <- forest(read_panel("fig6e"), FOCAL, "Coefficient on the m6A index (SD)")

fig <- panel_grid(list(a = a, b = b, c = cc, d = d, e = e), design = "ABC\nDDE", widths = c(100, 105, 161),
                  heights = c(100, 115))
save_fig(fig, "fig6", 215)
