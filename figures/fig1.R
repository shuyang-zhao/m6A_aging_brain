# Figure 1: resource, design and concordance of protein ageing between species
# reads source_data/fig1a.tsv, fig1c.tsv and fig1d.tsv; writes output/figures/fig1.png and .pdf
source("figures/style.R")

# a: age range (segment), donors of the largest unit (dot size at the midpoint of the range) and analysis units
# (right) of every dataset by species and molecular layer, rows from the top
ROW <- 5 # mm per row
PAD <- 0.7 # empty rows above and below the plot
GROUPS <- c("Human protein", "Human bulk RNA", "Human single-nucleus RNA", "Mouse protein", "Mouse RNA")
d <- read_panel("fig1a")
d[, group := factor(group, GROUPS)]
setorder(d, group, age_min, age_max)
d[, k := seq_len(.N), by = group]
offset <- d[, .(n = .N), by = group][, start := cumsum(c(0, head(n + 1.5, -1)))]
d <- merge(d, offset, by = "group")[, y := start + k]
heads <- offset[, .(group, y = start)]
AGE <- list(
  Human = list(lim = c(14, 107), breaks = seq(20, 100, 20), name = "Age (years)"),
  Mouse = list(lim = c(0, 35), breaks = seq(0, 30, 10), name = "Age (months)")
)
overview <- function(species) {
  age <- AGE[[species]]
  rows <- d[group %in% GROUPS[grepl(species, GROUPS)]]
  heads <- heads[group %in% rows$group]
  ggplot(rows, aes(colour = layer)) +
    geom_segment(aes(x = age_min, xend = age_max, y = y, yend = y), linewidth = 0.8, lineend = "round") +
    geom_point(aes((age_min + age_max) / 2, y, size = donors, fill = layer), shape = 21, colour = EDGE, stroke = 0.5) +
    geom_text(
      data = heads, aes(age$lim[1], y, label = group), inherit.aes = FALSE, hjust = 0, fontface = "bold", size = PT
    ) +
    annotate("text", x = age$lim[2], y = heads$y[1], label = "Units", hjust = -0.15, size = PT) +
    scale_y_reverse(
      breaks = rows$y, labels = rows$dataset, expand = expansion(add = PAD),
      sec.axis = dup_axis(breaks = rows$y, labels = rows$n_units, name = NULL)
    ) +
    scale_x_continuous(breaks = age$breaks, limits = age$lim, expand = expansion(0)) +
    scale_colour_manual(values = LAYER, guide = "none") +
    scale_fill_manual(values = LAYER, guide = "none") +
    scale_size(
      range = c(1.2, 5), limits = c(0, 300), breaks = c(10, 50, 100, 300), name = "Donors in the largest unit",
      guide = guide_legend(
        title.position = "top", nrow = 1, label.position = "bottom", override.aes = list(fill = GREY, colour = EDGE)
      )
    ) +
    coord_cartesian(clip = "off") +
    labs(x = age$name, y = NULL) +
    theme_fig() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(4.4, 11, 3, 2, "mm"))
}
# the panels are as high as their rows, so that a row is as tall in both plots; the donor legend takes the space
# under the mouse rows
row_span <- function(species) {
  g <- GROUPS[grepl(species, GROUPS)]
  d[group %in% g, max(y)] - heads[group %in% g, min(y)] + 2 * PAD
}
mouse <- overview("Mouse") / guide_area() +
  plot_layout(heights = unit(c(row_span("Mouse") * ROW, 1), c("mm", "null")), guides = "collect") +
  plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0), legend.justification = "top"))
a <- (overview("Human") + guides(size = "none") | mouse) +
  plot_layout(widths = c(62, 48), heights = unit(row_span("Human") * ROW, "mm")) + FLUSH

# c: age t of human protein (Koopmans 2026 DPFC L1-3) against mouse mRNA-adjusted protein (Takasugi 2024) on
# one-to-one orthologs, points coloured by their local density; the dashed line is the least-squares fit
genes <- read_panel("fig1c")
rho <- read_panel("fig1d")[
  human_unit == "Koopmans 2026 DPFC L1-3" & mouse_cohort == "Mouse cerebrum (Takasugi 2024)" &
    mouse_layer == "mRNA-adjusted protein"
]
XLIM <- c(-5.5, 5.5)
YLIM <- c(-10, 10)
note <- sprintf(
  "ρ = %.2f (%.2f to %.2f)\n%s\nn = %s genes", rho$rho, rho$ci_low, rho$ci_high, fmt_p(rho$p, "Permutation P"),
  format(nrow(genes), big.mark = ",")
)
# density of each point among its neighbours (two-dimensional kernel estimate on a 200 x 200 grid), scaled to 1; points
# are drawn in order of density so that the dense core sits on top
point_density <- function(x, y) {
  k <- MASS::kde2d(x, y, n = 200)
  d <- k$z[cbind(findInterval(x, k$x, all.inside = TRUE), findInterval(y, k$y, all.inside = TRUE))]
  d / max(d)
}
genes[, density := point_density(human_t, mouse_t)]
setorder(genes, density)
scatter <- ggplot(genes, aes(human_t, mouse_t)) +
  ZERO_X + ZERO_Y +
  geom_point(aes(colour = density), size = 0.8, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = INK, linetype = "22", linewidth = 0.8, xseq = XLIM) +
  annotate("text", x = XLIM[1], y = YLIM[2], label = note, hjust = 0, vjust = 1, size = PT, lineheight = 1) +
  scale_colour_gradientn(colours = c(FOCAL_LIGHT, FOCAL, FOCAL_DARK), guide = "none") +
  scale_x_continuous(labels = MINUS_LABELS) +
  scale_y_continuous(labels = MINUS_LABELS) +
  coord_cartesian(xlim = XLIM, ylim = YLIM) +
  labs(x = "Human protein age t\n(Koopmans 2026 DPFC L1-3)", y = "Mouse mRNA-adjusted protein\nage t (Takasugi 2024)") +
  theme_fig()

# d: Spearman rho of every human protein unit (columns) with the mouse age effects of three layers in two cohorts
# (rows); large dots, age permutation P < 0.05; units ordered by rho of the adjusted layer of the cerebral cohort
COHORTS <- c("Mouse cerebrum (Takasugi 2024)" = "Takasugi", "Mouse striatum (Langfelder 2016)" = "Langfelder")
d <- read_panel("fig1d")
d[, mouse_layer := factor(mouse_layer, c("Protein", "mRNA-adjusted protein", "mRNA"),
                          c("Protein", "mRNA-adjusted\nprotein", "mRNA"))]
d[, cohort := factor(COHORTS[mouse_cohort], rev(COHORTS))]
ranked <- d[mouse_layer == "mRNA-adjusted\nprotein" & cohort == "Takasugi"][order(-rho), human_unit]
d[, human_unit := factor(human_unit, ranked)]
heat <- ggplot(d, aes(human_unit, cohort)) +
  geom_point(aes(fill = rho, size = p < 0.05), shape = 21, colour = MATRIX_EDGE, stroke = 0.5) +
  facet_grid(mouse_layer ~ ., switch = "y") +
  scale_fill_div(0.3, "Spearman ρ") +
  scale_size_manual(values = c("FALSE" = 2.2, "TRUE" = 5.4), labels = c("FALSE" = "P ≥ 0.05", "TRUE" = "P < 0.05"),
                    name = NULL,
                    guide = guide_legend(order = 2, override.aes = list(fill = "#BDBDBD", colour = MATRIX_EDGE))) +
  scale_y_discrete(position = "right") +
  labs(x = "Human protein unit", y = "Mouse cohort") +
  theme_matrix() +
  BOXED +
  theme(
    axis.title = element_text(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 2, unit = "mm")),
    axis.text.y.right = element_text(margin = margin(l = 2, unit = "mm")), strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, lineheight = 0.95, margin = margin(r = 2, unit = "mm")),
    legend.position = "bottom"
  )

# b: discovery and validation tiers, mouse to human above and within human below
tier_boxes <- data.table(
  id = 1:4, x0 = c(14, 14, 14, 14), x1 = c(98, 98, 98, 98), y0 = c(76, 46, 26, 1), y1 = c(99, 70, 44, 19),
  fill = c(FOCAL_LIGHT, CONTRAST_LIGHT, FOCAL_LIGHT, CONTRAST_LIGHT),
  head = c("Discovery", "Validation", "Discovery", "Validation"),
  body = c("mouse units\n(protein, mRNA, neurons)",
           "human protein and mRNA\n(striatum, cortex, hippocampus,\ncerebellum, neurons)",
           "BrainGVEX", "other human cohorts")
)
tiers <- ggplot() +
  boxes(tier_boxes) +
  geom_segment(data = data.table(x = 56, y = c(75, 25), xend = 56, yend = c(71, 20.5)),
               aes(x, y, xend = xend, yend = yend), colour = EDGE, linewidth = 0.8,
               arrow = arrow(length = unit(2.2, "mm"), type = "closed")) +
  annotate("text", x = 5, y = c(72, 22), label = c("Mouse to human", "Within human"), size = PT, angle = 90) +
  frame()

fig <- panel_grid(list(a = a, b = tiers, c = scatter, d = heat), design = "AAC\nBDD", widths = c(56, 54, 73),
                  heights = c(78, 92))
save_fig(fig, "fig1", 340)
