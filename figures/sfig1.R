# Additional file 1 Fig. S1: resource and concordance, protein ageing beyond mRNA shared between species
# reads source_data/sfig1a and sfig1b; writes output/figures/sfig1.png and .pdf
source("figures/style.R")

# a: age t of human protein (Koopmans 2026 DPFC L1-3) against mouse mRNA-adjusted protein (Takasugi 2024) on
# one-to-one orthologs, points coloured by their local density; the dashed line is the least-squares fit
genes <- read_panel("sfig1a")
rho <- read_panel("sfig1b")[human_unit == "Koopmans 2026 DPFC L1-3" & mouse_layer == "mRNA-adjusted protein"]
XLIM <- c(-5.5, 5.5)
YLIM <- c(-10, 10)
note <- sprintf("ρ = %.2f (%.2f to %.2f)\n%s\nn = %s genes", rho$rho, rho$ci_low, rho$ci_high,
                fmt_p(rho$p, "Permutation P"), format(nrow(genes), big.mark = ","))
# density of each point among its neighbours (kernel estimate on a 200 x 200 grid), scaled to 1; the dense core is
# drawn last
point_density <- function(x, y) {
  k <- MASS::kde2d(x, y, n = 200)
  d <- k$z[cbind(findInterval(x, k$x, all.inside = TRUE), findInterval(y, k$y, all.inside = TRUE))]
  d / max(d)
}
genes[, density := point_density(human_t, mouse_t)]
setorder(genes, density)
a <- ggplot(genes, aes(human_t, mouse_t)) +
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

# b: Spearman rho of every human protein unit (columns) with the mouse age effects of three layers (rows) in the mouse
# cohort; large dots, age permutation P < 0.05; units ordered by rho of the mRNA-adjusted layer
d <- read_panel("sfig1b")
d[, mouse_layer := factor(mouse_layer, c("Protein", "mRNA-adjusted protein", "mRNA"),
                          c("Protein", "mRNA-adjusted\nprotein", "mRNA"))]
d[, human_unit := factor(human_unit, d[mouse_layer == "mRNA-adjusted\nprotein"][order(-rho), human_unit])]
b <- ggplot(d, aes(human_unit, mouse_layer)) +
  geom_point(aes(fill = rho, size = p < 0.05), shape = 21, colour = MATRIX_EDGE, stroke = 0.5) +
  scale_fill_div(0.3, "Spearman ρ") +
  scale_size_manual(values = c("FALSE" = 2.2, "TRUE" = 5.4), labels = c("FALSE" = "P ≥ 0.05", "TRUE" = "P < 0.05"),
                    name = NULL,
                    guide = guide_legend(order = 2, override.aes = list(fill = "#BDBDBD", colour = MATRIX_EDGE))) +
  scale_y_discrete(limits = rev) +
  labs(x = "Human protein unit", y = "Mouse, Takasugi 2024") +
  theme_matrix() +
  BOXED +
  theme(
    axis.title = element_text(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 2, unit = "mm")),
    axis.text.y = element_text(lineheight = 0.95, margin = margin(r = 2, unit = "mm")),
    legend.position = "bottom"
  )

fig <- panel_grid(list(a = a, b = b), design = "AB", widths = c(74, 109), heights = 78)
save_fig(fig, "sfig1", 156)
