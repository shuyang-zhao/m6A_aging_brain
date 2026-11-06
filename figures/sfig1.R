# Additional file 1 Fig. S1: protein aging beyond mRNA is shared between mouse and human
# reads source_data/sfig1a and sfig1b; writes output/figures/sfig1.png and .pdf
source("figures/style.R")

# a: age t of one-to-one orthologs, human protein (Koopmans 2026 DPFC L1-3) against mouse mRNA-adjusted protein
# (Takasugi 2024); points coloured by local density, dashed least-squares line
genes <- by_density(read_panel("sfig1a"), "human_t", "mouse_t")
rho <- read_panel("sfig1b")[human_unit == "Koopmans 2026 DPFC L1-3" & mouse_layer == "mRNA-adjusted protein"]
XLIM <- c(-5.5, 5.5)
YLIM <- c(-10, 10)
a <- ggplot(genes, aes(human_t, mouse_t)) +
  ref_line(xintercept = 0) +
  ref_line(yintercept = 0) +
  rasterise_pts(geom_point(aes(colour = .density), size = 0.5, stroke = 0)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linetype = "22", linewidth = lw(0.6),
              fullrange = TRUE) +
  note_corner(lab_stats(rho$rho, rho$p, nrow(genes), stat = "italic(rho)"), parse = TRUE) +
  note_corner(sprintf("95%% CI of ρ %s to %s", number(rho$ci_low), number(rho$ci_high)), vjust = 2.9) +
  scale_density(MOUSE) +
  scale_x_continuous(labels = MINUS_LABELS) +
  scale_y_continuous(labels = MINUS_LABELS) +
  coord_cartesian(xlim = XLIM, ylim = YLIM) +
  labs(x = "Human protein age t\n(Koopmans 2026 DPFC L1-3)", y = "Mouse mRNA-adjusted protein age t\n(Takasugi 2024)") +
  theme_pub()

# b: Spearman rho of every human protein unit with the age effects of the mouse cerebrum on three layers; dot size,
# permutation P; units ordered by rho of the mRNA-adjusted layer
d <- read_panel("sfig1b")
d[, layer := factor(mouse_layer, c("Protein", "mRNA-adjusted protein", "mRNA"),
                    c("Protein", "mRNA-adjusted\nprotein", "mRNA"))]
d[, unit := factor(human_unit, rev(d[mouse_layer == "mRNA-adjusted protein"][order(-rho), human_unit]))]
b <- ggplot(d, aes(layer, unit, fill = rho, size = -log10(p))) +
  geom_point(shape = 21, colour = EDGE, stroke = st(0.4)) +
  scale_fill_signed(0.3, "Spearman ρ") +
  scale_size_area(max_size = 3.2, name = expression(Permutation~italic(P)), breaks = -log10(c(0.5, 0.05, 0.005)),
                  labels = c("0.5", "0.05", "0.005"),
                  guide = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1)) +
  scale_x_discrete(position = "top") +
  labs(x = "Mouse cerebrum (Takasugi 2024)", y = "Human protein unit") +
  theme_tile() +
  theme(axis.title = element_text(), axis.title.y = element_text(margin = margin(r = 2, unit = "mm")),
        panel.grid.major = element_line(linewidth = lw(0.25), colour = LIGHT), legend.position = "bottom",
        legend.box = "vertical", legend.key.spacing.x = unit(1.5, "mm"))

fig <- fig_compose(a, b, guide_area(), design = "AB\nCB", widths = c(98, 85), heights = c(84, 26), guides = "collect")
save_fig(fig, "sfig1", 112)
