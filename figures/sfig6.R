# Additional file 1 Fig. S6: m6A against the other marks (Ψ, m5C, m7G)
# reads source_data/sfig6a and sfig6b; writes output/figures/sfig6.png and .pdf
source("figures/style.R")

MARKS <- c("m6A", "Ψ", "m5C", "m7G")
ROWS <- c("Mouse protein", "Mouse mRNA-adjusted protein", "Mouse mRNA", "Mouse neuronal mRNA", "Human protein",
          "Human striatum mRNA", "Human cortex mRNA", "Human hippocampus mRNA",
          "Human cerebellum mRNA", "Human neuronal mRNA")

# a: signature of the marked genes against unmarked genes (cameraPR z, combined over units; bold, q < 0.05)
d <- read_panel("sfig6a")[, `:=`(row = factor(row, rev(ROWS)), mark = factor(modification, MARKS))]
a <- ggplot(d, aes(mark, row, fill = z)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = cell_number(z), fontface = fifelse(q < 0.05, "bold", "plain")), size = PT, colour = INK) +
  scale_fill_div(4, "Signature of marked genes (cameraPR z)") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_heat() +
  theme(legend.position = "bottom")

# b: age effect of the modification load per context (z of the random-effects coefficient over cohorts), one row per
# mark; filled in the species colour where q < 0.05
GROUP <- c("Mouse, q < 0.05", "Human, q < 0.05", "q ≥ 0.05")
d <- read_panel("sfig6b")[, mark := factor(modification, MARKS)]
d[, group := factor(fifelse(q >= 0.05, GROUP[3], fifelse(grepl("^Mouse", row), GROUP[1], GROUP[2])), GROUP)]
b <- ggplot(d, aes(z, mark)) +
  annotate("rect", xmin = -1.96, xmax = 1.96, ymin = -Inf, ymax = Inf, fill = "#F3F3F3") +
  ZERO_X +
  geom_point(aes(fill = group, colour = group), shape = 21, size = 4, stroke = 0.5,
             position = position_jitter(height = 0.18, width = 0, seed = 1)) +
  scale_fill_manual(values = setNames(c(FOCAL, RED, "white"), GROUP), name = NULL) +
  scale_colour_manual(values = setNames(c(EDGE, EDGE, GREY), GROUP), name = NULL) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(labels = MINUS_LABELS) +
  labs(x = "Modification load z across contexts (grey band, |z| < 1.96)", y = NULL) +
  theme_rows() +
  theme(legend.position = "bottom")

fig <- panel_grid(list(a = a, b = b), design = "AB", widths = c(95, 88), heights = 80)
save_fig(fig, "sfig6", 160)
