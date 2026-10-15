# Additional file 1 Fig. S1: the other modifications (Ψ, m5C, m7G) beside m6A
# reads source_data/sfig1a and sfig1b; writes output/figures/sfig1.png and .pdf
source("figures/style.R")

MOD <- c(m6A = FOCAL, "Ψ" = RED, m5C = BLUE, m7G = CONTRAST)
ROWS <- c("Mouse protein", "Mouse protein given mRNA", "Mouse mRNA", "Mouse neuronal mRNA", "Human protein",
          "Human striatum mRNA", "Human cortex mRNA", "Human hippocampus mRNA", "Human cerebellum mRNA",
          "Human neuronal mRNA")

# a: signature of the marked genes against unmarked genes (cameraPR z, combined over units; bold, q < 0.05)
d <- read_panel("sfig1a")[, `:=`(row = factor(row, rev(ROWS)), mod = factor(modification, names(MOD)))]
a <- ggplot(d, aes(mod, row, fill = z)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = cell_number(z), fontface = fifelse(q < 0.05, "bold", "plain")), size = PT, colour = INK) +
  scale_fill_div(4, "Signature of marked genes (cameraPR z)") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_heat() +
  theme(legend.position = "bottom")

# b: age effect of the modification load per context (z of the random-effects coefficient over cohorts; filled,
# q < 0.05), one row per mark
d <- read_panel("sfig1b")[, mod := factor(modification, names(MOD))]
b <- ggplot(d, aes(z, mod)) +
  annotate("rect", xmin = -1.96, xmax = 1.96, ymin = -Inf, ymax = Inf, fill = "#F3F3F3") +
  ZERO_X +
  geom_point(aes(colour = I(fifelse(q < 0.05, EDGE, GREY)),
                 fill = I(fifelse(q < 0.05, MOD[as.character(mod)], "white"))),
             shape = 21, size = 4, stroke = 0.5, position = position_jitter(height = 0.18, width = 0, seed = 1)) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(labels = MINUS_LABELS) +
  labs(x = "Modification load z across contexts (grey band, |z| < 1.96)", y = NULL) +
  theme_rows()

fig <- panel_grid(list(a = a, b = b), design = "AB", widths = c(88, 95))
save_fig(fig, "sfig1", 176)
