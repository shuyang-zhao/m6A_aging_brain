# Additional file 1 Fig. S2: the m6A maps, the sensitivity of the m6A signature and the mouse load to the map, and the
# quality of the four mouse maps
# reads source_data/sfig2a-f; writes output/figures/sfig2.png and .pdf
source("figures/style.R")

ROWS <- c("mouse_protein|protein_only" = "Protein", "mouse_protein|mrna_adjusted" = "Protein given mRNA",
          "mouse_mrna|mrna_only" = "mRNA", "mouse_neurons|mrna_only" = "Neuronal mRNA",
          "human_protein|protein_only" = "Protein", "human_striatum|mrna_only" = "Striatum",
          "human_cortex|mrna_only" = "Cortex", "human_hippocampus|mrna_only" = "Hippocampus",
          "human_cerebellum|mrna_only" = "Cerebellum", "human_neurons|mrna_only" = "Neurons")
DEF <- c("brain peaks", "cell line", "union", "hippocampus cyto", "hippocampus syn", "cortex scDART", "neuron DART")
DEF_COL <- setNames(c(FOCAL, CONTRAST, GREY, FOCAL, RED, BLUE, CONTRAST), DEF)
rows_for <- function(d, species) {
  d <- d[grepl(paste0("^", species), context)]
  d[, `:=`(row = factor(ROWS[paste(context, layer, sep = "|")], unique(ROWS[grepl(species, names(ROWS))])),
           def = factor(definition, DEF))]
  droplevels(d[!is.na(row)])
}

# a: genes with an m6A site in the human brain map, the cell-line maps or both (Venn diagram)
v <- read_panel("sfig2a")
n <- setNames(v$genes, v$region)
circles <- data.table(x0 = c(-0.55, 0.55), y0 = 0, r = 1, fill = c(FOCAL_LIGHT, CONTRAST_LIGHT),
                      edge = c(FOCAL, CONTRAST))
a <- ggplot() +
  ggforce::geom_circle(data = circles, aes(x0 = x0, y0 = y0, r = r, fill = I(fill), colour = I(edge)), alpha = 0.7,
                       linewidth = 0.8) +
  annotate("text", x = c(-1, 0, 1), y = 0, size = PT, colour = INK,
           label = format(n[c("brain peaks only", "both", "cell line only")], big.mark = ",")) +
  annotate("text", x = c(-0.55, 0.55), y = 1.22, size = PT, colour = INK,
           label = c("Brain map", "Cell-line maps")) +
  coord_equal(xlim = c(-1.65, 1.65), ylim = c(-1.1, 1.4)) +
  labs(x = "Human genes with m6A sites") +
  theme_void(base_size = TEXT, base_family = FONT) +
  theme(axis.title.x = element_text(), plot.margin = margin(4.4, 5, 3, 2, "mm"))

# b: genes with an m6A site in each combination of the four mouse maps (UpSet; the eight largest combinations)
MOUSE <- c("hippocampus cyto", "hippocampus syn", "cortex scDART", "neuron DART")
u <- read_panel("sfig2b")[order(-genes)][1:8]
u[, combo := factor(seq_len(.N))]
b <- upset(
  u, setNames(nm = MOUSE),
  geom_col(width = 0.65, fill = FOCAL, colour = EDGE, linewidth = 0.5),
  geom_text(aes(label = format(genes, big.mark = ",")), vjust = -0.4, size = PT, colour = INK),
  labs(y = "Mouse genes with m6A sites"),
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)), labels = scales::label_comma()),
  dot_colour = \(d) DEF_COL[as.character(d$set)], dot_size = 3.2, heights = c(1.4, 1)
)

# c: human m6A signature by map (one dot per map; filled, q < 0.05); the human load by map is Figure 3d
d <- rows_for(read_panel("sfig2c"), "human")
cc <- ggplot(d, aes(z, row)) +
  ZERO_X +
  geom_point(data = d[q < 0.05], aes(fill = def), colour = EDGE, shape = 21, size = 4.4, stroke = 0.5,
             show.legend = TRUE) +
  geom_point(data = d[q >= 0.05], colour = GREY, fill = "white", shape = 21, size = 4.4, stroke = 0.5) +
  scale_fill_manual(values = DEF_COL, limits = levels(d$def), name = NULL) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(labels = MINUS_LABELS) +
  labs(x = "Human m6A signature (cameraPR z)", y = NULL) +
  theme_rows() +
  theme(legend.position = "bottom")

# d: mouse m6A load by map (estimate, 95% CI)
d <- rows_for(read_panel("sfig2d"), "mouse")
dodge <- position_dodge(0.7)
dd <- ggplot(d, aes(effect_pct, row, fill = def)) +
  ZERO_X +
  geom_linerange(aes(xmin = ci_low_pct, xmax = ci_high_pct), colour = CI, linewidth = 0.6, position = dodge) +
  geom_point(colour = EDGE, shape = 21, size = 3.4, stroke = 0.5, position = dodge) +
  scale_fill_manual(values = DEF_COL, name = NULL) +
  scale_y_discrete(limits = rev, labels = scales::label_wrap(10)) +
  scale_x_continuous(labels = MINUS_LABELS) +
  labs(x = "Age effect of the mouse m6A load (%)", y = NULL) +
  theme_rows() +
  theme(legend.position = "bottom", legend.key.size = unit(4.4, "mm")) +
  guides(fill = guide_legend(nrow = 3))

# e: mouse m6A signature by map (cameraPR z; bold, q < 0.05)
d <- rows_for(read_panel("sfig2e"), "mouse")
e <- ggplot(d, aes(row, def, fill = z)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = cell_number(z), fontface = fifelse(q < 0.05, "bold", "plain")), size = PT, colour = INK) +
  scale_fill_div(4, "Mouse m6A signature (cameraPR z)") +
  scale_x_discrete(position = "top", labels = scales::label_wrap(10), expand = c(0, 0)) +
  scale_y_discrete(limits = rev, expand = c(0, 0)) +
  theme_heat() +
  theme(legend.position = "bottom")

# f: the four mouse maps: share of sites in a DRACH motif, share placed on a canonical transcript, sites per marked gene
# and share at a GLORI position
q <- melt(read_panel("sfig2f"), id.vars = "map", measure.vars = c("drach", "placed", "sites_per_gene", "at_glori"))
q[, `:=`(variable = factor(variable, c("drach", "placed", "sites_per_gene", "at_glori"),
                           c("In a DRACH motif", "On a canonical transcript", "Sites per marked gene",
                             "At a GLORI position")),
         map = factor(map, rev(c("hippocampus cyto", "hippocampus syn", "cortex scDART", "neuron DART"))))]
f <- ggplot(q, aes(value, map, fill = map)) +
  geom_col(width = 0.7, colour = EDGE, linewidth = 0.5) +
  geom_text(aes(label = fifelse(variable == "Sites per marked gene", sprintf("%.1f", value),
                                sprintf("%.0f%%", 100 * value))), hjust = -0.15, size = PT, colour = INK) +
  facet_grid(. ~ variable, scales = "free_x") +
  scale_fill_manual(values = DEF_COL, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35)), breaks = NULL) +
  labs(x = NULL, y = NULL) +
  theme_rows() +
  BOXED +
  theme(panel.spacing = unit(4, "mm"))

fig <- panel_grid(list(a = a, b = b, c = cc, d = dd, e = e, f = f), design = "ABB\nCDE\nFFF", widths = c(58, 50, 75),
                  heights = c(56, 78, 34))
save_fig(fig, "sfig2", 336)
