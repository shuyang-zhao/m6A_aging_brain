# Additional file 1 Fig. S2: the m6A maps, the mouse load and signature by map, and the quality of the maps
# reads source_data/sfig2a-d; writes output/figures/sfig2.png and .pdf
source("figures/style.R")

ROWS <- c("mouse_protein|protein_only" = "Protein", "mouse_protein|mrna_adjusted" = "mRNA-adjusted protein",
          "mouse_mrna|mrna_only" = "mRNA", "mouse_neurons|mrna_only" = "Neuronal mRNA")
MAPS <- c("hippocampus cyto", "hippocampus syn", "union")
MAP_LABEL <- c("Hippocampus cyto", "Hippocampus syn", "Union")
MAP_COL <- c("hippocampus cyto" = FOCAL, "hippocampus syn" = BLUE, union = GREY, "human brain peaks" = RED)
by_map <- function(d) {
  d[, `:=`(row = factor(ROWS[paste(context, layer, sep = "|")], rev(ROWS)), map = factor(definition, MAPS))]
}

# a: mouse genes with an m6A site in the cytoplasmic map, the synaptic map or both (Venn diagram)
v <- read_panel("sfig2a")
n <- c(only_cyto = v[`hippocampus cyto` & !`hippocampus syn`, genes], both = v[`hippocampus cyto` & `hippocampus syn`, genes],
       only_syn = v[!`hippocampus cyto` & `hippocampus syn`, genes])
circles <- data.table(x0 = c(-0.55, 0.55), y0 = 0, r = 1, fill = c(FOCAL_LIGHT, BLUE_LIGHT), edge = c(FOCAL, BLUE))
a <- ggplot() +
  ggforce::geom_circle(data = circles, aes(x0 = x0, y0 = y0, r = r, fill = I(fill), colour = I(edge)), alpha = 0.7,
                       linewidth = 0.8) +
  annotate("text", x = c(-1, 0, 1), y = 0, size = PT, colour = INK, label = format(n, big.mark = ",")) +
  annotate("text", x = c(-0.9, 0.9), y = 1.2, size = PT, colour = INK, label = c("Hippocampus\ncyto", "Hippocampus\nsyn"),
           lineheight = 0.95) +
  coord_equal(xlim = c(-1.65, 1.65), ylim = c(-1.1, 1.7)) +
  labs(x = "Mouse genes with m6A sites") +
  theme_void(base_size = TEXT, base_family = FONT) +
  theme(axis.title.x = element_text(), plot.margin = margin(4.4, 5, 3, 2, "mm"))

# b: age effect of the mouse m6A load by map (estimate, 95% CI)
d <- by_map(read_panel("sfig2b"))[, map := factor(map, rev(MAPS))]
dodge <- position_dodge(0.7)
b <- ggplot(d, aes(effect_pct, row, fill = map)) +
  ZERO_X +
  geom_linerange(aes(xmin = ci_low_pct, xmax = ci_high_pct), colour = CI, linewidth = 0.6, position = dodge) +
  geom_point(colour = EDGE, shape = 21, size = 3.4, stroke = 0.5, position = dodge) +
  scale_fill_manual(values = MAP_COL, name = NULL, breaks = MAPS, labels = MAP_LABEL) +
  scale_x_continuous(labels = MINUS_LABELS) +
  scale_y_discrete(labels = scales::label_wrap(14)) +
  labs(x = "Age effect of the mouse m6A load (% per SD)", y = NULL) +
  theme_rows() +
  theme(legend.position = "bottom")

# c: mouse m6A signature of the marked genes by map (cameraPR z; bold, q < 0.05)
d <- by_map(read_panel("sfig2c"))
cc <- ggplot(d, aes(map, row, fill = z)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = cell_number(z), fontface = fifelse(q < 0.05, "bold", "plain")), size = PT, colour = INK) +
  scale_fill_div(5, "Mouse m6A signature (cameraPR z)") +
  scale_x_discrete(position = "top", labels = scales::label_wrap(10)(MAP_LABEL), expand = c(0, 0)) +
  scale_y_discrete(labels = scales::label_wrap(14), expand = c(0, 0)) +
  theme_heat() +
  theme(legend.position = "bottom")

# d: the three maps: sites (peaks), share in a DRACH motif, share placed on a canonical transcript, sites per gene
q <- read_panel("sfig2d")
q[, `:=`(map = factor(map, rev(names(MAP_COL)[c(1, 2, 4)])), sites = as.numeric(sites))]
q <- melt(q, id.vars = "map", measure.vars = c("sites", "drach", "placed", "sites_per_gene"))
q[, variable := factor(variable, c("sites", "drach", "placed", "sites_per_gene"),
                       c("Sites (peaks)", "In a DRACH motif", "On a canonical transcript", "Sites per gene"))]
q[, label := fcase(variable == "Sites (peaks)", scales::comma(value, accuracy = 1),
                   variable == "Sites per gene", sprintf("%.1f", value), default = sprintf("%.1f%%", 100 * value))]
dd <- ggplot(q, aes(value, map, fill = map)) +
  geom_col(width = 0.7, colour = EDGE, linewidth = 0.5) +
  geom_text(aes(label = label), hjust = -0.15, size = PT, colour = INK) +
  facet_grid(. ~ variable, scales = "free_x") +
  scale_fill_manual(values = MAP_COL, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.5)), breaks = NULL) +
  scale_y_discrete(labels = c("hippocampus cyto" = "Mouse hippocampus cyto", "hippocampus syn" = "Mouse hippocampus syn",
                              "human brain peaks" = "Human brain peaks")) +
  labs(x = NULL, y = NULL) +
  theme_rows() +
  BOXED +
  theme(panel.spacing = unit(4, "mm"))

fig <- panel_grid(list(a = a, b = b, c = cc, d = dd), design = "ABC\nDDD", widths = c(46, 70, 67), heights = c(60, 30))
save_fig(fig, "sfig2", 180)
