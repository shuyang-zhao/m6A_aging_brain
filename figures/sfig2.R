# Additional file 1 Fig. S2: the m6A maps, the mouse load and signature by map, and the quality of the maps
# reads source_data/sfig2a-d; writes output/figures/sfig2.png and .pdf
source("figures/style.R")

LAYER_ROWS <- c("mouse_protein|protein_only" = "Protein", "mouse_protein|mrna_adjusted" = "mRNA-adjusted protein",
                "mouse_mrna|mrna_only" = "mRNA", "mouse_neurons|mrna_only" = "Neuronal mRNA")
MAPS <- c("hippocampus cyto", "hippocampus syn", "union")
MAP_LABEL <- c("Hippocampus cyto", "Hippocampus syn", "Union")
MAP_FILL <- c("hippocampus cyto" = tint(MOUSE, 0.62), "hippocampus syn" = tint(MOUSE, 0.32), union = MOUSE,
              "human brain peaks" = HUMAN)
by_map <- function(d) {
  d[, `:=`(row = factor(LAYER_ROWS[paste(context, layer, sep = "|")], rev(LAYER_ROWS)), map = factor(definition, MAPS))]
}

# a: mouse genes with an m6A site in the cytoplasmic map, the synaptosomal map or both; circle areas proportional to
# the gene counts of the maps, the overlap to the genes in both
v <- read_panel("sfig2a")
n <- c(cyto = v[`hippocampus cyto` & !`hippocampus syn`, genes], both = v[`hippocampus cyto` & `hippocampus syn`, genes],
       syn = v[!`hippocampus cyto` & `hippocampus syn`, genes])
r <- 20 * sqrt(c(n[["cyto"]] + n[["both"]], n[["syn"]] + n[["both"]]) / (n[["cyto"]] + n[["both"]]))
lens <- function(g) {
  r[1]^2 * acos((g^2 + r[1]^2 - r[2]^2) / (2 * g * r[1])) + r[2]^2 * acos((g^2 + r[2]^2 - r[1]^2) / (2 * g * r[2])) -
    0.5 * sqrt((-g + r[1] + r[2]) * (g + r[1] - r[2]) * (g - r[1] + r[2]) * (g + r[1] + r[2]))
}
gap <- uniroot(\(g) lens(g) - pi * r[1]^2 * n[["both"]] / (n[["cyto"]] + n[["both"]]), c(abs(diff(r)) + 1e-6, sum(r) - 1e-6))$root
circles <- data.table(map = MAPS[1:2], x0 = c(-gap / 2, gap / 2), r = r)
set_total <- c(n[["cyto"]] + n[["both"]], n[["syn"]] + n[["both"]], sum(n))
a <- ggplot() +
  ggforce::geom_circle(data = circles, aes(x0 = x0, y0 = 0, r = r, fill = map), colour = EDGE, linewidth = lw(0.5),
                       alpha = 0.7) +
  note(0, 0, thousands(n[["both"]])) +
  note(-gap / 2 - r[1] - 1, 0, thousands(n[["cyto"]]), hjust = 1) +
  note(gap / 2 + r[2] + 1, 0, thousands(n[["syn"]]), hjust = 0) +
  key_rows(paste0(MAP_LABEL, ", ", thousands(set_total), " genes"), MAP_FILL[MAPS], x = -r[1] - 4, y = -r[1] - 6,
           step = 4.5, width = 3, square = TRUE) +
  scale_fill_manual(values = MAP_FILL, guide = "none") +
  coord_fixed(xlim = c(-gap / 2 - r[1] - 11, gap / 2 + r[2] + 11), ylim = c(-r[1] - 17, r[1] + 3), expand = FALSE) +
  theme_void_pub()

# b: age effect of the mouse m6A load by map (95% CI as the thin line)
b <- by_map(read_panel("sfig2b"))
dodge <- position_dodge(0.8, reverse = TRUE)
b <- ggplot(b, aes(effect_pct, row, fill = map)) +
  geom_vline(xintercept = 0, linewidth = lw(0.5)) +
  geom_col(width = 0.72, position = dodge, colour = EDGE, linewidth = lw(0.4)) +
  geom_linerange(aes(xmin = ci_low_pct, xmax = ci_high_pct), orientation = "y", position = dodge, colour = EDGE,
                 linewidth = lw(0.5)) +
  key_rows(MAP_LABEL, MAP_FILL[MAPS], x = 0.25, y = 4.3, step = 0.28, width = 0.18, square = TRUE) +
  scale_fill_manual(values = MAP_FILL, guide = "none") +
  scale_x_continuous(labels = MINUS_LABELS, limits = c(-2.5, 1.7), breaks = -2:1) +
  scale_y_discrete(labels = scales::label_wrap(14)) +
  labs(x = "Age effect of the mouse m6A load (% per SD)", y = NULL) +
  theme_rows()

# c: mouse m6A signature of the marked genes by map (cameraPR z; bold, q < 0.05)
c_ <- by_map(read_panel("sfig2c"))
c_ <- ggplot(c_, aes(map, row, fill = z)) +
  geom_tile(colour = "white", linewidth = lw(0.6)) +
  geom_text(aes(label = number(z, 1), fontface = fifelse(q < 0.05, "bold", "plain")), size = pt_size(6)) +
  scale_fill_signed(6, "Mouse m6A signature (cameraPR z)", light = TRUE) +
  scale_x_discrete(position = "top", labels = scales::label_wrap(10)(MAP_LABEL), expand = c(0, 0)) +
  scale_y_discrete(labels = scales::label_wrap(14), expand = c(0, 0)) +
  theme_tile() +
  theme(legend.position = "bottom", legend.margin = margin(t = 1, unit = "mm"))

# d: the maps by sites (peaks), share in a DRACH motif, share on a canonical transcript and sites per gene
q <- read_panel("sfig2d")
q[, `:=`(map = factor(map, rev(names(MAP_FILL)[c(1, 2, 4)])), sites = as.numeric(sites))]
q <- melt(q, id.vars = "map", measure.vars = c("sites", "drach", "placed", "sites_per_gene"))
q[, variable := factor(variable, c("sites", "drach", "placed", "sites_per_gene"),
                       c("Sites (peaks)", "In a DRACH motif", "On a canonical transcript", "Sites per gene"))]
q[, label := fcase(variable == "Sites (peaks)", thousands(value), variable == "Sites per gene", sprintf("%.1f", value),
                   default = sprintf("%.1f%%", 100 * value))]
dd <- ggplot(q, aes(value, map, fill = map)) +
  geom_col(width = 0.7, colour = EDGE, linewidth = lw(0.4)) +
  geom_text(aes(label = label), hjust = -0.15, size = pt_size(6)) +
  facet_wrap(~variable, nrow = 2, scales = "free_x") +
  scale_fill_manual(values = MAP_FILL, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.4)), breaks = NULL) +
  scale_y_discrete(labels = c("hippocampus cyto" = "Mouse cyto", "hippocampus syn" = "Mouse syn",
                              "human brain peaks" = "Human brain\npeaks")) +
  labs(x = NULL, y = NULL) +
  theme_box() +
  theme(axis.ticks.x = element_blank())

fig <- fig_compose(a, b, c_, dd, design = "AAABBBBBB\nCCCCDDDDD", heights = c(56, 60))
save_fig(fig, "sfig2", 120)
