# Figure 5: an FMRP-YTHDF reading axis tracks m6A-marked protein output in mouse
# reads source_data/fig5b, fig5c, fig5d and fig5e; writes output/figures/fig5.png and .pdf
source("figures/style.R")

ARMS <- c("Output", "Translation", "Stability")

# a: the coupling model, a schematic with simulated individuals: the slope of the score of m6A-linked transcripts on the
# level of a regulator (b_s) against the slope for bound transcripts without m6A (b_u); the three arms and the controls
set.seed(1)
sim <- data.table(level = rnorm(60))
sim <- rbind(sim[, .(level, score = 0.6 * level + rnorm(60, sd = 0.45), set = "m6A-linked transcripts")],
             sim[, .(level, score = 0.05 * level + rnorm(60, sd = 0.45), set = "Bound without m6A")])
SETS <- c("m6A-linked transcripts" = FOCAL, "Bound without m6A" = GREY)
model <- ggplot(sim, aes(level, score, colour = set)) +
  geom_point(aes(fill = set), colour = EDGE, shape = 21, size = 1.8, stroke = 0.5, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 1.2) +
  scale_colour_manual(values = SETS, name = NULL) +
  scale_fill_manual(values = SETS, name = NULL) +
  labs(x = "Regulator level in each individual", y = "Transcript score") +
  theme_fig() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "inside",
        legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1))
arm_boxes <- data.table(
  id = 1:4, x0 = 4, x1 = 96, y0 = c(75, 52, 29, 3), y1 = c(97, 74, 51, 25),
  fill = c(LAYER_LIGHT[c("Protein", "Protein given mRNA", "Bulk RNA")], NS),
  head = c("Output", "Translation", "Stability", "Control RBPs"),
  body = c("protein", "protein given mRNA", "mRNA", "same test; edge above\n95th percentile")
)
arms <- ggplot() + boxes(arm_boxes) + frame()
a <- (model | arms) + plot_layout(widths = c(1.15, 1)) + FLUSH

# b: age effect of every m6A regulator in mouse protein (random effects over cohorts, late versus early age group)
# against its q (Benjamini-Hochberg over all RNA-modification regulators of the context); colour, role of the regulator
v <- read_panel("fig5b")
v[, role := factor(role, c("reader", "writer", "eraser"), c("Reader", "Writer", "Eraser"))]
setorder(v, regulator)
b <- ggplot(v, aes(effect, -log10(q))) +
  geom_hline(yintercept = -log10(0.05), linewidth = 0.4, colour = INK, linetype = "22") +
  ZERO_X +
  geom_point(aes(fill = role, size = regulator == "YTHDF3"), colour = EDGE, shape = 21, stroke = 0.5) +
  ggrepel::geom_text_repel(
    data = v[q < 0.1 | regulator %in% c("YTHDF1", "YTHDF2", "METTL3", "FTO", "ALKBH5")],
    aes(label = regulator, fontface = fifelse(regulator == "YTHDF3", "bold", "plain"),
        point.size = fifelse(regulator == "YTHDF3", 6, 3.8)), size = PT, colour = INK,
    segment.size = 0.4, segment.colour = EDGE, min.segment.length = 0, box.padding = 0.3, nudge_y = 0.3, seed = 1,
    max.overlaps = Inf
  ) +
  scale_fill_manual(values = c(Reader = FOCAL, Writer = BLUE, Eraser = CONTRAST), name = NULL) +
  scale_size_manual(values = c("FALSE" = 3.8, "TRUE" = 6), guide = "none") +
  scale_x_continuous(labels = MINUS_LABELS) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = "Age effect in mouse protein (%, late versus early)", y = expression(-log[10] ~ italic(q))) +
  theme_fig() +
  theme(legend.position = "bottom")

# c: the reader network; regulators on a circle around the three arms; a line is an m6A-specific edge in at least one
# unit (thick, also shared across units, mashr; mustard, positive; teal, negative); node size, number of specific
# edges. The values of all edges are in Additional file 1 Fig. S5b
cells <- read_panel("fig5c")[, .(name = row, arm, value, spec = m6A_specific, shared)]
REGS <- c("YTHDF1", "YTHDF2", "YTHDF3", "YTHDC1", "YTHDC2", "FMR1", "FXR1", "FXR2", "PRRC2A", "HNRNPC", "RBMX",
          "IGF2BP2", "YTHDF module", "FMRP family", "Writer complex", "Erasers")
FAMILY <- c(`YTH readers` = FOCAL, `FMRP family` = RED, `Other readers` = BLUE, Modules = CONTRAST,
            `Writers, erasers` = GREY)
angle <- pi / 2 - 2 * pi * (seq_along(REGS) - 0.5) / length(REGS)
nodes <- data.table(name = REGS, family = factor(rep(names(FAMILY), c(5, 3, 4, 2, 2)), names(FAMILY)),
                    x = cos(angle), y = sin(angle))
hubs <- data.table(arm = ARMS, xend = c(0, -0.4, 0.4), yend = c(0.42, -0.24, -0.24))
links <- merge(cells[spec %in% TRUE, .(name, arm, shared, positive = value > 0)], nodes, by = "name")
links <- merge(links, hubs, by = "arm")
nodes[, degree := links[, .N, by = name][match(nodes$name, name), N]]
nodes[is.na(degree), degree := 0L]
nodes[, `:=`(lx = 1.16 * x, ly = 1.16 * y, hjust = fifelse(x > 0.1, 0, fifelse(x < -0.1, 1, 0.5)))]
cc <- ggplot() +
  geom_segment(data = links, aes(x, y, xend = xend, yend = yend, linewidth = shared,
                                 colour = I(fifelse(positive, CONTRAST, FOCAL_DARK))),
               alpha = 0.85, lineend = "round") +
  geom_point(data = nodes, aes(x, y, size = degree, fill = family), colour = EDGE, shape = 21, stroke = 0.5) +
  geom_text(data = nodes, aes(lx, ly, label = name, hjust = hjust), size = PT, colour = INK) +
  geom_label(data = hubs, aes(xend, yend, label = arm), size = PT, fill = "#F0F0F0", colour = INK, border.colour = NA,
             label.padding = unit(1.4, "mm"), label.r = unit(2, "mm")) +
  scale_linewidth_manual(values = c("FALSE" = 0.9, "TRUE" = 2.6),
                         labels = c("FALSE" = "m6A-specific", "TRUE" = "Shared"), name = NULL) +
  scale_fill_manual(values = FAMILY, guide = "none") +
  scale_size(range = c(3.2, 8.4), guide = "none") +
  coord_equal(xlim = c(-1.75, 1.6), ylim = c(-1.25, 1.25), clip = "off") +
  guides(linewidth = guide_legend(order = 1, override.aes = list(colour = CONTRAST))) +
  theme_void(base_size = TEXT, base_family = FONT) +
  theme(legend.text = element_text(size = TEXT), legend.position = "bottom", legend.box = "vertical",
        legend.key.size = unit(4.8, "mm"), legend.key.width = unit(8, "mm"), legend.spacing.y = unit(0.4, "mm"),
        plot.margin = margin(4.4, 5, 3, 2, "mm"))

# d: FMR1 edge to its m6A-marked targets (s; circles, filled above the 95th percentile of the control RBPs of the unit)
# and to its bound targets without m6A (u; grey dots) in every unit, by arm and brain CLIP list (human units through
# orthologs); triangles, the striatal mRNA unit refitted on its wild-type mice
ROWS5 <- c("Mouse Output", "Mouse Translation", "Mouse Stability", "Human Output", "Human Translation")
d <- read_panel("fig5d")
d[, `:=`(row = factor(paste(c(mouse = "Mouse", human = "Human")[species], arm), rev(ROWS5)),
         above = abs(b_s) > control_p95, wild = mice == "wild-type")]
jitter <- position_jitter(height = 0.22, width = 0, seed = 2)
dd <- ggplot(d, aes(y = row)) +
  ZERO_X +
  geom_point(aes(x = b_u), colour = GREY, size = 1.4, alpha = 0.8, position = jitter) +
  geom_point(data = d[wild == FALSE], aes(x = b_s, fill = above, colour = above), shape = 21, size = 3.4, stroke = 0.5,
             position = jitter) +
  geom_point(data = d[wild == TRUE], aes(x = b_s, fill = above, colour = above), shape = 24, size = 3.8, stroke = 0.5) +
  facet_wrap(~list, ncol = 1) +
  scale_fill_manual(values = c("TRUE" = FOCAL, "FALSE" = "white"), guide = "none") +
  scale_colour_manual(values = c("TRUE" = EDGE, "FALSE" = GREY), guide = "none") +
  scale_x_continuous(labels = MINUS_LABELS) +
  labs(x = "FMR1 edge (b); grey, without m6A", y = NULL) +
  theme_rows() +
  BOXED +
  theme(panel.spacing = unit(4, "mm"))

# e: percentile of the load edge of the YTHDF module and YTHDF3 among the load edges of the control RBPs (pooled over
# units), and in the striatal mRNA unit of the allelic series refitted on its wild-type mice
p <- read_panel("fig5e")
LEVELS <- c("Mouse Output", "Mouse Translation", "Mouse Stability", "Mouse Stability, wild-type mice", "Human Output",
            "Human Translation", "Human Stability")
p[, `:=`(set = factor(c(YTHDF = "YTHDF module", YTHDF3 = "YTHDF3")[regulator], c("YTHDF module", "YTHDF3")),
         panel = factor(paste(c(mouse = "Mouse", human = "Human")[species], arm), rev(LEVELS)))]
e <- ggplot(p, aes(pctile_signed, panel)) +
  geom_vline(xintercept = 95, linewidth = 0.4, colour = INK, linetype = "22") +
  geom_line(aes(group = panel), colour = EDGE, linewidth = 0.6) +
  geom_point(aes(fill = set), colour = EDGE, shape = 21, size = 5.2, stroke = 0.5) +
  scale_fill_manual(values = c("YTHDF module" = CONTRAST, YTHDF3 = FOCAL), name = NULL) +
  scale_y_discrete(labels = function(x) sub(", wild-type mice", ",\nwild-type mice", x)) +
  scale_x_continuous(breaks = c(0, 50, 95), expand = expansion(add = c(2, 0))) +
  coord_cartesian(xlim = c(0, 100), clip = "off") +
  labs(x = "Percentile of the load edge\namong control RBPs", y = NULL) +
  theme_rows() +
  theme(legend.position = "top", legend.margin = margin(0, 0, -2, 0), axis.text.y = element_text(lineheight = 0.9))

fig <- panel_grid(list(a = a, b = b, c = cc, d = dd, e = e), design = "AABB\nCDDE", widths = c(86, 14, 35, 48),
                  heights = c(52, 82))
save_fig(fig, "fig5", 268)
