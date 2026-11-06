# Figure 1: m6Brain, the cohorts and the design of the clock
# a schematic; b cohorts (source_data/fig1b, fig1b_donors); c design of the held-out predictions; d protein against mRNA
# of the same animals (source_data/fig1d); writes output/figures/fig1
source("figures/style.R")

arrow_to <- function(x, xend, y, yend = y) {
  annotate("segment", x = x, xend = xend, y = y, yend = yend, linewidth = lw(0.6), colour = EDGE,
           arrow = grid::arrow(length = unit(1, "mm"), type = "closed"))
}
FEATURES <- c("Orthrus RNA embedding", "RiboNN translation", "m6A sites, reader binding", "Transcript structure",
              "Expression, cell type, programs")
# a toy matrix of the feature blocks: a shade per gene and feature, the m6A block in the colour of m6Brain
toy <- CJ(row = seq_along(FEATURES), column = 1:9)
toy[, shade := (sin(row * 2.1 + column * 1.7) + 1) / 2]
toy[, fill := mapply(function(s, r) tint(if (r == 3) M6BRAIN else "#6E6E6E", 0.25 + 0.65 * s), shade, row)]
toy[, `:=`(x = 40 + column * 2.1, y = 41.5 - row * 6.6)]
# predictions of the toy genes, ranked, and the weights of the clock that they guide
predicted <- data.table(i = 1:12, z = c(-3.2, -2.5, -2, -1.4, -1, -0.5, -0.2, 0.3, 0.8, 1.3, 2.1, 2.9))
weights <- data.table(i = 1:12, w = c(-1.8, -1.1, 0, -0.5, 0, 0, 0, 0, 0.3, 0, 0.9, 1.6))
SITES_X <- c(5, 12, 17, 22, 27, 30, 32)

a <- canvas(183, 50) +
  annotate("segment", x = 8, xend = 136, y = 49, yend = 49, linewidth = lw(0.5)) +
  annotate("segment", x = 141, xend = 181, y = 49, yend = 49, linewidth = lw(0.5)) +
  note(72, 46.8, "m6Brain: gene-level model") +
  note(161, 46.8, "Sample-level clock") +
  # transcript
  note(18, 40.4, "m6A sites and readers") +
  annotate("rect", xmin = 2, xmax = 34, ymin = 21.6, ymax = 23, fill = tint(M6BRAIN, 0.7)) +
  annotate("rect", xmin = 9, xmax = 25, ymin = 20.4, ymax = 24.2, fill = tint(M6BRAIN, 0.55)) +
  annotate("segment", x = SITES_X, xend = SITES_X, y = 23, yend = 30, colour = M6BRAIN, linewidth = lw(0.5)) +
  annotate("point", x = SITES_X, y = 30.4, colour = M6BRAIN, size = 0.9) +
  annotate("point", x = c(17, 27), y = 33, shape = 21, fill = MOUSE, colour = EDGE, size = 2.1, stroke = st(0.3)) +
  note(17, 22.3, "CDS") +
  note(5.5, 18.6, "5′ UTR") +
  note(29.5, 18.6, "3′ UTR") +
  note(18, 11, "Transcript of each gene") +
  arrow_to(35, 40, 22.3) +
  # features
  geom_tile(data = toy, aes(x, y, fill = I(fill)), width = 2, height = 5.2, inherit.aes = FALSE) +
  note(rep(61.8, 5), 41.5 - seq_along(FEATURES) * 6.6, FEATURES, hjust = 0) +
  # TabPFN
  arrow_to(96, 100, 22.3) +
  annotate("rect", xmin = 100.5, xmax = 112, ymin = 7, ymax = 37.5, fill = tint(M6BRAIN, 0.8)) +
  note(106.2, 22.3, "TabPFN", angle = 90) +
  arrow_to(112.5, 116.5, 22.3) +
  # predictions
  note(127, 40.4, "Predicted age effect\nof each gene", lineheight = 0.95) +
  geom_rect(data = predicted, aes(xmin = 116 + i * 1.6, xmax = 117.2 + i * 1.6, ymin = 22.3, ymax = 22.3 + z * 3.2,
                                  fill = I(fifelse(z < 0, tint(MOUSE, 0.15), tint(HUMAN, 0.15)))), inherit.aes = FALSE) +
  annotate("segment", x = 117, xend = 137, y = 22.3, yend = 22.3, linewidth = lw(0.5)) +
  # clock
  arrow_to(138, 142.5, 22.3) +
  note(153.5, 40.4, "Sets the penalty of\nan elastic net of age", lineheight = 0.95) +
  geom_rect(data = weights, aes(xmin = 142 + i * 1.6, xmax = 143.2 + i * 1.6, ymin = 22.3, ymax = 22.3 + w * 4.2,
                                fill = I(M6BRAIN)), inherit.aes = FALSE) +
  annotate("segment", x = 143, xend = 163, y = 22.3, yend = 22.3, linewidth = lw(0.5)) +
  note(153.5, 11, "Clock weights") +
  arrow_to(164, 168, 25, 29) +
  arrow_to(164, 168, 19.6, 15.6) +
  annotate("rect", xmin = 168.2, xmax = 181, ymin = 25, ymax = 37.5, fill = tint(GREY, 0.6)) +
  annotate("rect", xmin = 168.2, xmax = 181, ymin = 7, ymax = 19.5, fill = tint(HUMAN, 0.7)) +
  note(174.6, 31.2, "Age\n25 cohorts", lineheight = 0.95) +
  note(174.6, 13.2, "Disease\nin years", lineheight = 0.95)

# b: the cohorts: age of every donor of the proteomes (ticks), age range of the transcriptomes (bars), donors of the
# largest unit and number of analysis units
cohorts <- read_panel("fig1b")
ages <- read_panel("fig1b_donors")
BLOCKS <- c("Mouse protein|discovery" = "Discovery", "Mouse protein|external" = "External", "Mouse RNA|discovery" = "RNA",
            "Human protein|validation" = "Protein", "Human bulk RNA|validation" = "Bulk RNA",
            "Human single-nucleus RNA|validation" = "Nuclei")
HUES <- c(Discovery = MOUSE, External = EXTERNAL, RNA = tint(MOUSE, 0.6), Protein = HUMAN, "Bulk RNA" = tint(HUMAN, 0.5),
          Nuclei = tint(HUMAN, 0.7))
cohorts[, block := factor(BLOCKS[paste(group, set, sep = "|")], BLOCKS)]
atlas <- function(xlab, limit, breaks, shown) {
  d <- cohorts[block %in% shown]
  d[, block := factor(block, shown)]
  setorder(d, block, -donors)
  d[, row := factor(paste(block, dataset, n_units), rev(paste(block, dataset, n_units)))]
  ticks <- ages[d[block %in% shown[shown %in% c("Discovery", "External", "Protein")], .(dataset, row, block)], on = "dataset", nomatch = NULL]
  numbers <- limit[2] + diff(limit) * c(0.14, 0.42)
  head <- data.table(block = factor(shown[1], shown), x = numbers, text = c("Donors", "Units"))
  ggplot(d) +
    geom_segment(aes(age_min, row, xend = age_max, yend = row, colour = block), linewidth = lw(2.4)) +
    geom_tile(data = ticks, aes(age, row), fill = "black", width = diff(limit) / 170, height = 0.62) +
    geom_text(aes(numbers[1], row, label = donors), hjust = 1, size = pt_size(6)) +
    geom_text(aes(numbers[2], row, label = n_units), hjust = 1, size = pt_size(6)) +
    geom_text(data = head, aes(x, Inf, label = text), hjust = 1, vjust = -0.5, size = pt_size(6)) +
    ggh4x::facet_grid2(
      block ~ ., scales = "free_y", space = "free_y", switch = "y",
      strip = ggh4x::strip_themed(background_y = ggh4x::elem_list_rect(fill = tint(HUES[shown], 0.3), colour = tint(HUES[shown], 0.3)))
    ) +
    scale_colour_manual(values = sapply(HUES, tint, 0.6), guide = "none") +
    scale_x_continuous(breaks = breaks, expand = expansion(mult = c(0.02, 0))) +
    scale_y_discrete(labels = setNames(d$dataset, d$row), expand = expansion(add = 0.6)) +
    coord_cartesian(xlim = c(limit[1], numbers[2]), clip = "off") +
    labs(x = xlab, y = NULL) +
    theme_pub() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), panel.spacing = unit(1.2, "mm"),
          strip.placement = "outside", strip.text.y.left = element_text(angle = 90, size = 6), strip.clip = "off",
          axis.title.x = element_text(hjust = 0.35), plot.margin = margin(4, 1.5, 1.5, 1.5, "mm"))
}
b <- atlas("Age (months)", c(0, 34), c(0, 10, 20, 30), c("Discovery", "External", "RNA")) +
  atlas("Age (years)", c(15, 105), c(25, 50, 75, 100), c("Protein", "Bulk RNA", "Nuclei"))

# c: held-out predictions: a gene is predicted by a model without its chromosome, a cohort by a clock without it, and
# the external cohorts by the clock of all discovery cohorts
grid <- CJ(row = 1:9, column = 1:13)
grid[, kind := fifelse(row >= 7, "External cohorts, scored only", fifelse(row == 3, "Held-out cohort", "Training cohorts"))]
grid[column == 9 & row < 7, kind := "Held-out chromosome"]
KINDS <- c("Training cohorts" = tint(MOUSE, 0.6), "Held-out cohort" = MOUSE, "Held-out chromosome" = "white",
           "External cohorts, scored only" = tint(EXTERNAL, 0.3))
c_ <- ggplot(grid, aes(column, -row, fill = kind)) +
  geom_tile(colour = "white", linewidth = lw(0.5), width = 0.94, height = 0.9) +
  annotate("rect", xmin = 8.5, xmax = 9.5, ymin = -6.5, ymax = -0.5, fill = NA, colour = EDGE, linetype = "22", linewidth = lw(0.5)) +
  note(7, 0.1, "Chromosomes (genes)", vjust = 0) +
  note(0.2, -3.5, "Cohorts (samples)", angle = 90, vjust = 0) +
  key_rows(names(KINDS), KINDS, x = 1, y = -10.4, step = 1.25, width = 0.7, square = TRUE) +
  annotate("rect", xmin = 1, xmax = 1.7, ymin = -13.2, ymax = -12.6, fill = NA, colour = EDGE, linetype = "22", linewidth = lw(0.4)) +
  scale_fill_manual(values = KINDS, guide = "none") +
  coord_cartesian(xlim = c(-0.2, 14), ylim = c(-15, 1.1), expand = FALSE, clip = "off") +
  theme_void_pub()

# d: protein against mRNA age t in the mouse cerebrum, the same animals
d <- read_panel("fig1d")
rho <- d[, cor(protein_t, mrna_t, method = "spearman")]
n <- nrow(d)
shown <- d[, .(x = quantile(mrna_t, c(0.001, 0.999)), y = quantile(protein_t, c(0.001, 0.999)))]
d <- by_density(d, "mrna_t", "protein_t")
d <- ggplot(d, aes(mrna_t, protein_t, colour = .density)) +
  ref_line(yintercept = 0) +
  ref_line(xintercept = 0) +
  rasterise_pts(geom_point(size = 0.4, stroke = 0)) +
  note_corner(lab_stats(rho, stat = "rho", n = n), parse = TRUE) +
  scale_density(MOUSE) +
  scale_x_continuous(labels = MINUS_LABELS) +
  scale_y_continuous(labels = MINUS_LABELS) +
  coord_cartesian(xlim = shown$x, ylim = shown$y) +
  labs(x = "mRNA age t", y = "Protein age t") +
  theme_pub()

fig <- fig_panels(a, b, c_, d, design = "AAAAAAAAAAAA\nBBBBBBBBCCCC\nBBBBBBBBDDDD", heights = c(50, 60, 56))
save_fig(fig, "fig1", 168)
