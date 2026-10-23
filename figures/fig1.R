# Figure 1: m6Brain, the cohorts and the design of the clock
# a framework schematic, b cohorts (source_data/fig1b.tsv), c design diagram; writes output/figures/fig1.png and .pdf
source("figures/style.R")

LINE <- 4.9 # mm per line of 12 pt text in the schematics

# rounded box with a bold head line and body lines, text block centred in the box; coordinates in mm
node <- function(x0, y0, x1, y1, fill, head, body = "", head_y = NA, body_y = NA) {
  n <- if (nzchar(body)) lengths(strsplit(body, "\n", fixed = TRUE)) else 0L
  cy <- (y0 + y1) / 2
  data.table(x0, y0, x1, y1, fill, head, body, id = paste(x0, y0), cx = (x0 + x1) / 2,
             head_y = fifelse(is.na(head_y), cy + n * LINE / 2, head_y),
             body_y = fifelse(is.na(body_y), cy - LINE / 2, body_y))
}
shapes <- function(b) {
  corners <- b[, .(x = c(x0, x1, x1, x0), y = c(y0, y0, y1, y1)), by = .(id, fill)]
  list(
    ggforce::geom_shape(data = corners, aes(x, y, group = id, fill = I(fill)), radius = unit(2.4, "mm")),
    geom_text(data = b[head != "m6A"], aes(cx, head_y, label = head), size = PT, fontface = "bold", colour = INK),
    # plotmath so that the digit is set in Arial and bold
    geom_text(data = b[head == "m6A"], aes(cx, head_y, label = "bold(m^6*A)"), parse = TRUE, size = PT, colour = INK),
    geom_text(data = b[nzchar(body)], aes(cx, body_y, label = body), size = PT, colour = INK, lineheight = 0.95)
  )
}
arrows <- function(x, y, xend, yend) {
  geom_segment(data = data.table(x, y, xend, yend), aes(x, y, xend = xend, yend = yend), colour = EDGE,
               linewidth = 0.7, arrow = arrow(length = unit(2.4, "mm"), type = "closed"))
}
canvas <- function(w, h) {
  list(coord_cartesian(xlim = c(0, w), ylim = c(0, h), expand = FALSE),
       theme_void(base_size = TEXT, base_family = FONT), theme(plot.margin = margin(0, 0, 0, 0)))
}

# a: framework, 366 x 92 mm
WA <- 366
HA <- 92
feat <- rbind(
  node(50, 65.9, 124, 78, BLUE_LIGHT, "Orthrus", "RNA foundation model"),
  node(50, 51.3, 124, 63.4, BLUE_LIGHT, "RiboNN", "translation model"),
  node(50, 36.7, 124, 48.8, FOCAL_LIGHT, "m6A", "sites and reader binding"),
  node(50, 22.1, 124, 34.2, BLUE_LIGHT, "Structure", "transcript structure features"),
  node(50, 7.5, 124, 19.6, BLUE_LIGHT, "Context", "expression, cell type, programs")
)
steps <- rbind(
  node(4, 7.5, 38, 78, LIGHT, "Transcripts", "of each\ngene"),
  feat,
  node(136, 7.5, 180, 78, CONTRAST_LIGHT, "TabPFN", "tabular\nfoundation\nmodel"),
  node(192, 28, 238, 58, FOCAL_LIGHT, "Prediction", "age effect\nof each gene"),
  node(250, 28, 296, 58, FOCAL_LIGHT, "m6Brain clock", "elastic net\nof age, guided\nby the prediction"),
  node(308, 46, 362, 72, RED_LIGHT, "Age", "r and MAE in\nheld-out and\nexternal cohorts"),
  node(308, 14, 362, 40, RED_LIGHT, "Disease", "AD, PD, ALS, HD\nin years")
)
brackets <- data.table(x = c(50, 250), xend = c(238, 362), label = c("m6Brain: gene-level model", "Sample-level clock"))
a <- ggplot() +
  shapes(steps) +
  arrows(38, feat[, (y0 + y1) / 2], 50, feat[, (y0 + y1) / 2]) +
  arrows(124, feat[, (y0 + y1) / 2], 136, feat[, (y0 + y1) / 2]) +
  arrows(c(180, 238, 296, 296), rep(43, 4), c(192, 250, 308, 308), c(43, 43, 59, 27)) +
  geom_segment(data = brackets, aes(x, 83, xend = xend, yend = 83), colour = INK, linewidth = 0.45) +
  geom_text(data = brackets, aes((x + xend) / 2, 88, label = label), size = PT, colour = INK) +
  canvas(WA, HA)

# b: age range (segment), donors of the largest unit (dot size at the midpoint) and analysis units (right) of every
# dataset, grouped by set and layer; mouse (months) and human (years)
d <- read_panel("fig1b")
GROUPS <- c("Mouse protein", "Mouse RNA", "Human protein", "Human bulk RNA", "Human single-nucleus RNA")
SETS <- c(discovery = "Discovery", external = "External", validation = "Validation")
d[, `:=`(species = fifelse(grepl("^Mouse", group), "Mouse", "Human"), group = factor(group, GROUPS),
         set = factor(set, names(SETS)))]
setorder(d, group, set, age_min, age_max)
d[, head := paste0(SETS[as.character(set)], ": ", tolower(substr(group, 1, 1)), substring(group, 2))]
d[, head := factor(head, unique(head))]
d[, col := fcase(species == "Human" & layer == "Protein", RED, species == "Human", RED_LIGHT,
                 set == "external", CONTRAST, layer == "Protein", FOCAL, default = FOCAL_LIGHT)]
d[, k := seq_len(.N), by = head]
offset <- d[, .(n = .N, species = species[1]), by = head][order(head)]
offset[, start := cumsum(c(0, head(n + 1.4, -1))), by = species]
d <- merge(d, offset[, .(head, start)], by = "head")[, y := start + k]
YMAX <- max(d$y) + 0.7
AGE <- list(Mouse = list(lim = c(0, 35), breaks = seq(0, 30, 10), name = "Age (months)"),
            Human = list(lim = c(14, 107), breaks = seq(20, 100, 20), name = "Age (years)"))
HEAD_X <- c(Mouse = -19.5, Human = -57) # heading start in age units, left of the dataset names
overview <- function(sp) {
  age <- AGE[[sp]]
  rows <- d[species == sp]
  heads <- offset[species == sp][, y := start]
  ggplot(rows) +
    geom_segment(aes(x = age_min, xend = age_max, y = y, yend = y, colour = I(col)), linewidth = 0.9,
                 lineend = "round") +
    geom_point(aes((age_min + age_max) / 2, y, size = donors, fill = I(col)), shape = 21, colour = EDGE,
               stroke = 0.5) +
    geom_text(data = heads, aes(HEAD_X[[sp]], y, label = head), hjust = 0, fontface = "bold", size = PT) +
    annotate("text", x = age$lim[2], y = heads$y[1], label = "Units", hjust = -0.2, size = PT) +
    scale_y_reverse(breaks = rows$y, labels = rows$dataset, limits = c(YMAX, -0.3), expand = expansion(0),
                    sec.axis = dup_axis(breaks = rows$y, labels = rows$n_units, name = NULL)) +
    scale_x_continuous(breaks = age$breaks, expand = expansion(0)) +
    scale_size(range = c(1.4, 5.4), limits = c(0, 300), breaks = c(10, 50, 100, 300),
               name = "Donors in the largest unit",
               guide = guide_legend(nrow = 1, override.aes = list(fill = GREY, colour = EDGE))) +
    coord_cartesian(xlim = age$lim, clip = "off") +
    labs(x = age$name, y = NULL) +
    theme_rows() +
    theme(plot.margin = margin(8, 9, 3, 2, "mm"))
}
b <- wrap_plots(overview("Mouse"), plot_spacer(), overview("Human"), nrow = 1, widths = c(116, 8, 116),
                guides = "collect") &
  theme(legend.position = "bottom")

# c: design of the clock: the gene-level model and the clock refitted without each cohort, and external cohorts scored
# with the clock of all cohorts; the gene-level model is fitted on held-out chromosomes
WC <- 126
HC <- 170
squares <- function(x0, x1, y0, y1, n, fill) {
  w <- (x1 - x0 - (n - 1) * 1.6) / n
  data.table(xmin = x0 + (seq_len(n) - 1) * (w + 1.6), xmax = x0 + (seq_len(n) - 1) * (w + 1.6) + w,
             ymin = y0, ymax = y1, fill = fill)
}
n_held <- d[group == "Mouse protein" & set == "discovery", .N]
n_external <- d[group == "Mouse protein" & set == "external", .N]
cohorts <- rbind(squares(4, 62, 143, 151, n_held, c(rep(FOCAL_LIGHT, 3), FOCAL, rep(FOCAL_LIGHT, n_held - 4))),
                 squares(64, 122, 143, 151, n_held, FOCAL_LIGHT),
                 squares(29, 37, 28, 35, 1, FOCAL),
                 squares(69, 117, 28, 35, n_external, CONTRAST))
chrom <- squares(14, 112, 99, 107, 11, c(rep(BLUE_LIGHT, 6), NA, rep(BLUE_LIGHT, 4)))
design <- rbind(
  node(4, 84, 122, 120, LIGHT, "Held-out chromosomes", "a gene is predicted by a model that\nhas not seen its chromosome",
       head_y = 115.5, body_y = 91),
  node(4, 56, 62, 76, FOCAL_LIGHT, "Clock", "refitted without\nthe held-out cohort"),
  node(64, 56, 122, 76, FOCAL_LIGHT, "Clock", "of all cohorts"),
  node(4, 24, 62, 50, FOCAL_LIGHT, "Predicted age", "held-out cohort", head_y = 45.5, body_y = 40.5),
  node(64, 24, 122, 50, CONTRAST_LIGHT, "Predicted age", "external cohorts", head_y = 45.5, body_y = 40.5),
  node(4, 2, 122, 14, RED_LIGHT, "Pearson r with age per cohort")
)
cc <- ggplot() +
  shapes(design) +
  geom_rect(data = cohorts, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = I(fill)),
            colour = EDGE, linewidth = 0.4) +
  geom_rect(data = chrom[!is.na(fill)], aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = I(fill)),
            colour = EDGE, linewidth = 0.4) +
  geom_rect(data = chrom[is.na(fill)], aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white",
            colour = EDGE, linewidth = 0.8, linetype = "22") +
  geom_text(data = data.table(x = c(33, 93), label = c("Held-out cohort", "External cohorts")),
            aes(x, 162, label = label), size = PT, fontface = "bold") +
  geom_text(data = data.table(x = c(33, 93), label = c("one cohort held out\nin turn", "all cohorts of the set")),
            aes(x, 134.5, label = label), size = PT, lineheight = 0.95) +
  arrows(c(33, 93, 33, 93, 33, 93, 33, 93), c(128, 128, 84, 84, 56, 56, 24, 24),
         c(33, 93, 33, 93, 33, 93, 33, 93), c(121, 121, 77, 77, 51, 51, 15, 15)) +
  canvas(WC, HC)

fig <- panel_grid(list(a = a, b = b, c = cc), design = "AAA\nBBC", widths = c(120, 120, WC), heights = c(HA, 172))
save_fig(fig, "fig1", HA + 172)
