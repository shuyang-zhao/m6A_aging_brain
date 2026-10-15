# Shared style of the figures, modelled on the panels of Pang et al. 2026 (RovHer, Nature Genetics): Arial, all text
# black at 12 pt, no panel titles, thin black axes, faceted panels boxed, grey intervals, one palette, the panel layout,
# the save function and the panels that two figures share (sites_panel, upset)
# the figures are drawn at twice the printed size of a Nature double column (366 mm wide), so that 12 pt text prints at
# 6 pt when the figure is set at 183 mm; every size below is at that drawn size
# run from the repository root; the figure scripts read source_data/ and write output/figures/<name>.png and .pdf
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

FONT <- "Arial"
stopifnot("Arial is not installed" = FONT %in% systemfonts::system_fonts()$family)
WIDTH <- 366
MAX_HEIGHT <- 340
DPI <- 300
TEXT <- 12 # every text of every figure, in pt
PT <- TEXT / .pt # geom_text size of 12 pt text

# palette of the reference: deep teal (m6A, protein and the focal results), brick red (protein given mRNA, a second
# family), mustard (a contrasting series) and light blue (further series), each with a light tone for areas; open
# grey points and grey sets, the comparisons and non-significant results. Heatmaps run from teal (lower in old age)
# through white to brick red (higher)
FOCAL <- "#087C9A"
FOCAL_DARK <- "#05566B"
FOCAL_LIGHT <- "#A6CFDA"
RED <- "#C3535E"
RED_LIGHT <- "#EDB8BC"
CONTRAST <- "#EAC64B"
CONTRAST_LIGHT <- "#F6E3A4"
BLUE <- "#8DB6D6"
BLUE_LIGHT <- "#D2E3F0"
DIVERGING <- c("#087C9A", "#3B96AE", "#7EBAC9", "#C3DFE6", "#F7F7F7", "#F3CFCF", "#E3999E", "#D06E76", "#B5404C")
INK <- "#000000" # all text, axes and ticks
EDGE <- "#4D4D4D" # outlines of points and bars, arrows and connectors
CI <- "#8C8C8C" # confidence intervals (EDGE where they cross bars)
GREY <- "#BBBBBB" # comparison sets
NS <- "#DADADA" # not significant
LIGHT <- "#EBEBEB"
NA_FILL <- "#F0F0F0" # not quantified
# each molecular layer keeps one hue across the figures (solid for points, lines and bars, light for areas)
LAYER <- c(Protein = FOCAL, "Protein given mRNA" = RED, "Bulk RNA" = BLUE, "Single-nucleus RNA" = CONTRAST)
LAYER_LIGHT <- c(Protein = FOCAL_LIGHT, "Protein given mRNA" = RED_LIGHT, "Bulk RNA" = BLUE_LIGHT,
                 "Single-nucleus RNA" = CONTRAST_LIGHT)

# the four contexts of the load models with the layer of each, and the five targets of the gene-level models
CONTEXTS <- c("Mouse protein", "Mouse protein given mRNA", "Human protein", "Human striatum mRNA")
KEYS <- setNames(c("Protein", "Protein given mRNA", "Protein", "Bulk RNA"), CONTEXTS)
TARGETS <- c("Mouse protein", "Mouse protein given mRNA", "Human protein", "Human mRNA", "Human neuronal mRNA")
# outline and fill of a point by its layer: filled in the layer's colour where P < 0.05, open and grey otherwise
p_colours <- function(p, key) {
  list(edge = fifelse(p < 0.05, EDGE, GREY), inside = fifelse(p < 0.05, LAYER[key], "white"))
}

theme_fig <- function() {
  theme_classic(base_size = TEXT, base_family = FONT) +
    theme(
      geom = element_geom(ink = INK, family = FONT, fontsize = TEXT, linewidth = 0.5),
      text = element_text(colour = INK, size = TEXT),
      line = element_line(colour = INK, linewidth = 0.45),
      axis.line = element_line(colour = INK, linewidth = 0.45),
      axis.ticks = element_line(colour = INK, linewidth = 0.45),
      axis.ticks.length = unit(1.4, "mm"),
      axis.text = element_text(size = TEXT, colour = INK),
      axis.title = element_text(size = TEXT, colour = INK),
      strip.background = element_blank(),
      strip.text = element_text(size = TEXT, colour = INK, margin = margin(0, 0, 2, 0, "mm")),
      strip.clip = "off",
      legend.text = element_text(size = TEXT, colour = INK),
      legend.title = element_text(size = TEXT, colour = INK),
      legend.key.size = unit(4.8, "mm"),
      legend.key.spacing.y = unit(0.8, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(3, "mm"),
      legend.background = element_blank(),
      plot.title = element_blank(),
      plot.margin = margin(4.4, 5, 3, 2, "mm"),
      panel.spacing = unit(6, "mm")
    )
}
# faceted panels: each facet in a thin black box
BOXED <- theme(panel.border = element_rect(colour = INK, fill = NA, linewidth = 0.45), axis.line = element_blank())

# heatmaps (tables): no axes
theme_heat <- function() {
  theme_fig() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
          panel.spacing = unit(1.6, "mm"))
}
# estimates by category, rows labelled without a y axis line
theme_rows <- function() {
  theme_fig() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
}
# dot matrices in place of heatmaps: a faint grid, no axes; the colour of a dot is its value, its size the evidence
theme_matrix <- function() {
  theme_fig() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
          panel.grid.major = element_line(colour = "#EDEDED", linewidth = 0.5), panel.spacing = unit(3, "mm"))
}
MATRIX_EDGE <- "#9E9E9E" # outline of the dots, so that values near zero stay visible
# nested patchworks keep no margin of their own
FLUSH <- plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
# panels are separate elements (neighbours do not align their axes) in cells without margins; each is lettered in
# bold lower case at the top left corner of its cell
panel_grid <- function(panels, ...) {
  tag <- theme(
    plot.tag = element_text(family = FONT, size = TEXT, face = "bold", colour = INK, hjust = 0, vjust = 1,
                            margin = margin(t = 2, unit = "mm")),
    plot.tag.location = "plot", plot.tag.position = c(0, 1), plot.margin = margin(0, 0, 0, 0)
  )
  wrap_plots(Map(function(p, letter) wrap_elements(full = p) + labs(tag = letter) + tag, panels, names(panels)), ...) +
    FLUSH
}

MINUS <- function(x) gsub("-", "−", x, fixed = TRUE)
# number with its sign (and a true minus) and two decimals, for annotations of effects
signed <- function(x) MINUS(sprintf("%+.2f", x))
MINUS_LABELS <- scales::label_number(style_negative = "minus", drop0trailing = TRUE)
cell_number <- function(x) MINUS(sprintf("%.1f", round(x, 1) + 0))

# P or q as in the text, two significant digits and bounded below at 0.0001
fmt_p <- function(p, label = "P") {
  exact <- formatC(signif(p, 2), format = "fg", flag = "#", digits = 2)
  ifelse(p < 1e-4, paste(label, "< 0.0001"), paste(label, "=", exact))
}

# diverging fill scale of a heatmap, squished at the limit (the end labels say so), its colour bar first in the legend
scale_fill_div <- function(limit, name) {
  ends <- function(x) {
    out <- MINUS_LABELS(x)
    out[x == -limit] <- paste0("≤", out[x == -limit])
    out[x == limit] <- paste0("≥", out[x == limit])
    out
  }
  scale_fill_gradientn(
    colours = DIVERGING, limits = c(-limit, limit), oob = scales::squish, na.value = NA_FILL, name = name,
    breaks = c(-limit, 0, limit), labels = ends,
    guide = guide_colourbar(
      theme = theme(legend.key.width = unit(26, "mm"), legend.key.height = unit(3.2, "mm")),
      title.position = "top", title.hjust = 0.5, ticks = FALSE, frame.colour = NA, order = 1
    )
  )
}

# reference lines dashed, data lines solid
ZERO_X <- geom_vline(xintercept = 0, linewidth = 0.4, colour = INK, linetype = "22")
ZERO_Y <- geom_hline(yintercept = 0, linewidth = 0.4, colour = INK, linetype = "22")

read_panel <- function(name) fread(file.path("source_data", paste0(name, ".tsv")))

# schematics: rounded boxes (ggforce) with a bold first line, on a 0-100 frame
boxes <- function(b) {
  corners <- b[, .(x = c(x0, x1, x1, x0), y = c(y0, y0, y1, y1)), by = .(id, fill)]
  list(
    ggforce::geom_shape(data = corners, aes(x, y, group = id, fill = I(fill)), radius = unit(2.4, "mm")),
    geom_text(data = b, aes((x0 + x1) / 2, y1 - 3.5, label = head), size = PT, fontface = "bold", colour = INK,
              vjust = 1),
    geom_text(data = b, aes((x0 + x1) / 2, (y0 + y1) / 2 - 3, label = body), size = PT, colour = INK,
              lineheight = 0.95)
  )
}
frame <- function() {
  list(coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE), theme_void(base_size = TEXT, base_family = FONT),
       theme(plot.margin = margin(4.4, 5, 3, 2, "mm")))
}

# age effect of genes by number of m6A sites (Figures 2b and 3a): protein, mean over units with a 95% interval across
# units (number, the load coefficient over cohorts as cited in the text); protein given mRNA, one line per unit with its
# 95% CI; `open` names the units drawn with a white fill, `per` is added to the unit of the coefficient
sites_panel <- function(name, species, axes, shapes = c(21, 22, 24), open = NULL, per = "") {
  panels <- paste(species, c("protein", "protein\ngiven mRNA"))
  mean_ci <- read_panel(name)[, .(sites, effect, low = ci_low, high = ci_high, line = "Mean over units")]
  units <- read_panel(paste0(name, "_given_mrna"))[, .(sites, effect, low = ci_low, high = ci_high, line = unit)]
  d <- rbind(mean_ci[, panel := panels[1]], units[, panel := panels[2]])
  line_order <- c("Mean over units", sort(unique(units$line)))
  d[, `:=`(sites = factor(sites, c("0", "1", "2-3", ">=4"), c("0", "1", "2–3", "≥4")),
           panel = factor(panel, panels), line = factor(line, line_order),
           key = fifelse(grepl("given", panel), "Protein given mRNA", "Protein"))]
  d[, inside := fifelse(line %in% open, "white", LAYER[key])]
  coef <- read_panel("fig2d")[panel == panels[1] & step == "programmes"]
  coef[, `:=`(panel = factor(panel, panels), label = paste0(signed(effect), "%", per, "\nper log unit\n", fmt_p(p)))]
  dodge <- position_dodge(0.3)
  key_fill <- fifelse(line_order %in% open, "white", GREY)
  ggplot(d, aes(sites, effect, group = line, shape = line)) +
    ZERO_Y +
    geom_linerange(aes(ymin = low, ymax = high), linewidth = 0.6, colour = CI, position = dodge) +
    geom_line(aes(colour = I(LAYER[key])), linewidth = 1.1, position = dodge) +
    geom_point(aes(fill = I(inside)), colour = EDGE, size = 3.6, stroke = 0.5, position = dodge) +
    geom_text(data = coef, aes(x = Inf, y = Inf, label = label), hjust = 1.06, vjust = 1.2, size = PT,
              lineheight = 0.95, colour = INK, inherit.aes = FALSE) +
    facet_wrap(~panel, nrow = 1, scales = "free_y") +
    scale_shape_manual(values = shapes, name = NULL,
                       guide = guide_legend(override.aes = list(colour = EDGE, fill = key_fill))) +
    scale_y_continuous(labels = MINUS_LABELS, expand = expansion(mult = c(0.08, 0.1))) +
    axes +
    theme_fig() +
    BOXED +
    theme(legend.position = "bottom", strip.text = element_text(lineheight = 0.95), panel.spacing = unit(5, "mm"))
}

# UpSet panel (Figures S2b and S6b): bars of the genes in each combination of sets (rows of u, largest first; `combo`,
# the position of the combination, and one logical column per set, named in `sets` with the label of the set as value)
# above a matrix of dots that marks the sets of each combination; `...` are the layers of the bars and dot_colour(d)
# gives the colour of the marked dots from the long table d
upset <- function(u, sets, ..., dot_colour, dot_size, heights) {
  dots <- melt(u, id.vars = setdiff(names(u), names(sets)), measure.vars = names(sets), variable.name = "set",
               value.name = "member")
  dots[, set := factor(sets[as.character(set)], rev(sets))]
  bars <- ggplot(u, aes(combo, genes)) +
    theme_fig() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(4.4, 5, 0, 2, "mm")) +
    labs(x = NULL, y = "Genes") +
    list(...)
  dot_matrix <- ggplot(dots, aes(combo, set)) +
    geom_line(data = dots[member == TRUE], aes(group = combo), colour = EDGE, linewidth = 0.8) +
    geom_point(aes(colour = I(fifelse(member, dot_colour(dots), "#E3E3E3"))), size = dot_size) +
    labs(x = NULL, y = NULL) +
    theme_fig() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
          plot.margin = margin(0, 5, 3, 2, "mm"))
  (bars / dot_matrix + plot_layout(heights = heights)) + FLUSH
}

# modification names (m6A, m5C, m7G) with the position as a superscript wherever a figure writes them: such texts
# become plotmath, so that the digit is set in Arial (the font has no superscript digits); several lines are stacked
# with atop
MARK <- "m[0-9][ACGU]"
superscript_marks <- function(g) {
  as_math <- function(s) {
    lines <- lapply(strsplit(s, "\n", fixed = TRUE)[[1]], function(line) {
      found <- gregexpr(MARK, line)
      parts <- regmatches(line, found, invert = TRUE)[[1]]
      marks <- regmatches(line, found)[[1]]
      pieces <- list(parts[1])
      for (i in seq_along(marks)) {
        mark <- bquote(m^.(substr(marks[i], 2, 2)) * .(as.name(substr(marks[i], 3, 3))))
        pieces <- c(pieces, mark, parts[i + 1])
      }
      as.call(c(as.name("paste"), pieces))
    })
    Reduce(function(a, b) call("atop", a, b), lines)
  }
  # the other labels of the same text become plotmath too (one baseline), unless they have several lines
  if (inherits(g, "text") && is.character(g$label) && any(grepl(MARK, g$label))) {
    convert <- function(s) !is.na(s) && (grepl(MARK, s) || (nzchar(s) && !grepl("\n", s)))
    g$label <- as.expression(lapply(g$label, function(s) if (convert(s)) as_math(s) else s))
  }
  if (inherits(g, "gTree")) g$children[] <- lapply(g$children, superscript_marks)
  if (inherits(g, "gtable")) g$grobs <- lapply(g$grobs, superscript_marks)
  g
}

# the figure is drawn once as a 300-dpi image (600 dpi when set at 183 mm) and the PDF holds that image, so that no
# viewer can substitute or stack the glyphs of a font
save_fig <- function(fig, name, height) {
  stopifnot("figure is too tall" = height <= MAX_HEIGHT)
  file <- file.path("output/figures", name)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  png(paste0(file, ".png"), width = WIDTH, height = height, units = "mm", res = DPI, type = "cairo", bg = "white")
  grid::grid.draw(superscript_marks(patchworkGrob(fig)))
  invisible(dev.off())
  pdf(paste0(file, ".pdf"), width = WIDTH / 25.4, height = height / 25.4)
  grid::grid.raster(png::readPNG(paste0(file, ".png")), interpolate = FALSE)
  invisible(dev.off())
  message(file, ".pdf: ", WIDTH, " x ", height, " mm at ", DPI, " dpi")
}
