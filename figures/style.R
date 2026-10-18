# Shared style of the figures, modelled on the panels of Pang et al. 2026 (RovHer, Nature Genetics): Arial, all text
# black at 12 pt, no panel titles, thin black axes, faceted panels boxed, grey intervals, one palette, the panel layout,
# the save function and the helpers that several figures share
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

# palette of the reference: deep teal (m6A, protein and the focal results), brick red (a second family), mustard (a
# contrasting series) and light blue (further series), each with a light tone for areas; open grey points and grey
# sets, the comparisons and non-significant results. Heatmaps run from teal (lower in old age) through white to brick
# red (higher)
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
# the five targets of the gene-level models
TARGETS <- c("Mouse protein", "Mouse mRNA-adjusted protein", "Human protein", "Human mRNA", "Human neuronal mRNA")
SPECIES_COL <- c(mouse = FOCAL, human = RED)
SPECIES_NAME <- c(mouse = "Mouse", human = "Human")

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
# bold lower case letter at the top left corner of a panel or cell
PANEL_TAG <- theme(
  plot.tag = element_text(family = FONT, size = TEXT, face = "bold", colour = INK, hjust = 0, vjust = 1,
                          margin = margin(t = 2, unit = "mm")),
  plot.tag.location = "plot", plot.tag.position = c(0, 1)
)
# panels are separate elements (neighbours do not align their axes) in cells without margins, each lettered
panel_grid <- function(panels, ...) {
  cell <- function(p, letter) {
    wrap_elements(full = p) + labs(tag = letter) + PANEL_TAG + theme(plot.margin = margin(0, 0, 0, 0))
  }
  wrap_plots(Map(cell, panels, names(panels)), ...) + FLUSH
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

# "1st", "22nd", "93rd", "11th"
ordinal <- function(x) {
  x <- round(x)
  paste0(x, fifelse(x %% 100 %in% 11:13, "th", c("th", "st", "nd", "rd", rep("th", 6))[x %% 10 + 1]))
}

# "Mouse protein", "Human mRNA-adjusted protein", ... as a factor of the species-layer pairs of the data, mouse first
# and the layers in their usual order
LAYERS <- c("protein", "mRNA-adjusted protein", "mRNA")
layer_rows <- function(species, layer) {
  row <- paste(SPECIES_NAME[species], layer)
  rank <- order(match(species, names(SPECIES_NAME)), match(layer, LAYERS))
  factor(row, unique(row[rank]))
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
