# Shared style of the figures: the theme, recipes and composition of figures/pub at a Nature double column, the colour
# dictionary of the paper and the helpers that several figures share
# every character is black Arial (6 pt ticks, 7 pt axis titles, 8 pt bold letters)
# run from the repository root; the figure scripts read source_data/ and write output/figures/<name>.pdf, .png, .proof.png
suppressPackageStartupMessages(library(data.table))
source("figures/pub/pub.R")
pub_journal("nature")
stopifnot("Arial is not installed" = PUB_FONT == "Arial")

# one colour per entity across all figures; marks only, text stays black
MOUSE <- "#087C9A" # mouse discovery cohorts, held-out mouse
EXTERNAL <- "#EAC64B" # external mouse proteomes
HUMAN <- "#C3535E" # human cohorts and human disease
M6BRAIN <- "#2B3A67" # the m6Brain model and clock
LOAD <- "#8DB6D6" # the net guided by the m6A load
EDGE <- "#4D4D4D" # outlines of points and bars, connectors
CI <- "#8C8C8C" # confidence intervals
GREY <- "#BDBDBD" # controls, comparators, not significant
LIGHT <- "#E6E6E6" # context, bands
SET <- c("Held-out mouse" = MOUSE, "External mouse" = EXTERNAL, "Held-out human" = HUMAN)
SPECIES <- c(Mouse = MOUSE, Human = HUMAN)
CELL <- c(Neuron = "#4B3F8F", Astrocyte = "#EDA14D", Oligodendrocyte = "#A6D5A0", Microglia = "#B0417F",
          Ubiquitous = "#8C8C8C", "No cell-type data" = "#DEDEDE")
# number of m6A sites: an ordered ramp from grey to the hue of the cohort set
SITES <- c("0", "1", "2–3", "≥4")
sites_ramp <- function(hue) setNames(c(GREY, tint(hue, 0.6), tint(hue, 0.25), shade(hue, 0.15)), SITES)
# signed effects: teal, lower in old age, through white to brick, higher; the light ramp keeps black numbers legible
PUB_DIV <- c("#087C9A", "#3B96AE", "#7EBAC9", "#C3DFE6", "#F7F7F7", "#F3CFCF", "#E3999E", "#D06E76", "#B5404C")
PUB_DIV_LIGHT <- c("#6FB2C4", "#9CCBD7", "#CFE6EB", "#F7F7F7", "#F6D8D9", "#EBB0B4", "#DD8B92")
# the five targets of the gene-level models and the models themselves
TARGETS <- c("Mouse protein", "Mouse mRNA-adjusted protein", "Human protein", "Human mRNA", "Human neuronal mRNA")
SPECIES_NAME <- c(mouse = "Mouse", human = "Human")

theme_tile <- function() {
  theme_pub() + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
}
theme_rows <- function() theme_pub() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
theme_box <- function() theme_pub(border = TRUE) + theme(panel.spacing = unit(1.5, "mm"))
# an empty canvas in millimetres, for schematics
canvas <- function(w, h) {
  ggplot() + coord_fixed(xlim = c(0, w), ylim = c(0, h), expand = FALSE, clip = "off") + theme_void_pub()
}

# panels composed without aligning their axes, so that the row names of one panel do not squeeze the others
fig_panels <- function(..., design, heights) {
  fig_compose(lapply(list(...), function(p) wrap_elements(full = p)), design = design, heights = heights)
}

# text at the tick size; parse = TRUE for labels from lab_stats() and fmt_p()
note <- function(x, y, label, hjust = 0.5, vjust = 0.5, ...) {
  annotate("text", x = x, y = y, label = label, hjust = hjust, vjust = vjust, size = pt_size(6), ...)
}
note_corner <- function(label, x = -Inf, y = Inf, hjust = -0.06, vjust = 1.25, ...) {
  note(x, y, label, hjust = hjust, vjust = vjust, ...)
}
# a key inside a panel: one row per entry from (x, y) downwards, a short line (or square) and the name in black
key_rows <- function(labels, colours, x, y, step, width, gap = width * 0.35, square = FALSE) {
  rows <- y - (seq_along(labels) - 1) * step
  mark <- if (square) {
    annotate("point", x = x + width / 2, y = rows, colour = colours, shape = 15, size = 1.8)
  } else {
    annotate("segment", x = x, xend = x + width, y = rows, yend = rows, colour = colours, linewidth = lw(1))
  }
  list(mark, note(x + width + gap, rows, labels, hjust = 0))
}
# named points in a stacked column at x = `at`, from `top` downwards in steps of `step`, ordered as the points, so that
# the leaders cannot cross each other or a name; hjust = 0 puts the names right of the column, 1 left of it; `slant`
# moves each name along x, a staircase that keeps the leaders of the names above clear of the names below
label_column <- function(d, x, y, label, at, top, step, hjust = 0, fontface = "plain", slant = 0) {
  d <- as.data.table(d)[order(-get(y))]
  d[, `:=`(.at = at + (seq_len(.N) - 1) * slant, .y = top - (seq_len(.N) - 1) * step)]
  d[, .name := .at + (1 - 2 * hjust) * abs(at) * 0.025]
  list(
    geom_segment(data = d, aes(x = .data[[x]], y = .data[[y]], xend = .at, yend = .y), colour = "grey50",
                 linewidth = lw(0.3), inherit.aes = FALSE),
    geom_text(data = d, aes(.name, .y, label = .data[[label]]), hjust = hjust, fontface = fontface, size = pt_size(6),
              inherit.aes = FALSE)
  )
}
MINUS <- function(x) gsub("-", "−", x, fixed = TRUE)
# rounded half away from zero, as in the text (a binary 2.55 would otherwise print as 2.5)
half_up <- function(x, digits) round(x + sign(x) * 1e-9, digits) + 0
signed <- function(x, digits = 2) MINUS(sprintf(paste0("%+.", digits, "f"), half_up(x, digits)))
number <- function(x, digits = 2) MINUS(sprintf(paste0("%.", digits, "f"), half_up(x, digits)))
MINUS_LABELS <- scales::label_number(style_negative = "minus", drop0trailing = TRUE)
thousands <- function(x) format(as.integer(x), big.mark = ",", trim = TRUE)
# "r = -0.68, P = 4.5 x 10^-6, n = 18" as plotmath for parse = TRUE; the letter of r is set in `stat`
lab_stats <- function(value = NULL, p = NULL, n = NULL, stat = "italic(r)", digits = 2) {
  parts <- c(
    if (!is.null(value)) sprintf("%s == '%s'", stat, number(value, digits)),
    if (!is.null(p)) fmt_p(p),
    if (!is.null(n)) sprintf("italic(n) == '%s'", thousands(n))
  )
  if (length(parts) > 1) sprintf("list(%s)", paste(parts, collapse = ", ")) else parts
}

# fill of a signed heat table, squished at the limit (the end labels say so), the colour bar laid flat
scale_fill_signed <- function(limit, name, light = FALSE) {
  ends <- function(x) {
    out <- MINUS_LABELS(x)
    out[x == -limit] <- paste0("≤", out[x == -limit])
    out[x == limit] <- paste0("≥", out[x == limit])
    out
  }
  scale_fill_gradientn(
    colours = if (light) PUB_DIV_LIGHT else PUB_DIV, limits = c(-limit, limit), oob = scales::squish, name = name,
    na.value = "#F0F0F0", breaks = c(-limit, 0, limit), labels = ends,
    guide = guide_colourbar(
      direction = "horizontal", title.position = "top", title.hjust = 0.5, ticks = FALSE,
      theme = theme(legend.key.width = unit(16, "mm"), legend.key.height = unit(1.8, "mm"))
    )
  )
}

# colour of a dense scatter: points by the local density, from a tint of the hue to a shade of it
by_density <- function(d, x, y) {
  density <- MASS::kde2d(d[[x]], d[[y]], n = 120)
  d[, .density := density$z[cbind(findInterval(get(x), density$x), findInterval(get(y), density$y))]]
  d[order(.density)]
}
scale_density <- function(hue) {
  scale_colour_gradientn(colours = c(tint(hue, 0.72), hue, shade(hue, 0.5)), guide = "none")
}

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

# the figure at its printed size: vector PDF with embedded Arial, 600-dpi PNG and a proof at printed size
save_fig <- function(fig, name, height) {
  cairo_pdf(tempfile(fileext = ".pdf"), family = PUB_FONT) # the layout measures Arial, not the default device's font
  grob <- superscript_marks(if (inherits(fig, "patchwork")) patchworkGrob(fig) else ggplotGrob(fig))
  invisible(dev.off())
  save_pub(grob, file.path("output/figures", paste0(name, ".pdf")), width = "double", height = height)
}
