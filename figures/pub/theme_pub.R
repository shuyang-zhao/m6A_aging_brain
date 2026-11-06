# Theme, journal sizes, save and small helpers of the figures; every figure takes its style from here

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

# --- Fonts -------------------------------------------------------------------
# fontconfig answers "Helvetica" with a lookalike and never errors, so resolve
# against faces actually installed and say what was picked.
pub_font <- function(prefer = c("Arial", "Helvetica", "Helvetica Neue",
                                "Liberation Sans", "Nimbus Sans", "DejaVu Sans"),
                     quiet = FALSE) {
  if (!requireNamespace("systemfonts", quietly = TRUE)) return(prefer[[1]])
  have <- tolower(unique(systemfonts::system_fonts()$family))
  hit  <- prefer[tolower(prefer) %in% have]
  if (!length(hit)) return("sans")
  if (!identical(hit[[1]], "Arial"))
    warning("Arial is not installed; using '", hit[[1]], "'. Install it (msttcorefonts) before drawing the figures.",
            call. = FALSE)
  hit[[1]]
}
PUB_FONT <- pub_font(quiet = TRUE)
# Every character on a figure is black Arial; colour belongs to marks
PUB_INK <- "#000000"


# --- Journals ----------------------------------------------------------------
# From the author guides of the journals; mm, type in pt.
#   text: c(min, axis title, tick/legend); tag: panel-label size and case.
PUB_JOURNALS <- list(
  nature  = list(cols = c(single = 89, onehalf = 120, double = 183), max_height = 170,
                 text = c(min = 5, title = 7, tick = 6), tag = 8, case = "a", line_min = 0.25),
  # NG, NBT, NCB, NSMB, NMI, NMeth share one guide: 5–7 pt, max width 180 mm.
  natfam  = list(cols = c(single = 89, double = 180), max_height = 170,
                 text = c(min = 5, title = 7, tick = 6), tag = 8, case = "a", line_min = 0.25),
  # Cell Press (Cell, Mol Cell, Cell Systems, Cell Genomics): 6–8 pt Arial, 0.5–1.5 pt strokes.
  cell    = list(cols = c(single = 85, onehalf = 114, double = 174), max_height = 200,
                 text = c(min = 6, title = 7, tick = 6), tag = 9, case = "A", line_min = 0.5),
  science = list(cols = c(single = 57, double = 121, triple = 184), max_height = 230,
                 text = c(min = 5, title = 7, tick = 6), tag = 10, case = "A", line_min = 0.5),
  sciadv  = list(cols = c(single = 90, onehalf = 144, double = 184), max_height = 199,
                 text = c(min = 6, title = 8, tick = 7), tag = 9, case = "A", line_min = 0.28),
  nar     = list(cols = c(single = 84, double = 178), max_height = 230,
                 text = c(min = 6, title = 7, tick = 6), tag = 9, case = "A", line_min = 0.5)
)
PUB_ALIAS <- c(ng = "natfam", nbt = "natfam", ncb = "natfam", nsmb = "natfam", nmi = "natfam",
               nmeth = "natfam", natgen = "natfam", molcell = "cell", cellsys = "cell",
               cellgen = "cell", cellsystems = "cell", cellgenomics = "cell")

pub_spec <- function(journal = pub_journal()) {
  j <- tolower(journal)
  if (j %in% names(PUB_ALIAS)) j <- PUB_ALIAS[[j]]
  s <- PUB_JOURNALS[[j]]
  if (is.null(s)) stop("unknown journal '", journal, "'; known: ",
                       paste(c(names(PUB_JOURNALS), names(PUB_ALIAS)), collapse = ", "), call. = FALSE)
  s
}

# Set once at the top of a project's style file; theme, tags and save follow it.
pub_journal <- function(journal = NULL) {
  if (is.null(journal)) return(getOption("pub.journal", "nature"))
  pub_spec(journal)
  options(pub.journal = journal)
  invisible(journal)
}

pub_width <- function(journal = pub_journal(), col = "single") {
  w <- pub_spec(journal)$cols
  if (!col %in% names(w)) stop("'", journal, "' has no '", col, "' width; have ",
                               paste(names(w), collapse = ", "), call. = FALSE)
  unname(w[[col]])
}
PUB_WIDTHS <- lapply(PUB_JOURNALS, function(s) c(s$cols, max_height = s$max_height))  # v1 name


# --- Units ---------------------------------------------------------------------
# geom_text(size=) is in mm, so size = 7 is ~20 pt.
pt_size <- function(size_pt) size_pt / .pt
# ggplot linewidth -> PDF points: lwd = linewidth * .pt, and R's lwd 1 = 0.75 pt.
lw <- function(pt) pt / (.pt * 0.75)
# Point outlines use a different unit (stroke * .stroke / 2): st(0.4) = a 0.4 pt outline.
st <- function(pt) pt / 0.75 * 2 / .stroke

# theme_void() that keeps the house font for axis text *and* geom text; plain
# theme_void() lets annotate()/geom_text() fall back to "sans" (Liberation/DejaVu).
theme_void_pub <- function(size = 6) {
  th <- theme_void(base_size = size, base_family = PUB_FONT)
  if (utils::packageVersion("ggplot2") >= "4.0.0") th <- th + theme(geom = element_geom(family = PUB_FONT, fontsize = size))
  th
}


# --- Theme -----------------------------------------------------------------------
# Open L axes, outward ticks, no grid, one sans face, axis titles one step above tick labels, black text, small
# borderless legends; `grid = "y"` adds faint reference rules for reading values off horizontal bars.
theme_pub <- function(base_size = NULL, base_family = PUB_FONT, border = FALSE,
                      journal = pub_journal(), grid = c("none", "x", "y", "xy")) {
  s <- pub_spec(journal)
  grid <- match.arg(grid)
  title <- if (is.null(base_size)) s$text[["title"]] else base_size
  tick  <- if (is.null(base_size)) s$text[["tick"]] else max(base_size - 1, s$text[["min"]])
  rule  <- element_line(linewidth = lw(0.25), colour = "grey88")
  ln    <- element_line(linewidth = lw(max(0.5, s$line_min)), colour = "black")
  th <- theme_bw(base_size = title, base_family = base_family) +
    theme(
      text              = element_text(colour = "black", size = title),
      axis.text         = element_text(size = tick, colour = "black"),
      axis.title        = element_text(size = title),
      axis.ticks        = ln,
      axis.ticks.length = unit(1, "mm"),
      axis.line         = if (border) element_blank() else ln,
      panel.border      = if (border) element_rect(linewidth = lw(max(0.5, s$line_min)),
                                                   colour = "black", fill = NA) else element_blank(),
      panel.grid        = element_blank(),
      panel.grid.major.x = if (grid %in% c("x", "xy")) rule else element_blank(),
      panel.grid.major.y = if (grid %in% c("y", "xy")) rule else element_blank(),
      panel.background  = element_blank(),
      plot.background   = element_blank(),
      legend.key        = element_blank(),
      legend.background = element_blank(),
      legend.title      = element_text(size = tick),
      legend.text       = element_text(size = tick),
      legend.key.size   = unit(2.6, "mm"),
      legend.margin     = margin(0, 0, 0, 0),
      legend.box.spacing = unit(1.5, "mm"),
      strip.background  = element_blank(),
      strip.text        = element_text(size = title, colour = "black", margin = margin(1, 1, 1, 1, "mm")),
      plot.title        = element_text(size = title, face = "plain", hjust = 0, colour = "black"),
      plot.tag          = element_text(size = s$tag, face = "bold", colour = "black"),
      plot.margin       = margin(1.5, 1.5, 1.5, 1.5, unit = "mm")
    )
  # ggplot2 >= 4 carries geom defaults in the theme, so every geom gets 0.75 pt
  # lines and 6 pt text without repeating it at each call.
  if (utils::packageVersion("ggplot2") >= "4.0.0")
    th <- th + theme(geom = element_geom(ink = "black", linewidth = lw(0.75), borderwidth = lw(0.5),
                                         pointsize = 0.8, family = base_family, fontsize = tick))
  th
}


# --- Colour constants ------------------------------------------------------------
okabe_ito <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
okabe_ito_named <- c(black = "#000000", orange = "#E69F00", skyblue = "#56B4E9", green = "#009E73",
                     yellow = "#F0E442", blue = "#0072B2", vermillion = "#D55E00", purple = "#CC79A7")
# Greys of not significant results, controls, background and context
pub_grey   <- "#BDBDBD"
pub_greys  <- c(dark = "#4D4D4D", mid = "#9E9E9E", ns = "#BDBDBD", light = "#D9D9D9", faint = "#EDEDED")
pub_accent <- c(down = "#2166AC", up = "#B2182B", ns = pub_grey)
# Muted diverging (RdBu) with a near-white centre; contact-map ramp ("fall", cooltools).
PUB_DIV  <- c("#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#EF8A62", "#B2182B")
# Same map stopping at mid-saturation: black text stays readable on every cell.
PUB_DIV_LIGHT <- c("#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D")
PUB_FALL <- c("#FFFFFF", "#F2F1BF", "#F6E08A", "#F6B95A", "#F18C43", "#E7442C", "#D42126", "#961224", "#120609")

scale_colour_pub <- function(...) scale_colour_manual(values = unname(okabe_ito), ...)
scale_fill_pub   <- function(...) scale_fill_manual(values = unname(okabe_ito), ...)
scale_color_pub  <- scale_colour_pub

# Symmetric limits so white sits on zero and two panels can share one key.
lim_sym <- function(x, q = 1) { m <- stats::quantile(abs(x), q, na.rm = TRUE); c(-m, m) }
scale_fill_div <- function(limits = NULL, name = waiver(), ...)
  scale_fill_gradientn(colours = PUB_DIV, limits = limits, oob = scales::squish, name = name, ...)
scale_colour_div <- function(limits = NULL, name = waiver(), ...)
  scale_colour_gradientn(colours = PUB_DIV, limits = limits, oob = scales::squish, name = name, ...)


# --- Save --------------------------------------------------------------------------
# .pdf vector with embedded editable fonts (submit), .png 600 dpi (slides),
# .proof.png at ~printed size (LOOK at this one), optional .source.csv.
save_pub <- function(p, file, width = "single", height = 60, journal = pub_journal(),
                     dpi = 600, proof = TRUE, proof_dpi = 150, svg = FALSE,
                     data = NULL, verbose = TRUE) {
  journal <- journal %||% pub_journal()
  if (is.character(width)) width <- pub_width(journal, width)
  mx <- pub_spec(journal)$max_height
  if (height > mx) warning(sprintf("height %.0f mm exceeds %s's %.0f mm; production will shrink every label.",
                                   height, journal, mx), call. = FALSE)
  if (!grepl("\\.pdf$", file, ignore.case = TRUE)) stop("`file` must end in .pdf", call. = FALSE)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  # cairo_pdf embeds and subsets the real font; showtext would outline glyphs,
  # which Nature rejects.
  # symbolfamily: plotmath draws minus/times/Greek from a Symbol font, which
  # fontconfig resolves to DejaVu; point it at the house face.
  pdf_dev <- function(filename, ...) grDevices::cairo_pdf(filename, ...,
    symbolfamily = grDevices::cairoSymbolFont(PUB_FONT, usePUA = FALSE))
  ggsave(file, p, width = width, height = height, units = "mm", device = pdf_dev)
  # Pass "png" + type, not a wrapper function: ggsave hands units/res only to
  # devices it knows, and a wrapper gets a ~2-pixel canvas.
  png_save <- function(f, res) {
    if (requireNamespace("ragg", quietly = TRUE))
      ggsave(f, p, width = width, height = height, units = "mm", dpi = res, bg = "white", device = ragg::agg_png)
    else ggsave(f, p, width = width, height = height, units = "mm", dpi = res, bg = "white", device = "png", type = "cairo")
  }
  out <- c(pdf = file)
  out["png"] <- sub("\\.pdf$", ".png", file, ignore.case = TRUE)
  png_save(out[["png"]], dpi)
  if (proof) {
    out["proof"] <- sub("\\.pdf$", ".proof.png", file, ignore.case = TRUE)
    png_save(out[["proof"]], proof_dpi)
  }
  if (svg && requireNamespace("svglite", quietly = TRUE))
    ggsave(sub("\\.pdf$", ".svg", file), p, width = width, height = height, units = "mm", device = svglite::svglite)
  if (!is.null(data)) {
    out["data"] <- sub("\\.pdf$", ".source.csv", file, ignore.case = TRUE)
    utils::write.csv(data, out[["data"]], row.names = FALSE)
  }
  if (verbose) message(sprintf("save_pub: %.0f x %.0f mm -> %s  (look at %s)", width, height,
                               basename(file), basename(out[["proof"]] %||% out[["png"]])))
  invisible(out)
}
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a


# --- Rasterise dense layers ----------------------------------------------------------
# Tens of thousands of vector points bloat the PDF and stall Illustrator. The
# layer is drawn into a cairo bitmap at the panel's printed size at draw time;
# axes and text stay vector. No ggrastr/Cairo needed (neither builds here).
rasterise_pts <- function(layer, dpi = 600) {
  parent <- layer$geom
  layer$geom <- ggproto(NULL, parent,
    draw_panel = function(self, data, panel_params, coord, ...) {
      g <- ggproto_parent(parent, self)$draw_panel(data, panel_params, coord, ...)
      grob(child = g, dpi = dpi, cl = "pub_raster")
    })
  layer
}
rasterise_layer <- rasterise_pts

drawDetails.pub_raster <- function(x, recording) {
  w <- convertWidth(unit(1, "npc"), "in", valueOnly = TRUE)
  h <- convertHeight(unit(1, "npc"), "in", valueOnly = TRUE)
  f <- tempfile(fileext = ".png"); dev <- grDevices::dev.cur()
  grDevices::png(f, width = max(1, round(w * x$dpi)), height = max(1, round(h * x$dpi)),
                 res = x$dpi, type = "cairo", bg = "transparent")
  grid.newpage(); pushViewport(viewport())  # a fresh device has no root viewport yet
  grid.draw(x$child)
  grDevices::dev.off(); grDevices::dev.set(dev)
  grid.raster(png::readPNG(f), interpolate = FALSE)
  unlink(f)
}


# --- Small helpers ---------------------------------------------------------------------
# Bars sit on the axis; the default expansion floats them over a white sliver.
expand_bar_y <- function(...) scale_y_continuous(expand = expansion(mult = c(0, 0.05)), ...)

# Exact P for plotmath (parse = TRUE): "italic(P) == 3.4 %*% 10^-5".
fmt_p <- function(p, digits = 2, prefix = "italic(P)") {
  vapply(p, function(v) {
    if (is.na(v)) return("")
    if (v < 1e-300) return(paste(prefix, "< 10^-300"))
    if (v >= 0.001) return(sprintf("%s == '%s'", prefix, formatC(signif(v, digits), format = "fg", digits = digits)))
    e <- floor(log10(v)); m <- signif(v / 10^e, digits)
    if (m >= 10) { m <- m / 10; e <- e + 1 }
    sprintf("%s == %s %%*%% 10^%d", prefix, format(m), e)
  }, "")
}

invisible(NULL)
