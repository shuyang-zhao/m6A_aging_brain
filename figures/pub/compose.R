# Assembly of panels into one figure at its printed size


# patchwork with journal panel labels (a, b, c for Nature family; A, B, C for
# Cell Press / Science) at the journal's tag size. A recipe that already returns
# a stacked patchwork (tracks, metagene + cartoon, tornado, butterfly) is wrapped
# so it gets one label, not one per sub-plot.
#   widths/heights: relative numbers, or unit(c(40, 80), "mm") for fixed sizes
#   design: patchwork layout string, e.g. "AAB\nCDE"
fig_compose <- function(..., design = NULL, widths = NULL, heights = NULL, ncol = NULL, nrow = NULL,
                        guides = "keep", journal = pub_journal(), tags = TRUE) {
  s <- pub_spec(journal)
  ps <- list(...)
  if (length(ps) == 1 && is.list(ps[[1]]) && !inherits(ps[[1]], "gg")) ps <- ps[[1]]
  ps <- lapply(ps, function(p) if (inherits(p, "patchwork")) patchwork::wrap_elements(full = p) else p)
  out <- patchwork::wrap_plots(ps, design = design, widths = widths, heights = heights,
                               ncol = ncol, nrow = nrow, guides = guides)
  if (tags) out <- out + patchwork::plot_annotation(tag_levels = s$case) &
    theme(plot.tag = element_text(size = s$tag, face = "bold", family = PUB_FONT),
          plot.tag.position = c(0, 1))
  out
}

# Reserved space for a schematic, micrograph, structure render or blot that will
# be drawn elsewhere (BioRender, ChimeraX, Illustrator) and dropped in. Keeps the
# grid honest: the figure is composed at its true size from the start.
slot <- function(label = "schematic", fill = "grey97") {
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = fill, colour = "grey70",
             linetype = "22", linewidth = lw(0.5)) +
    annotate("text", x = 0.5, y = 0.5, label = label, colour = "black", size = pt_size(6)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) + theme_void_pub()
}

# A finished raster (PNG) placed as a panel — a rendered schematic or image.
# Export it at >= 300 dpi at the size it will occupy.
img_panel <- function(file) {
  patchwork::wrap_elements(full = grid::rasterGrob(png::readPNG(file), interpolate = TRUE))
}

invisible(NULL)
