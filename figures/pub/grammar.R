# Layers for any panel: reference lines, labels on the data, P brackets, design strips and transcript cartoons


# --- Reference inside the panel ----------------------------------------------------
# Null, identity, chance, control median: the comparison the reader needs, drawn
# where the data is, not described in the legend.
ref_line <- function(yintercept = NULL, xintercept = NULL, slope = NULL, intercept = 0,
                     colour = "grey55", pt = 0.5) {
  args <- list(linetype = "22", linewidth = lw(pt), colour = colour)
  if (!is.null(slope))      return(do.call(geom_abline, c(list(slope = slope, intercept = intercept), args)))
  if (!is.null(yintercept)) return(do.call(geom_hline, c(list(yintercept = yintercept), args)))
  do.call(geom_vline, c(list(xintercept = xintercept), args))
}

# A shaded zone (no-effect band, noise floor, expected range). Add it first so
# it sits behind the data.
ref_band <- function(ymin = -Inf, ymax = Inf, xmin = -Inf, xmax = Inf, fill = "grey92")
  annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill, colour = NA)

# Vertical highlight through stacked tracks: same call on every track.
band_x <- function(xmin, xmax, fill = "#F3E3B5", alpha = 0.5)
  annotate("rect", xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill, alpha = alpha, colour = NA)


# --- Labels on the data ---------------------------------------------------------------
# Series named at the end of their line: no legend to decode. Text is always
# black Arial (house rule); the line itself carries the colour.
# Leave room on the right: scale_x_continuous(expand = expansion(mult = c(.02, .25))).
label_ends <- function(data, x, y, group, label = group, size = 6) {
  d <- as.data.frame(data)
  d <- do.call(rbind, lapply(split(d, d[[group]]), function(s) s[which.max(s[[x]]), , drop = FALSE]))
  ggrepel::geom_text_repel(data = d, aes(.data[[x]], .data[[y]], label = .data[[label]]),
                           colour = "black", hjust = 0, direction = "y", nudge_x = diff(range(data[[x]], na.rm = TRUE)) * 0.03,
                           size = pt_size(size), segment.size = lw(0.3), segment.colour = "grey60", family = PUB_FONT,
                           min.segment.length = 0.5, box.padding = 0.05, show.legend = FALSE, seed = 1)
}

# Selected points named with thin leaders (black text; the point carries the
# colour); everything else stays unlabelled. Pass only the rows to label.
label_points <- function(data, x, y, label, size = 6, ...) {
  m <- aes(.data[[x]], .data[[y]], label = .data[[label]])
  ggrepel::geom_text_repel(data = data, mapping = m, size = pt_size(size), segment.size = lw(0.3),
                           colour = "black", segment.colour = "grey40",
                           min.segment.length = 0, box.padding = 0.25, max.overlaps = Inf, family = PUB_FONT,
                           show.legend = FALSE, seed = 1, ...)
}

# Category labels carrying n — "WT\nn = 12" — for scale_x_discrete(labels = ...).
lab_n <- function(x, sep = "\n") {
  t <- table(x)
  stats::setNames(sprintf("%s%sn = %s", names(t), sep, format(as.integer(t), big.mark = ",", trim = TRUE)), names(t))
}


# --- Statistics where the eye is -------------------------------------------------------
# Exact P over a short bracket between two discrete positions (1-based)
p_bracket <- function(x1, x2, y, p, tip = NULL, size = 6) {
  colour <- "black"
  tip <- tip %||% (abs(y) * 0.03)
  list(annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = lw(0.5), colour = colour),
       annotate("segment", x = c(x1, x2), xend = c(x1, x2), y = y, yend = y - tip, linewidth = lw(0.5), colour = colour),
       annotate("text", x = (x1 + x2) / 2, y = y + tip * 0.6, label = fmt_p(p), parse = TRUE,
                vjust = 0, size = pt_size(size), colour = colour),
       expand_limits(y = y + tip * 5))  # headroom so the exponent is not clipped
}


# --- The design and the molecule as axes ----------------------------------------------
# +/- matrix under bars instead of rotated condition names. `design` has one row
# per x level (in plot order) and one column per factor; TRUE/FALSE become +/−,
# anything else is printed. Stack: (p + theme(axis.text.x = element_blank())) /
# design_strip(d) + plot_layout(heights = c(1, 0.1 * ncol(d)))
design_strip <- function(design, size = 6) {
  d <- as.data.frame(design)
  lv <- seq_len(nrow(d))
  long <- do.call(rbind, lapply(names(d), function(f) {
    v <- d[[f]]
    data.frame(x = lv, factor = f, lab = if (is.logical(v)) ifelse(v, "+", "−") else as.character(v))
  }))
  long$factor <- factor(long$factor, levels = rev(names(d)))
  ggplot(long, aes(x, factor, label = lab)) +
    geom_text(size = pt_size(size)) +
    scale_x_continuous(limits = c(0.5, nrow(d) + 0.5), expand = c(0, 0)) +
    theme_void_pub() +
    theme(axis.text.y = element_text(size = size, hjust = 1, margin = margin(r = 1.5, unit = "mm")),
          plot.margin = margin(0, 1.5, 0, 1.5, "mm"))
}

# Transcript cartoon for a metagene x axis: thin UTRs, thick CDS (light fill so
# its label stays black), UTR labels underneath.
# Expects the metagene on [0, length(regions)] with one unit per region.
tx_cartoon <- function(regions = c("5′UTR", "CDS", "3′UTR"), thick = c(FALSE, TRUE, FALSE),
                       fill = "grey75", size = 6) {
  n <- length(regions)
  d <- data.frame(x0 = seq_len(n) - 1, x1 = seq_len(n), h = ifelse(thick, 0.45, 0.2), lab = regions, thick = thick)
  ggplot(d) +
    geom_rect(aes(xmin = x0, xmax = x1, ymin = -h, ymax = h), fill = fill, colour = NA) +
    geom_text(data = d[d$thick, ], aes((x0 + x1) / 2, 0, label = lab), colour = "black", size = pt_size(size)) +
    geom_text(data = d[!d$thick, ], aes((x0 + x1) / 2, -0.7, label = lab), size = pt_size(size)) +
    scale_x_continuous(limits = c(0, n), expand = c(0, 0)) + scale_y_continuous(limits = c(-1, 0.5)) +
    theme_void_pub() + theme(plot.margin = margin(0, 1.5, 0, 1.5, "mm"))
}

# Gene model as layers on a genomic x: exons (type "cds" thick, "utr" thin) over
# an intron line. `exons`: data.frame(start, end, type).
gene_layers <- function(exons, y = 0, colour = "black", name = NULL, size = 6) {
  rng <- range(c(exons$start, exons$end))
  ly <- list(annotate("segment", x = rng[1], xend = rng[2], y = y, yend = y, linewidth = lw(0.5), colour = colour),
             annotate("rect", xmin = exons$start, xmax = exons$end,
                      ymin = y - ifelse(exons$type == "cds", 0.35, 0.18),
                      ymax = y + ifelse(exons$type == "cds", 0.35, 0.18), fill = colour, colour = NA))
  if (!is.null(name)) ly <- c(ly, list(annotate("text", x = rng[1], y = y - 0.75, label = name, hjust = 0,
                                                fontface = "italic", size = pt_size(size))))
  ly
}

# Short corner arrows instead of numbered axes on an embedding.
embed_axes <- function(x, y, labels = c("UMAP1", "UMAP2"), frac = 0.14, size = 6) {
  xr <- range(x, na.rm = TRUE); yr <- range(y, na.rm = TRUE)
  x0 <- xr[1]; y0 <- yr[1]; dx <- diff(xr) * frac; dy <- diff(yr) * frac
  a <- grid::arrow(length = unit(0.9, "mm"), type = "closed")
  list(annotate("segment", x = x0, xend = x0 + dx, y = y0, yend = y0, arrow = a, linewidth = lw(0.5)),
       annotate("segment", x = x0, xend = x0, y = y0, yend = y0 + dy, arrow = a, linewidth = lw(0.5)),
       annotate("text", x = x0 + dx / 2, y = y0 - diff(yr) * 0.035, label = labels[1], size = pt_size(size)),
       annotate("text", x = x0 - diff(xr) * 0.035, y = y0 + dy / 2, label = labels[2], angle = 90, size = pt_size(size)))
}

# One-line key above a panel instead of a boxed legend: a colour square, then the
# name in black. Stack: key_inline(cols) / p + plot_layout(heights = c(.06, 1)).
key_inline <- function(colours, size = 6, gap = 3) {
  w <- nchar(names(colours)) + gap
  d <- data.frame(lab = names(colours), x = cumsum(c(0, utils::head(w, -1))), col = unname(colours))
  ggplot(d, aes(x, 0)) +
    geom_point(aes(colour = I(col)), shape = 15, size = 1.8) +
    geom_text(aes(x = x + 1.2, label = lab), hjust = 0, size = pt_size(size), colour = "black") +
    scale_x_continuous(limits = c(0, sum(w)), expand = c(0, 0)) +
    theme_void_pub() + theme(plot.margin = margin(0, 1.5, 0, 1.5, "mm"))
}

invisible(NULL)
