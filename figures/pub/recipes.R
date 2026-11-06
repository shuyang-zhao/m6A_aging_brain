# Panel types: each takes a data frame and column names and returns a ggplot (stacked ones a patchwork);
# `colours` is a named vector from the colour dictionary of style.R


.cols <- function(lv, colours = NULL) {
  lv <- as.character(lv)
  if (!is.null(colours)) {
    miss <- setdiff(lv, names(colours))
    if (length(miss)) stop("no colour for: ", paste(miss, collapse = ", "), call. = FALSE)
    return(colours[lv])
  }
  # control grey and one accent for two groups; otherwise Okabe-Ito
  if (length(lv) == 2) return(stats::setNames(c("#8C8C8C", "#D55E00"), lv))
  stats::setNames(rep_len(c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"), length(lv)), lv)
}
.fac <- function(v) if (is.factor(v)) droplevels(v) else factor(v, levels = unique(v))
.msd <- function(v) data.frame(y = mean(v), ymin = mean(v) - stats::sd(v), ymax = mean(v) + stats::sd(v))
.mse <- function(v) { s <- stats::sd(v) / sqrt(length(v)); data.frame(y = mean(v), ymin = mean(v) - s, ymax = mean(v) + s) }
.maybe_raster <- function(layer, n, at = 5000) if (n > at) rasterise_pts(layer) else layer


# ==== Comparisons between groups (wet lab and general) =========================

# Mean bar (tint) + every replicate (full colour) + SD; y from 0. Cell Press
# requires the individual points. `design`: one row per x level -> +/- matrix
# under the axis instead of rotated condition names.
fig_bar_dots <- function(df, x, y, fill = x, colours = NULL, err = c("sd", "se"), design = NULL,
                         ylab = y, xlab = NULL, dodge = 0.75) {
  err <- match.arg(err); d <- as.data.frame(df)
  d[[x]] <- .fac(d[[x]]); d[[fill]] <- .fac(d[[fill]])
  cols <- .cols(levels(d[[fill]]), colours)
  pd <- position_dodge(width = dodge)
  g <- ggplot(d, aes(.data[[x]], .data[[y]], fill = .data[[fill]], colour = .data[[fill]])) +
    stat_summary(fun = mean, geom = "col", width = 0.62, linewidth = lw(0.5), position = pd) +
    stat_summary(fun.data = if (err == "sd") .msd else .mse, geom = "errorbar", width = 0.22,
                 linewidth = lw(0.5), colour = "black", position = pd) +
    geom_point(position = position_jitterdodge(jitter.width = 0.18, dodge.width = dodge, seed = 1),
               size = 1.1, stroke = 0) +
    scale_fill_manual(values = tint(cols, 0.6), guide = if (fill == x) "none" else "legend") +
    scale_colour_manual(values = cols, guide = "none") +
    expand_bar_y() + labs(x = xlab, y = ylab, fill = NULL) + theme_pub()
  if (is.null(design)) return(g)
  (g + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())) /
    design_strip(design) + patchwork::plot_layout(heights = c(1, 0.09 * ncol(design)))
}

# Capless box (tint fill, hue outline) with every point; n under each label.
# violin = TRUE for large n; points drop out above `max_points` per group.
fig_box <- function(df, x, y, colour = x, colours = NULL, violin = FALSE, points = TRUE,
                    n_labels = TRUE, ref = NULL, ylab = y, xlab = NULL, max_points = 1500) {
  d <- as.data.frame(df); d[[x]] <- .fac(d[[x]]); d[[colour]] <- .fac(d[[colour]])
  cols <- .cols(levels(d[[colour]]), colours)
  g <- ggplot(d, aes(.data[[x]], .data[[y]], colour = .data[[colour]], fill = .data[[colour]]))
  if (!is.null(ref)) g <- g + ref_line(yintercept = ref)
  if (violin) g <- g + geom_violin(colour = NA, scale = "width", width = 0.8) +
    geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", linewidth = lw(0.5), staplewidth = 0)
  else g <- g + geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = lw(0.6), staplewidth = 0)
  if (points && max(table(d[[x]])) <= max_points)
    g <- g + .maybe_raster(ggbeeswarm::geom_quasirandom(width = if (violin) 0.3 else 0.18, size = 0.6,
                                                        stroke = 0, alpha = 0.7), nrow(d))
  g + scale_colour_manual(values = cols, guide = "none") +
    scale_fill_manual(values = tint(cols, if (violin) 0.7 else 0.75), guide = "none") +
    (if (n_labels) scale_x_discrete(labels = lab_n(d[[x]])) else NULL) +
    labs(x = xlab, y = ylab) + theme_pub()
}

# Matched samples joined by grey lines; group mean as a heavy segment.
fig_paired <- function(df, id, x, y, colours = NULL, ylab = y, xlab = NULL) {
  d <- as.data.frame(df); d[[x]] <- .fac(d[[x]]); cols <- .cols(levels(d[[x]]), colours)
  ggplot(d, aes(.data[[x]], .data[[y]])) +
    geom_line(aes(group = .data[[id]]), colour = "grey75", linewidth = lw(0.5)) +
    geom_point(aes(colour = .data[[x]]), size = 1.3) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.35, linewidth = lw(0.8), colour = "black",
                 fun.min = mean, fun.max = mean) +
    scale_colour_manual(values = cols, guide = "none") + labs(x = xlab, y = ylab) + theme_pub()
}

# ECDF family; ordinal = TRUE maps an ordered stratifier to a one-hue ramp. The legend sits inside and carries n, with
# a dashed cross at the median line and the reference value.
fig_ecdf <- function(df, value, group, colours = NULL, ordinal = FALSE, test = TRUE, ref = NULL,
                     xlab = value, ylab = "Cumulative fraction", hue = "#08519C") {
  d <- as.data.frame(df); d <- d[!is.na(d[[value]]) & !is.na(d[[group]]), ]
  d[[group]] <- .fac(d[[group]]); lv <- levels(d[[group]])
  cols <- if (ordinal) stats::setNames(grDevices::colorRampPalette(c(tint(hue, 0.75), hue, shade(hue, 0.45)))(length(lv)), lv)
          else .cols(lv, colours)
  n <- table(d[[group]])
  g <- ggplot(d, aes(.data[[value]], colour = .data[[group]])) +
    ref_line(yintercept = 0.5) + (if (!is.null(ref)) ref_line(xintercept = ref) else NULL) +
    stat_ecdf(geom = "step", linewidth = lw(0.9)) +
    scale_colour_manual(values = cols, labels = sprintf("%s (n = %s)", lv, format(as.integer(n[lv]), big.mark = ",")),
                        name = NULL) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.01)) +
    labs(x = xlab, y = ylab) + theme_pub() +
    theme(legend.position = "inside", legend.position.inside = c(0.98, 0.03),
          legend.justification = c(1, 0), legend.key.height = unit(2.4, "mm"))
  if (test && length(lv) == 2) {
    sp <- split(d[[value]], d[[group]])
    p <- suppressWarnings(stats::ks.test(sp[[1]], sp[[2]])$p.value)
    g <- g + annotate("text", x = -Inf, y = Inf, label = fmt_p(max(p, 1e-300)), parse = TRUE,
                      hjust = -0.08, vjust = 1.4, size = pt_size(6))
  }
  g
}


# ==== Effects across many categories — the alternatives to a forest plot ======

# Dot matrix: colour = signed effect (fixed symmetric scale), size = significance;
# non-significant cells as small open grey circles so the pattern reads at once.
fig_dotmatrix <- function(df, x, y, effect, sig, sig_cut = 0.05, limits = NULL, row_group = NULL,
                          fill_lab = "Effect", size_lab = expression(-log[10]~FDR), cap = 10, xlab = NULL, ylab = NULL) {
  d <- as.data.frame(df); d$.s <- pmin(-log10(d[[sig]]), cap); d$.sig <- d[[sig]] < sig_cut
  d[[x]] <- .fac(d[[x]]); d[[y]] <- .fac(d[[y]])
  g <- ggplot(d, aes(.data[[x]], .data[[y]])) +
    geom_point(data = d[!d$.sig, ], shape = 1, size = 0.9, colour = "grey70", stroke = st(0.5)) +
    geom_point(data = d[d$.sig, ], aes(size = .data$.s, fill = .data[[effect]]), shape = 21,
               colour = "grey20", stroke = st(0.4)) +
    scale_fill_div(limits = limits %||% lim_sym(d[[effect]]), name = fill_lab) +
    scale_size_area(max_size = 3.4, name = size_lab) +
    labs(x = xlab, y = ylab) + theme_pub() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          panel.grid.major = element_line(linewidth = lw(0.25), colour = "grey93"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.key.width = unit(2.2, "mm"))
  if (!is.null(row_group))
    g <- g + facet_grid(rows = vars(.data[[row_group]]), scales = "free_y", space = "free_y", switch = "y") +
      theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
            panel.spacing.y = unit(1.2, "mm"))
  g
}

# Two-sided counts or effects per category (up vs down, condition A vs B):
# mirrored bars with the category names in the gutter between them.
fig_butterfly <- function(df, label, left, right, colours = c(left = "#2166AC", right = "#B2182B"),
                          titles = c(left, right), order = right, xlab = NULL) {
  d <- as.data.frame(df); d[[label]] <- factor(d[[label]], levels = d[[label]][order(d[[order]])])
  side <- function(col, fill, mirror) {
    ggplot(d, aes(.data[[col]], .data[[label]])) +
      geom_col(fill = fill, width = 0.7) +
      scale_x_continuous(trans = if (mirror) "reverse" else "identity", expand = expansion(mult = c(0, 0.05))) +
      scale_y_discrete(position = if (mirror) "right" else "left") +
      labs(x = xlab, y = NULL, title = if (mirror) titles[1] else titles[2]) + theme_pub() +
      theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
            plot.title = element_text(hjust = if (mirror) 1 else 0, face = "bold"),
            axis.text.y.right = element_text(hjust = 0.5, margin = margin(l = 1, r = 1, unit = "mm")),
            axis.text.y.left = if (mirror) element_text() else element_blank())
  }
  side(left, colours[["left"]], TRUE) + side(right, colours[["right"]], FALSE) +
    patchwork::plot_layout(widths = c(1, 1))
}

# Every item ranked (S-curve), null/ns in grey, the tails coloured and named.
# For screen scores, per-gene effects, enrichment over thousands of terms.
fig_scurve <- function(df, score, label = NULL, sig = NULL, top_n = 6, colours = pub_accent,
                       ylab = score, xlab = "Rank") {
  d <- as.data.frame(df); d <- d[order(d[[score]]), ]; d$.rank <- seq_len(nrow(d))
  hit <- if (is.null(sig)) d$.rank <= top_n | d$.rank > nrow(d) - top_n else d[[sig]]
  d$.cls <- ifelse(!hit, "ns", ifelse(d[[score]] > 0, "up", "down"))
  g <- ggplot(d, aes(.data$.rank, .data[[score]], colour = .data$.cls)) +
    ref_line(yintercept = 0) +
    .maybe_raster(geom_point(size = 0.5, stroke = 0), nrow(d)) +
    scale_colour_manual(values = colours, guide = "none") +
    scale_x_continuous(expand = expansion(mult = 0.02)) + labs(x = xlab, y = ylab) + theme_pub()
  if (!is.null(label)) {
    tails <- rbind(utils::head(d[d$.cls == "down", ], top_n), utils::tail(d[d$.cls == "up", ], top_n))
    g <- g + label_points(tails, ".rank", score, label, fontface = "italic",
                          nudge_x = ifelse(tails$.cls == "up", -nrow(d) * 0.08, nrow(d) * 0.08))
  }
  g
}

# Table + estimate in one object: text columns (n, OR, P) aligned row-by-row with
# a cap-less CI on a log axis, sorted, coloured by direction, banded rows.
fig_tableforest <- function(df, label, est, lo, hi, text = NULL, null = 1, log = TRUE,
                            colours = pub_accent, xlab = "Odds ratio", widths = NULL) {
  d <- as.data.frame(df); d <- d[order(d[[est]]), ]; d$.y <- seq_len(nrow(d))
  d$.dir <- ifelse(d[[lo]] > null, "up", ifelse(d[[hi]] < null, "down", "ns"))
  band <- geom_rect(data = d[d$.y %% 2 == 0, ], aes(ymin = .data$.y - 0.5, ymax = .data$.y + 0.5),
                    xmin = -Inf, xmax = Inf, fill = "grey95", inherit.aes = FALSE)
  yscale <- scale_y_continuous(breaks = d$.y, labels = d[[label]], expand = c(0, 0.6))
  ci <- ggplot(d, aes(.data[[est]], .data$.y, colour = .data$.dir)) + band +
    ref_line(xintercept = null) +
    geom_segment(aes(x = .data[[lo]], xend = .data[[hi]], yend = .data$.y), linewidth = lw(0.9)) +
    geom_point(size = 1.4) + scale_colour_manual(values = colours, guide = "none") +
    (if (log) scale_x_log10() else NULL) + yscale + labs(x = xlab, y = NULL) + theme_pub() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
  if (is.null(text)) return(ci)
  # column widths follow the longest entry, so numbers never run into each other
  hdr <- names(text) %||% text
  w <- vapply(seq_along(text), function(i) max(nchar(c(hdr[i], as.character(d[[text[[i]]]])))), 1)
  pos <- cumsum(w + 1.5)
  long <- do.call(rbind, lapply(seq_along(text), function(i)
    data.frame(.y = d$.y, col = pos[i], txt = as.character(d[[text[[i]]]]))))
  tb <- ggplot(long, aes(.data$col, .data$.y, label = .data$txt)) + band +
    geom_text(size = pt_size(6), hjust = 1) +
    scale_x_continuous(breaks = pos, labels = hdr, position = "top", limits = c(0, max(pos) + 0.3), expand = c(0, 0)) +
    scale_y_continuous(breaks = NULL, expand = c(0, 0.6)) + labs(x = NULL, y = NULL) + theme_pub() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x.top = element_text(hjust = 1))
  ci + tb + patchwork::plot_layout(widths = widths %||% c(1.4, 0.075 * max(pos)))
}

# Heat table with the value printed in each tile: exact numbers for a few dozen
# cells, colour for the pattern. Fixed symmetric limits, grey for missing. The
# ramp stops at mid-saturation so black numbers stay legible on every tile.
fig_heat_table <- function(df, x, y, value, label = NULL, digits = 2, limits = NULL, diverging = TRUE,
                           hue = "#08519C", fill_lab = value, xlab = NULL, ylab = NULL) {
  d <- as.data.frame(df); d[[x]] <- .fac(d[[x]]); d[[y]] <- .fac(d[[y]])
  lim <- limits %||% (if (diverging) lim_sym(d[[value]]) else range(d[[value]], na.rm = TRUE))
  ramp <- if (diverging) PUB_DIV_LIGHT else c("#FFFFFF", tint(hue, 0.7), tint(hue, 0.35))
  d$.lab <- if (is.null(label)) ifelse(is.na(d[[value]]), "", formatC(d[[value]], format = "f", digits = digits)) else d[[label]]
  ggplot(d, aes(.data[[x]], .data[[y]], fill = .data[[value]])) +
    geom_tile(colour = "white", linewidth = lw(0.6)) +
    geom_text(aes(label = .data$.lab), colour = "black", size = pt_size(6)) +
    scale_fill_gradientn(colours = ramp, limits = lim, oob = scales::squish, na.value = "grey90", name = fill_lab) +
    scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(limits = rev, expand = c(0, 0)) +
    labs(x = xlab, y = ylab) + theme_pub() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.key.width = unit(2.2, "mm"))
}

# Bars that start at the null (OR = 1, log2FC = 0), exact P typed at the bar end.
fig_null_bars <- function(df, label, est, p = NULL, null = 1, log = TRUE, colours = pub_accent,
                          xlab = "Odds ratio") {
  d <- as.data.frame(df); d <- d[order(d[[est]]), ]; d$.y <- seq_len(nrow(d))
  d$.dir <- ifelse(d[[est]] > null, "up", "down")
  g <- ggplot(d) +
    geom_rect(aes(xmin = pmin(null, .data[[est]]), xmax = pmax(null, .data[[est]]),
                  ymin = .data$.y - 0.36, ymax = .data$.y + 0.36, fill = .data$.dir)) +
    ref_line(xintercept = null, colour = "black") +
    scale_fill_manual(values = colours, guide = "none") +
    scale_y_continuous(breaks = d$.y, labels = d[[label]], expand = c(0, 0.6)) +
    (if (log) scale_x_log10(expand = expansion(mult = c(0.05, 0.35))) else
       scale_x_continuous(expand = expansion(mult = c(0.05, 0.35)))) +
    coord_cartesian(clip = "off") + labs(x = xlab, y = NULL) + theme_pub() +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(1.5, 12, 1.5, 1.5, "mm"))
  if (!is.null(p))
    g <- g + geom_text(aes(x = ifelse(.data[[est]] > null, .data[[est]], null), y = .data$.y,
                           label = fmt_p(.data[[p]])), parse = TRUE, hjust = -0.1, size = pt_size(6))
  g
}

# Volcano: grey cloud, hits in the two direction colours, counts as coloured
# corner labels instead of a legend, names on leaders in the hit colour.
fig_volcano <- function(df, lfc = "log2FoldChange", p = "padj", label = NULL, lfc_cut = 1, p_cut = 0.05,
                        top_n = 10, colours = pub_accent, rasterise = TRUE,
                        xlab = expression(log[2]~"fold change"), ylab = expression(-log[10]~italic(P)[adj])) {
  d <- as.data.frame(df); d <- d[!is.na(d[[lfc]]) & !is.na(d[[p]]), ]
  y <- -log10(d[[p]]); fm <- max(y[is.finite(y)]); y[!is.finite(y)] <- fm * 1.05; d$.y <- y
  d$.dir <- ifelse(d[[p]] < p_cut & d[[lfc]] >= lfc_cut, "up",
                   ifelse(d[[p]] < p_cut & d[[lfc]] <= -lfc_cut, "down", "ns"))
  pts <- geom_point(size = 0.45, stroke = 0)
  g <- ggplot(d[order(d$.dir != "ns"), ], aes(.data[[lfc]], .data$.y, colour = .data$.dir)) +
    ref_line(xintercept = c(-lfc_cut, lfc_cut)) + ref_line(yintercept = -log10(p_cut)) +
    (if (rasterise) rasterise_pts(pts) else pts) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = pt_size(6),
             label = sprintf("Down %s", format(sum(d$.dir == "down"), big.mark = ","))) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.3, size = pt_size(6),
             label = sprintf("Up %s", format(sum(d$.dir == "up"), big.mark = ","))) +
    scale_colour_manual(values = colours, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +  # room for the corner counts
    labs(x = xlab, y = ylab) + theme_pub()
  if (!is.null(label) && top_n > 0) {
    s <- d[d$.dir != "ns", ]; s <- utils::head(s[order(-(s$.y * abs(s[[lfc]]))), ], top_n)
    g <- g + label_points(s, lfc, ".y", label, fontface = "italic",
                          ylim = c(-Inf, max(d$.y)))  # keep names out of the count corners
  }
  g
}

# GO/pathway: lollipop on fold enrichment, colour = -log10 padj (v1 recipe).
fig_go_dot <- function(df, term = "term", x = "fold_enrichment", size = "count", colour = "padj",
                       top_n = 15, xlab = "Fold enrichment", wrap = 32) {
  d <- utils::head(as.data.frame(df)[order(df[[colour]]), ], top_n)
  d[[term]] <- vapply(as.character(d[[term]]), function(s) paste(strwrap(s, wrap), collapse = "\n"), "")
  d[[term]] <- factor(d[[term]], levels = d[[term]][order(d[[x]])])
  d$.s <- -log10(d[[colour]])
  ggplot(d, aes(.data[[x]], .data[[term]])) +
    geom_segment(aes(xend = 0, yend = .data[[term]]), linewidth = lw(0.5), colour = "grey80") +
    geom_point(aes(size = .data[[size]], colour = .data$.s)) +
    scale_colour_gradientn(colours = c("#C6DBEF", "#4292C6", "#08306B"), name = expression(-log[10]~italic(P)[adj])) +
    scale_size_continuous(range = c(0.6, 2.6), name = "Genes") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08)), limits = c(0, NA)) +
    labs(x = xlab, y = NULL) + theme_pub() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())
}


# ==== Signal along a coordinate =================================================

# Track stack: unframed filled areas, track name at the left, data range as
# "[0–max]" top-left, optional highlight band through every track, gene model
# last. `df`: long table (track, x, y). `genes`: data.frame(start, end, type[, gene]).
fig_tracks <- function(df, x = "pos", y = "signal", track = "track", colours = NULL, genes = NULL,
                       highlight = NULL, shared_y = FALSE, xlab = "Position", heights = NULL) {
  d <- as.data.frame(df); d[[track]] <- .fac(d[[track]]); lv <- levels(d[[track]])
  cols <- .cols(lv, colours); xr <- range(d[[x]])
  ymx <- tapply(d[[y]], d[[track]], max, na.rm = TRUE); if (shared_y) ymx[] <- max(ymx)
  one <- function(t) {
    s <- d[d[[track]] == t, ]
    ggplot(s, aes(.data[[x]], .data[[y]])) +
      (if (!is.null(highlight)) band_x(highlight[1], highlight[2]) else NULL) +
      geom_area(fill = cols[[t]], colour = NA) +
      annotate("text", x = xr[1], y = ymx[[t]], label = sprintf("[0–%s]", signif(ymx[[t]], 2)),
               hjust = 0, vjust = 1, size = pt_size(6)) +
      scale_x_continuous(limits = xr, expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, ymx[[t]]), expand = c(0, 0)) + labs(x = NULL, y = t) + theme_pub() +
      theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
            axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
            plot.margin = margin(0.4, 1.5, 0.4, 1.5, "mm"))
  }
  ps <- lapply(lv, one)
  if (!is.null(genes)) {
    gm <- ggplot() + (if (!is.null(highlight)) band_x(highlight[1], highlight[2]) else NULL) +
      gene_layers(genes, name = if ("gene" %in% names(genes)) genes$gene[1] else NULL) +
      scale_x_continuous(limits = xr, expand = c(0, 0)) + scale_y_continuous(limits = c(-1.1, 0.5)) +
      labs(x = xlab, y = NULL) + theme_pub() +
      theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
    ps <- c(ps, list(gm))
  } else ps[[length(ps)]] <- ps[[length(ps)]] + labs(x = xlab) +
      theme(axis.line.x = element_line(linewidth = lw(0.5)), axis.text.x = element_text(),
            axis.ticks.x = element_line(linewidth = lw(0.5)))
  patchwork::wrap_plots(ps, ncol = 1, heights = heights %||% c(rep(1, length(lv)), if (!is.null(genes)) 0.8))
}
fig_coverage <- function(df, x = "pos", y = "cov", group = NULL, ...) {  # v1 name
  df <- as.data.frame(df)
  if (is.null(group)) { df$.t <- "signal"; group <- ".t" }
  fig_tracks(df, x = x, y = y, track = group, ...)
}

# Metagene on scaled coordinates (one unit per region) with the transcript drawn
# as the x axis. null = column with a background curve, drawn as a grey fill.
fig_metagene <- function(df, x = "position", y = "value", group = NULL, ymin = NULL, ymax = NULL, null = NULL,
                         regions = c("5′UTR", "CDS", "3′UTR"), colours = NULL, ylab = "Density",
                         cartoon = TRUE) {
  d <- as.data.frame(df); nb <- length(regions)
  g <- ggplot(d, if (is.null(group)) aes(.data[[x]], .data[[y]]) else
    aes(.data[[x]], .data[[y]], colour = .data[[group]], fill = .data[[group]]))
  if (!is.null(null)) g <- g + geom_area(aes(y = .data[[null]]), fill = "grey88", colour = NA, position = "identity")
  if (!is.null(ymin)) g <- g + geom_ribbon(aes(ymin = .data[[ymin]], ymax = .data[[ymax]]), alpha = 0.22, colour = NA)
  g <- g + ref_line(xintercept = seq_len(nb - 1)) + geom_line(linewidth = lw(0.9)) +
    scale_x_continuous(limits = c(0, nb), breaks = NULL, expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + labs(x = NULL, y = ylab) + theme_pub() +
    theme(axis.line.x = element_blank(), legend.position = "inside", legend.position.inside = c(0.99, 0.99),
          legend.justification = c(1, 1), legend.title = element_blank())
  if (!is.null(group)) { cols <- .cols(levels(.fac(d[[group]])), colours)
    g <- g + scale_colour_manual(values = cols) + scale_fill_manual(values = cols) }
  if (!cartoon) return(g + scale_x_continuous(breaks = seq_len(nb) - 0.5, labels = regions, expand = c(0, 0)))
  g / tx_cartoon(regions) + patchwork::plot_layout(heights = c(1, 0.13))
}

# Signal heatmap around anchors (tornado): aggregate profile on top, rows sorted
# by signal within each group, white gaps between groups, one-hue ramp, the
# 99th percentile as ceiling so a few hot rows do not flatten the rest.
fig_signal_heatmap <- function(df, row = "id", pos = "pos", value = "signal", group = NULL, colours = NULL,
                               hue = "#08306B", cap_q = 0.99, xlab = "Distance to anchor (kb)", fill_lab = "Signal") {
  d <- as.data.frame(df); if (is.null(group)) { d$.g <- "all"; group <- ".g" }
  d[[group]] <- .fac(d[[group]]); cols <- .cols(levels(d[[group]]), colours)
  o <- tapply(d[[value]], d[[row]], mean); d[[row]] <- factor(d[[row]], levels = names(sort(o)))
  cap <- stats::quantile(d[[value]], cap_q, na.rm = TRUE)
  prof <- stats::aggregate(d[[value]], list(pos = d[[pos]], g = d[[group]]), mean)
  top <- ggplot(prof, aes(.data$pos, .data$x, colour = .data$g)) + geom_line(linewidth = lw(0.9)) +
    scale_colour_manual(values = cols, name = NULL) + scale_x_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = "Mean") + theme_pub() +
    theme(axis.text.x = element_blank(), legend.position = "inside", legend.position.inside = c(1, 1),
          legend.justification = c(1, 1))
  hm <- ggplot(d, aes(.data[[pos]], .data[[row]], fill = pmin(.data[[value]], cap))) +
    rasterise_pts(geom_raster()) +
    facet_grid(rows = vars(.data[[group]]), scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradientn(colours = c("#FFFFFF", tint(hue, 0.55), hue), name = fill_lab) +
    scale_x_continuous(expand = c(0, 0)) + labs(x = xlab, y = NULL) + theme_pub() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          strip.text.y.left = element_text(angle = 90), panel.spacing.y = unit(1, "mm"),
          legend.position = "bottom", legend.key.height = unit(1.6, "mm"), legend.key.width = unit(6, "mm"))
  top / hm + patchwork::plot_layout(heights = c(0.28, 1))
}

# Two matrices in one square, split on the diagonal (A upper, B lower) on one
# shared colour scale with its range printed — the 3D-genome convention, equally
# good for two correlation or co-occurrence matrices.
fig_split_matrix <- function(A, B, labels = c("A", "B"), log = TRUE, colours = PUB_FALL, limits = NULL,
                             fill_lab = "Contacts") {
  n <- nrow(A); g <- expand.grid(i = seq_len(n), j = seq_len(n))
  g$v <- ifelse(g$j >= g$i, A[cbind(g$i, g$j)], B[cbind(g$i, g$j)])
  lim <- limits %||% range(g$v[g$v > 0 | !log], na.rm = TRUE)
  ggplot(g, aes(.data$j, .data$i, fill = .data$v)) + rasterise_pts(geom_raster()) +
    annotate("segment", x = 0.5, y = 0.5, xend = n + 0.5, yend = n + 0.5, colour = "grey40", linewidth = lw(0.3)) +
    annotate("text", x = n, y = 1, label = labels[1], hjust = 1.1, vjust = 1.3, fontface = "bold", size = pt_size(6.5)) +
    annotate("text", x = 1, y = n, label = labels[2], hjust = -0.1, vjust = -0.4, fontface = "bold", size = pt_size(6.5)) +
    scale_fill_gradientn(colours = colours, limits = lim, trans = if (log) "log10" else "identity",
                         labels = if (log) scales::label_log() else waiver(), oob = scales::squish, name = fill_lab) +
    scale_y_reverse(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) + coord_fixed() +
    theme_pub() + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
                        axis.title = element_blank(), legend.key.width = unit(1.8, "mm"))
}


# ==== Cells and embeddings =======================================================

# Embedding with names on the clusters instead of a legend, small rasterised
# points, corner arrows instead of axes; `highlight` greys everything else.
fig_embedding <- function(df, x, y, group, colours = NULL, highlight = NULL, label = TRUE,
                          size = 0.35, axes = c("UMAP1", "UMAP2")) {
  d <- as.data.frame(df); d[[group]] <- .fac(d[[group]]); lv <- levels(d[[group]])
  cols <- .cols(lv, colours); if (!is.null(highlight)) cols[setdiff(lv, highlight)] <- "#E3E3E3"
  d <- d[order(d[[group]] %in% (highlight %||% lv)), ]
  ctr <- stats::aggregate(d[c(x, y)], list(g = d[[group]]), stats::median)
  if (!is.null(highlight)) ctr <- ctr[ctr$g %in% highlight, ]
  g <- ggplot(d, aes(.data[[x]], .data[[y]], colour = .data[[group]])) +
    rasterise_pts(geom_point(size = size, stroke = 0)) +
    scale_colour_manual(values = cols, guide = if (label) "none" else "legend") +
    embed_axes(d[[x]], d[[y]], axes) + theme_pub() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(),
          axis.title = element_blank())
  if (label) g <- g + ggrepel::geom_text_repel(data = ctr, aes(.data[[x]], .data[[y]], label = .data$g),
                                               colour = "black", size = pt_size(6),
                                               fontface = "bold", segment.size = lw(0.3), box.padding = 0.1,
                                               bg.colour = "white", bg.r = 0.12, family = PUB_FONT,
                                               inherit.aes = FALSE, seed = 1)
  g
}

# Marker dot plot: size = fraction expressing, colour = scaled mean (one hue);
# a strip of identity colours replaces the group name legend.
fig_markerdot <- function(df, gene, group, frac, expr, colours = NULL, hue = "#08306B",
                          fill_lab = "Mean expr.", size_lab = "% cells") {
  d <- as.data.frame(df); d[[group]] <- .fac(d[[group]]); d[[gene]] <- .fac(d[[gene]])
  cols <- .cols(levels(d[[group]]), colours)
  key <- data.frame(g = factor(levels(d[[group]]), levels = levels(d[[group]])), col = unname(cols))
  # points first: the first layer fixes x as discrete, then the key strip sits at x = 0.2
  ggplot(d, aes(.data[[gene]], .data[[group]])) +
    geom_point(aes(size = .data[[frac]], colour = .data[[expr]])) +
    geom_tile(data = key, aes(x = 0.2, y = .data$g, fill = I(.data$col)), width = 0.4, height = 0.85, inherit.aes = FALSE) +
    scale_colour_gradientn(colours = c("#EDEDED", tint(hue, 0.45), hue), name = fill_lab) +
    scale_size_area(max_size = 3.2, name = size_lab) + scale_y_discrete(limits = rev) +
    labs(x = NULL, y = NULL) + theme_pub() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
          legend.key.width = unit(2, "mm"))
}

# 100% composition bars: n above each bar, segment names at the right of the
# last bar instead of a legend. pct = TRUE prints percentages (black) inside
# segments ≥ label_min — use light fills if you do.
fig_composition <- function(df, x, fill, colours = NULL, pct = FALSE, label_min = 0.08, xlab = NULL, ylab = "Fraction") {
  d <- as.data.frame(df); d[[x]] <- .fac(d[[x]]); d[[fill]] <- .fac(d[[fill]])
  t <- as.data.frame(table(x = d[[x]], f = d[[fill]])); t$p <- stats::ave(t$Freq, t$x, FUN = function(v) v / sum(v))
  t$f <- factor(t$f, levels = rev(levels(d[[fill]])))
  t <- t[order(t$x, t$f), ]; t$top <- stats::ave(t$p, t$x, FUN = cumsum); t$mid <- t$top - t$p / 2
  cols <- .cols(levels(d[[fill]]), colours); n <- table(d[[x]]); last <- t[t$x == levels(t$x)[nlevels(t$x)], ]
  ggplot(t, aes(.data$x, .data$p, fill = .data$f)) +
    geom_col(width = 0.7, colour = "white", linewidth = lw(0.4)) +
    (if (pct) geom_text(data = t[t$p >= label_min, ], aes(y = 1 - .data$mid, label = sprintf("%.0f", 100 * .data$p)),
                        size = pt_size(6), colour = "black") else NULL) +
    annotate("text", x = seq_along(n), y = 1.03, label = format(as.integer(n), big.mark = ","), vjust = 0, size = pt_size(6)) +
    annotate("text", x = nlevels(t$x) + 0.42, y = 1 - last$mid, label = last$f, hjust = 0,
             colour = "black", size = pt_size(6)) +
    scale_fill_manual(values = cols, guide = "none") +
    scale_y_continuous(breaks = c(0, 0.5, 1), expand = expansion(mult = c(0, 0.08))) +
    coord_cartesian(clip = "off") + labs(x = xlab, y = ylab) + theme_pub() +
    theme(plot.margin = margin(1.5, 14, 1.5, 1.5, "mm"))
}


# ==== Models, benchmarks, time ==================================================

# Benchmark table (funky-heatmap style): circle area = min-max scaled score per
# metric, fill = rank within metric; methods ordered by mean; one method bold.
fig_bench_table <- function(df, method, metric, value, higher_better = TRUE, highlight = NULL,
                            hue = "#08519C", xlab = NULL) {
  d <- as.data.frame(df)
  d$.s <- stats::ave(d[[value]], d[[metric]], FUN = function(v) {
    s <- (v - min(v, na.rm = TRUE)) / diff(range(v, na.rm = TRUE)); if (higher_better) s else 1 - s })
  d$.r <- stats::ave(-d$.s, d[[metric]], FUN = rank)
  ord <- names(sort(tapply(d$.s, d[[method]], mean)))
  d[[method]] <- factor(d[[method]], levels = ord); d[[metric]] <- .fac(d[[metric]])
  ov <- data.frame(m = factor(ord, levels = ord), s = as.numeric(sort(tapply(d$.s, d[[method]], mean))))
  nm <- nlevels(d[[metric]])
  ggplot(d, aes(.data[[metric]], .data[[method]])) +
    geom_point(aes(size = .data$.s, fill = .data$.r), shape = 21, colour = "grey30", stroke = st(0.4)) +
    geom_rect(data = ov, aes(xmin = nm + 0.7, xmax = nm + 0.7 + .data$s * 1.6, ymin = as.numeric(.data$m) - 0.32,
                             ymax = as.numeric(.data$m) + 0.32), fill = "grey55", inherit.aes = FALSE) +
    annotate("text", x = nm + 0.7, y = length(ord) + 0.75, label = "Mean", hjust = 0, size = pt_size(6)) +
    geom_text(data = ov, aes(x = 0.35, y = .data$m, label = .data$m,
                             fontface = ifelse(.data$m %in% highlight, "bold", "plain")),
              hjust = 1, size = pt_size(6), inherit.aes = FALSE) +
    scale_fill_gradientn(colours = c(hue, tint(hue, 0.5), "#F2F2F2"), name = "Rank") +
    scale_size_area(max_size = 3.6, guide = "none") +
    scale_x_discrete(position = "top", expand = expansion(add = c(0.6, 2.4))) +
    scale_y_discrete(expand = expansion(add = c(0.6, 1))) + coord_cartesian(clip = "off") +
    labs(x = xlab, y = NULL) + theme_pub() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(),
          axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0),
          plot.margin = margin(1.5, 1.5, 1.5, 16, "mm"))
}

# Predicted vs observed (or replicate vs replicate): points coloured by local
# density, dashed identity line, r and RMSE in the corner, equal axes.
fig_pred_obs <- function(df, x, y, density = TRUE, xlab = x, ylab = y, colour = "#08519C", size = 0.5) {
  d <- as.data.frame(df); d <- d[is.finite(d[[x]]) & is.finite(d[[y]]), ]
  if (density) {
    k <- MASS::kde2d(d[[x]], d[[y]], n = 100)
    d$.den <- k$z[cbind(findInterval(d[[x]], k$x), findInterval(d[[y]], k$y))]
    d <- d[order(d$.den), ]
  }
  lim <- range(c(d[[x]], d[[y]])); r <- stats::cor(d[[x]], d[[y]]); rmse <- sqrt(mean((d[[x]] - d[[y]])^2))
  pts <- if (density) geom_point(aes(colour = .data$.den), size = size, stroke = 0) else
    geom_point(colour = colour, size = size, stroke = 0, alpha = 0.6)
  ggplot(d, aes(.data[[x]], .data[[y]])) + ref_line(slope = 1) + .maybe_raster(pts, nrow(d), 2000) +
    scale_colour_gradientn(colours = c("#C6DBEF", "#4292C6", "#08306B", "#081D58"), guide = "none") +
    annotate("text", x = lim[1], y = lim[2], hjust = 0, vjust = 1, size = pt_size(6),
             label = sprintf("r = %.2f\nRMSE = %.2f\nn = %s", r, rmse, format(nrow(d), big.mark = ","))) +
    coord_equal(xlim = lim, ylim = lim) + labs(x = xlab, y = ylab) + theme_pub()
}

# Lines over time / dose / model size / threshold, named at their ends, optional
# ribbon and reference. For ROC/PR put the AUC in the group label.
fig_lines <- function(df, x, y, group, colours = NULL, ymin = NULL, ymax = NULL, ref = NULL, points = TRUE,
                      logx = FALSE, xlab = x, ylab = y) {
  d <- as.data.frame(df); d[[group]] <- .fac(d[[group]]); cols <- .cols(levels(d[[group]]), colours)
  g <- ggplot(d, aes(.data[[x]], .data[[y]], colour = .data[[group]], fill = .data[[group]]))
  if (!is.null(ref)) g <- g + ref_line(yintercept = ref)
  if (!is.null(ymin)) g <- g + geom_ribbon(aes(ymin = .data[[ymin]], ymax = .data[[ymax]]), alpha = 0.18, colour = NA)
  g + geom_line(linewidth = lw(0.9)) + (if (points) geom_point(size = 0.9) else NULL) +
    label_ends(d, x, y, group) +
    scale_colour_manual(values = cols, guide = "none") + scale_fill_manual(values = cols, guide = "none") +
    (if (logx) scale_x_log10(expand = expansion(mult = c(0.02, 0.28))) else
       scale_x_continuous(expand = expansion(mult = c(0.02, 0.28)))) +
    labs(x = xlab, y = ylab) + theme_pub()
}

invisible(NULL)
