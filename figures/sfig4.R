# Additional file 1 Fig. S4: details of the clock; a, r of every weight set in every cohort and pooled per set; b, the
# label-free mode; c, unseen proteins; d, the error across cohorts against the training-mean null
# reads source_data/sfig4a-d, sfig4b_null and sfig4c_null (scripts 52 and 55); writes output/figures/sfig4
source("figures/style.R")

SETS <- c("held-out|mouse" = "Held-out mouse", "external|mouse" = "External mouse", "held-out|human" = "Held-out human")
GROUPS <- c("m6Brain clock and references" = "m6Brain clock\nand references",
            "Proteomic clocks" = "Proteomic\nclocks", "Transcriptome clocks" = "Transcriptome\nclocks",
            "Signatures and half-lives" = "Signatures and\nhalf-lives",
            "Label-free and gene sets" = "Label-free and\ngene sets")
tag_set <- function(d) d[, set := factor(SETS[paste(set, species, sep = "|")], SETS)]
# row label as plotmath: m6A with its superscript, the m6Brain clock in bold
math_label <- function(s, bold) {
  parts <- strsplit(gsub(paste0("(", MARK, ")"), "\r\\1\r", s), "\r")[[1]]
  pieces <- lapply(parts, function(p) {
    if (grepl(paste0("^", MARK, "$"), p)) bquote(m^.(substr(p, 2, 2)) * .(as.name(substr(p, 3, 3)))) else p
  })
  # every row holds a superscript, seen or not, so that all rows sit on one baseline
  label <- as.call(c(as.name("paste"), quote(phantom(""^6)), pieces))
  if (bold) bquote(bold(.(label))) else label
}
# arrows dropping onto the axis of a histogram from the height `top`, each named left of its tail
drop_arrows <- function(d) {
  if (!"end" %in% names(d)) d <- copy(d)[, end := 0]
  list(geom_segment(data = d, aes(x = x, xend = x, y = top, yend = end, colour = I(colour)), linewidth = lw(1),
                    arrow = grid::arrow(length = unit(1.3, "mm"), type = "closed")),
       geom_text(data = d, aes(x = x - 0.012, y = top, label = label, vjust = vjust), hjust = 1, size = pt_size(6),
                 lineheight = 0.9))
}
# one break at zero and one at the leading digit of the largest count
count_breaks <- function(room) function(limits) {
  e <- 10^floor(log10(limits[2] / room))
  c(0, e * floor(limits[2] / room / e))
}

# a: Pearson r of every weight set with age in every cohort (columns ordered by the m6Brain clock within the set) and
# pooled over the cohorts of each set; rows grouped as in Fig. 5c, then the label-free weight sets, ordered by the
# mean pooled r; bold, the m6Brain clock
d <- tag_set(read_panel("sfig4a"))
d[, pooled := kind == "pooled"]
order_of <- d[weights == "m6brain_clock" & !pooled][order(set, -r), cohort]
d[, mean_r := mean(r[pooled]), by = label]
rows <- unique(d[order(match(group, names(GROUPS)), -(weights == "m6brain_clock"), -mean_r),
                 .(label, group, bold = weights == "m6brain_clock")])
rows[, y := -(seq_len(.N) + 0.5 * (match(group, names(GROUPS)) - 1))]
d <- rows[d[, !"group"], on = "label"]
d[, `:=`(block = factor(fifelse(pooled, "Pooled", as.character(set)), c(SETS, "Pooled")),
         column = factor(fifelse(pooled, as.character(set), cohort), c(unique(order_of), SETS)))]
groups <- rows[, .(top = max(y) + 0.45, bottom = min(y) - 0.45, label = GROUPS[group[1]]), by = group]
groups[, `:=`(middle = (top + bottom) / 2, block = factor("Pooled", c(SETS, "Pooled")))]
rows <- rows[order(y)]
strip_fill <- c(tint(unname(SET), 0.45), LIGHT)
a <- ggplot(d, aes(column, y)) +
  geom_tile(aes(fill = r), colour = "white", linewidth = lw(0.6)) +
  geom_text(aes(label = sub("^(−?)0[.]", "\\1.", number(r)), fontface = fifelse(bold, "bold", "plain")),
            size = pt_size(5)) +
  geom_segment(data = groups, aes(x = 3.65, xend = 3.65, y = top, yend = bottom), linewidth = lw(0.5)) +
  geom_text(data = groups, aes(x = 3.8, y = middle, label = label), hjust = 0, size = pt_size(6), lineheight = 0.9) +
  ggh4x::facet_grid2(. ~ block, scales = "free_x", space = "free_x",
                     strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = strip_fill,
                                                                                      colour = strip_fill))) +
  scale_fill_signed(1, "Pearson r of the score with age", light = TRUE) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_continuous(breaks = rows$y, labels = as.expression(Map(math_label, rows$label, rows$bold)),
                     expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_tile() +
  theme(legend.position = "inside", legend.position.inside = c(0, 1), legend.justification = c(1, 0),
        legend.margin = margin(r = 3, unit = "mm"), strip.placement = "outside",
        strip.text.x = element_text(size = 6), panel.spacing.x = unit(1, "mm"),
        axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5),
        plot.margin = margin(1.5, 16, 1.5, 1.5, "mm"))

# b: label-free mode, pooled r of 1,000 random weights (histogram) and of the m6Brain predictions used as weights
# (coloured arrow, with its P against the random weights) and of the three comparison weights (grey arrows) per set
BW <- 0.03
ROOM <- 4.3 # height of a facet in units of the tallest bar
lf <- tag_set(read_panel("sfig4b"))
random <- tag_set(read_panel("sfig4b_null"))
peak <- random[, .(peak = max(hist(r, breaks = seq(-1, 1, BW), plot = FALSE)$counts)), by = set]
lf <- peak[lf, on = "set"][order(weights == "model")]
lf[, `:=`(x = r, top = peak * (1.5 + 0.62 * (frank(r, ties.method = "first") - 1)), vjust = 0.5,
          colour = fifelse(weights == "model", SET[as.character(set)], EDGE)), by = set]
b <- ggplot(random, aes(r)) +
  geom_histogram(binwidth = BW, boundary = 0, fill = GREY, colour = NA) +
  geom_hline(yintercept = 0, linewidth = lw(0.5)) +
  drop_arrows(lf) +
  geom_text(data = lf[weights == "model"], aes(x = r + 0.03, y = top, label = fmt_p(signif(p_random, 2))), hjust = 0,
            parse = TRUE, size = pt_size(6)) +
  geom_text(data = peak, aes(x = -Inf, y = Inf, label = as.character(set)), hjust = -0.06, vjust = 1.4,
            size = pt_size(6), inherit.aes = FALSE) +
  geom_blank(data = peak, aes(y = ROOM * peak), inherit.aes = FALSE) +
  facet_grid(set ~ ., scales = "free_y") +
  scale_x_continuous(labels = MINUS_LABELS, breaks = c(-0.5, 0, 0.5, 1)) +
  scale_y_continuous(breaks = count_breaks(ROOM), expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(-0.7, 1.15)) +
  labs(x = "Pooled Pearson r with age", y = "Random draws") +
  theme_box() +
  theme(strip.text = element_blank())

# c: unseen proteins, pooled r of 1,000 permutations of the weights among the same genes (histogram) and r of the
# label-free score from the weights of the genes without a protein age effect in any discovery cohort, in the two
# external deep-DIA cohorts and pooled (arrows, with the P of the pooled r against the permutations)
u <- read_panel("sfig4c")
perm <- read_panel("sfig4c_null")
peak <- max(hist(perm$r, breaks = seq(-1, 1, BW), plot = FALSE)$counts)
u[, `:=`(x = r, top = peak * (1.5 + 0.95 * (frank(r, ties.method = "first") - 1)), pooled = kind == "pooled",
         end = fifelse(cohort == "Chen 2025", peak * 0.5, 0))]
u[, `:=`(colour = fifelse(pooled, EXTERNAL, EDGE), vjust = fifelse(pooled, 0, 0.5),
         label = fifelse(pooled, sprintf("Pooled, %s genes", thousands(n_genes)),
                         sprintf("%s, n = %d,\n%d genes", cohort, n, n_genes)))]
cc <- ggplot(perm, aes(r)) +
  geom_histogram(binwidth = BW, boundary = 0, fill = GREY, colour = NA) +
  geom_hline(yintercept = 0, linewidth = lw(0.5)) +
  drop_arrows(u) +
  geom_text(data = u[pooled == TRUE], aes(x = r - 0.012, y = top, label = fmt_p(signif(p_random, 2))), hjust = 1,
            vjust = 1.3, parse = TRUE, size = pt_size(6)) +
  scale_x_continuous(labels = MINUS_LABELS, breaks = c(-0.2, 0, 0.2, 0.4)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  coord_cartesian(xlim = c(-0.42, 0.45)) +
  labs(x = "r of the unseen-gene score", y = "Permutations") +
  theme_pub()

# d: MAE of the m6Brain clock (arrow head) against the training-mean null (arrow tail), one column per cohort ordered by
# the offset of its mean age from the mean age of its training samples (printed after the name)
m <- tag_set(read_panel("sfig4d"))
PANELS <- c("Mouse (months)", "Human (years)")
m[, panel := factor(PANELS[(species == "human") + 1], PANELS)]
m <- m[order(panel, age_offset)]
m[, column := factor(sprintf("%s (%s)", cohort, signed(age_offset, 1)), sprintf("%s (%s)", cohort, signed(age_offset, 1)))]
dd <- ggplot(m, aes(column)) +
  geom_segment(data = m[mae != mae_null], aes(xend = column, y = mae_null, yend = mae, colour = set), linewidth = lw(1),
               arrow = grid::arrow(length = unit(1.1, "mm"), type = "closed")) +
  geom_point(data = m[mae == mae_null], aes(y = mae, colour = set), size = 1) +
  ggh4x::facet_grid2(. ~ panel, scales = "free", space = "free_x", independent = "y") +
  scale_colour_manual(values = SET, guide = "none") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)), breaks = scales::breaks_pretty(4)) +
  labs(x = NULL, y = "MAE (null to clock)") +
  theme_box() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5), panel.spacing.x = unit(5, "mm"))

# taller than a journal figure: this page of Additional file 1 holds a table of 38 rows
fig <- fig_panels(a, b, cc, dd, design = "AAAAAAAAAAAAAAAAAA\nBBBBBBBCCCCDDDDDDD", heights = c(142, 73))
save_fig(fig, "sfig4", 215)
