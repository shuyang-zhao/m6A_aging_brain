# Figure 6: the m6Brain clock in neurodegeneration
# reads source_data/fig6a-f (script 53); writes output/figures/fig6
source("figures/style.R")

# a: age acceleration of cases (brick) and controls (grey) by study (box, median and quartiles; dots, donors), with the
# coefficient of the disease adjusted for age, sex and batch and its P; above, the two AD studies combined
a <- read_panel("fig6a")
tests <- read_panel("fig6a_tests")
STUDIES <- c("Ping 2020", "Merrihew 2023", "Blumenreich 2024", "Umoh 2018")
a[, `:=`(i = match(study, STUDIES), case = group != "Control")]
a[, x := (i - 1) * 3.7 + 1 + 1.35 * case]
heads <- tests[study %in% STUDIES, .(study, est, p)][, i := match(study, STUDIES)]
heads[, x := (i - 1) * 3.7 + 1.675]
ad <- tests[study == "AD combined"]
ticks <- unique(a[, .(x, group)])[order(x)]
a <- ggplot(a, aes(x, acceleration)) +
  ref_line(yintercept = 0) +
  geom_boxplot(aes(group = x, fill = case, colour = case), width = 0.8, outlier.shape = NA, linewidth = lw(0.6), staplewidth = 0) +
  ggbeeswarm::geom_quasirandom(aes(group = x, colour = case), size = 0.85, stroke = 0, width = 0.22) +
  geom_text(data = heads, aes(x, 32.6, label = study), size = pt_size(6)) +
  geom_text(data = heads, aes(x, 28.8, label = paste0(signed(est, 1), " years")), size = pt_size(6)) +
  geom_text(data = heads, aes(x, 25, label = fmt_p(p)), parse = TRUE, size = pt_size(6)) +
  annotate("segment", x = 0.6, xend = 6.45, y = 34.8, yend = 34.8, linewidth = lw(0.5)) +
  note(0.6, 35.6, sprintf("list('AD combined %s years', %s)", signed(ad$est, 1), fmt_p(ad$p)), hjust = 0, vjust = 0, parse = TRUE) +
  scale_fill_manual(values = c("TRUE" = tint(HUMAN, 0.7), "FALSE" = tint(GREY, 0.6)), guide = "none") +
  scale_colour_manual(values = c("TRUE" = HUMAN, "FALSE" = "#8C8C8C"), guide = "none") +
  scale_x_continuous(breaks = ticks$x, labels = ticks$group, expand = expansion(add = 0.75)) +
  scale_y_continuous(labels = MINUS_LABELS, breaks = seq(-10, 20, 10)) +
  coord_cartesian(ylim = c(-13, 42), clip = "off") +
  labs(x = NULL, y = "Age acceleration (years)") +
  theme_pub()

# b: HD knock-in mice at 6 months: acceleration of each mouse by CAG length of the genotype (dots; grey, 7 and 20
# repeats), genotype means (bars) and the least-squares line
b <- read_panel("fig6b")
hd <- tests[disease == "HD"]
means <- b[, .(acceleration = mean(acceleration)), by = cag]
b <- ggplot(b, aes(cag, acceleration)) +
  ref_line(yintercept = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = shade(MOUSE, 0.3), fill = tint(MOUSE, 0.7), linewidth = lw(0.8)) +
  geom_point(aes(colour = cag >= 36), size = 0.8, stroke = 0, position = position_jitter(width = 3, height = 0, seed = 1)) +
  geom_segment(data = means, aes(x = cag - 7, xend = cag + 7, yend = acceleration), linewidth = lw(0.9)) +
  note_corner(sprintf("%s months per 100 CAG", signed(hd$est))) +
  note_corner(fmt_p(hd$p), vjust = 2.6, parse = TRUE) +
  scale_colour_manual(values = c("TRUE" = MOUSE, "FALSE" = "#8C8C8C"), guide = "none") +
  scale_x_continuous(breaks = c(0, 100, 200), expand = expansion(add = 10)) +
  scale_y_continuous(labels = MINUS_LABELS, expand = expansion(mult = c(0.04, 0.26))) +
  labs(x = "CAG repeats", y = "Age acceleration (months)") +
  theme_pub()

# c: acceleration by Braak stage and CERAD score in the two AD cohorts (dots, donors; line, least squares), with the
# slope adjusted for age, sex and batch and its P
c_ <- read_panel("fig6c")
slopes <- read_panel("fig6c_tests")
COHORT <- c("Ping 2020" = HUMAN, "Merrihew 2023" = shade(HUMAN, 0.45))
severity <- function(column, name, ylab = NULL) {
  s <- slopes[contrast == paste("per", name)]
  s[, text := signed(est, 1)]
  ggplot(c_, aes(.data[[column]], acceleration, colour = study)) +
    ref_line(yintercept = 0) +
    geom_point(size = 0.8, stroke = 0, position = position_jitter(width = 0.09, height = 0, seed = 1)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = lw(0.8)) +
    key_rows(paste(sub(" .*", "", s$study), s$text), COHORT[s$study], x = -0.1, y = 38.5, step = 4.8, width = max(c_[[column]]) * 0.06) +
    scale_colour_manual(values = COHORT, guide = "none") +
    scale_x_continuous(breaks = 0:6) +
    scale_y_continuous(labels = MINUS_LABELS, breaks = seq(-10, 30, 10)) +
    coord_cartesian(ylim = c(-13, 41)) +
    labs(x = name, y = ylab) +
    theme_pub()
}
c_ <- severity("braak", "Braak stage", "Age acceleration (years)") + severity("cerad", "CERAD score")

# d: the acceleration splits exactly over groups of proteins: each bar adds to the running total; share of the
# acceleration and of the summed absolute weights on the right
d <- read_panel("fig6d")
parts <- d[partition == "cell source x direction"]
parts[, c("cell", "direction") := tstrsplit(group, " | ", fixed = TRUE)]
parts[, cell := fifelse(cell == "no cell-type data", "No cell-type data", paste0(toupper(substr(cell, 1, 1)), substring(cell, 2)))]
waterfall <- function(name, xlab, label) {
  w <- parts[study == name][order(-est)]
  total <- d[study == name & partition == "total", est]
  w[, `:=`(row = .I, end = cumsum(est))][, start := end - est]
  rows <- nrow(w) + 1
  span <- range(c(0, w$start, w$end))
  right <- span[2] + diff(span) * c(0.2, 0.42)
  ggplot(w) +
    geom_vline(xintercept = 0, linewidth = lw(0.5)) +
    geom_segment(aes(x = end, xend = end, y = row, yend = row + 1), colour = CI, linewidth = lw(0.3)) +
    geom_rect(aes(xmin = pmin(start, end), xmax = pmax(start, end), ymin = row - 0.38, ymax = row + 0.38, fill = cell)) +
    annotate("rect", xmin = 0, xmax = total, ymin = rows - 0.38, ymax = rows + 0.38, fill = EDGE) +
    geom_text(aes(right[1], row, label = MINUS(sprintf("%.0f", 100 * share_of_total))), hjust = 1, size = pt_size(6)) +
    geom_text(aes(right[2], row, label = sprintf("%.0f", 100 * weight_share)), hjust = 1, size = pt_size(6)) +
    note(right, 0.42, c("% of\nacceleration", "% of\nweights"), hjust = 1, vjust = 0, lineheight = 0.9) +
    note(right[1], rows, label, hjust = 1) +
    scale_fill_manual(values = CELL, guide = "none") +
    scale_y_reverse(breaks = seq_len(rows), labels = c(paste0(w$cell, ", ", w$direction), "Acceleration"),
                    expand = expansion(add = c(0.6, 2.2))) +
    scale_x_continuous(labels = MINUS_LABELS, breaks = scales::breaks_pretty(4)(span), expand = expansion(mult = c(0.03, 0.02))) +
    coord_cartesian(xlim = c(span[1], right[2]), clip = "off") +
    labs(x = xlab, y = NULL) +
    theme_rows()
}
total_of <- function(name) d[study == name & partition == "total", est]
d_plot <- waterfall("AD combined", "AD, years", sprintf("%s", signed(total_of("AD combined"), 1))) +
  waterfall("Langfelder 2016", "HD, months per 100 CAG", sprintf("%s", signed(total_of("Langfelder 2016"), 2)))

# e: acceleration against the m6A index of every donor (cases and expanded alleles coloured, controls grey), Pearson r
# over the donors of a study
e <- read_panel("fig6e")
r <- read_panel("fig6e_tests")[donors == "all"]
ORDER <- c("Ping 2020", "Merrihew 2023", "Blumenreich 2024", "Umoh 2018", "Langfelder 2016")
STRIPS <- c("AD\nPing 2020", "AD\nMerrihew 2023", "PD\nBlumenreich 2024", "ALS\nUmoh 2018", "HD (months)\nLangfelder 2016")
e[, `:=`(strip = factor(STRIPS[match(study, ORDER)], STRIPS),
         mark = fifelse(group == "Control", "Control", fifelse(disease_name == "HD", "Mouse", "Human")))]
r[, `:=`(strip = factor(STRIPS[match(study, ORDER)], STRIPS), text = lab_stats(r), p_text = fmt_p(p, digits = 1)), by = study]
e <- ggplot(e, aes(m6a_index, acceleration, colour = mark)) +
  geom_point(size = 0.9, stroke = 0) +
  geom_text(data = r, aes(x = Inf, y = Inf, label = text), parse = TRUE, hjust = 1.06, vjust = 1.3, size = pt_size(6), colour = "black") +
  geom_text(data = r, aes(x = Inf, y = Inf, label = p_text), parse = TRUE, hjust = 1.06, vjust = 2.7, size = pt_size(6), colour = "black") +
  facet_wrap(~strip, nrow = 1, scales = "free") +
  scale_colour_manual(values = c(Control = "#8C8C8C", Human = HUMAN, Mouse = MOUSE), guide = "none") +
  scale_x_continuous(labels = MINUS_LABELS, breaks = scales::breaks_pretty(3)) +
  scale_y_continuous(labels = MINUS_LABELS, expand = expansion(mult = c(0.05, 0.62))) +
  labs(x = "m6A index (SD)", y = "Age acceleration (years)") +
  theme_box() +
  theme(strip.text = element_text(size = 6))

# f: coefficient of the disease on the m6A index (SD; HD per 100 CAG repeats; bars, 95% CI lines; faded, P >= 0.05)
f <- read_panel("fig6f")
ROWS <- c("Ping 2020" = "AD, Ping 2020", "Merrihew 2023" = "AD, Merrihew 2023", "AD combined" = "AD combined",
          "Blumenreich 2024" = "PD", "Umoh 2018" = "ALS", "Langfelder 2016" = "HD per 100 CAG")
f[, `:=`(row = factor(ROWS[study], rev(ROWS)), fill = fifelse(disease == "HD", MOUSE, HUMAN), solid = p < 0.05)]
f <- ggplot(f, aes(est, row)) +
  geom_col(aes(fill = I(fill), alpha = solid), width = 0.68, orientation = "y") +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = row), colour = EDGE, linewidth = lw(0.4)) +
  geom_vline(xintercept = 0, linewidth = lw(0.5)) +
  geom_text(aes(x = 0.1, label = fmt_p(p)), parse = TRUE, hjust = 0, size = pt_size(6)) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4), guide = "none") +
  scale_x_continuous(labels = MINUS_LABELS, breaks = c(-2, -1, 0)) +
  coord_cartesian(xlim = c(-2.3, 2.5)) +
  labs(x = "Effect on the m6A index (SD)", y = NULL) +
  theme_rows()

# g: model
transcript <- function(x, y, length, sites, reader) {
  flags <- data.table(x = x + length * sites)
  list(
    annotate("rect", xmin = x, xmax = x + length, ymin = y - 0.55, ymax = y + 0.55, fill = tint(M6BRAIN, 0.75)),
    annotate("segment", x = flags$x, xend = flags$x, y = y + 0.55, yend = y + 2.2, colour = M6BRAIN, linewidth = lw(0.5)),
    annotate("point", x = flags$x, y = y + 2.4, colour = M6BRAIN, size = 0.8),
    if (reader) annotate("point", x = flags$x[c(2, 4)], y = y + 3.7, shape = 21, fill = MOUSE, colour = EDGE, size = 2, stroke = st(0.3))
  )
}
arrow <- function(x, xend, y, yend = y) {
  annotate("segment", x = x, xend = xend, y = y, yend = yend, linewidth = lw(0.6), colour = EDGE,
           arrow = grid::arrow(length = unit(1, "mm"), type = "closed"))
}
bars <- function(x, y, heights, fill) {
  n <- seq_along(heights)
  annotate("rect", xmin = x + (n - 1) * 2.1, xmax = x + (n - 1) * 2.1 + 1.6, ymin = y, ymax = y + heights, fill = fill)
}
g <- canvas(183, 30) +
  transcript(3, 17, 34, c(0.12, 0.3, 0.5, 0.68, 0.82, 0.93), TRUE) +
  transcript(3, 6, 34, 0.55, FALSE) +
  note(3, 26, "Many m6A sites", hjust = 0) +
  note(3, 11.2, "Few m6A sites", hjust = 0) +
  annotate("point", x = 40.5, y = 26, shape = 21, fill = MOUSE, colour = EDGE, size = 2, stroke = st(0.3)) +
  note(42.5, 26, "Readers (mouse)", hjust = 0) +
  arrow(38.5, 44, 17.5) +
  arrow(38.5, 44, 6.5) +
  bars(45.5, 14, c(7, 5.6, 4.2, 2.8), MOUSE) +
  bars(45.5, 3.4, c(6, 6, 6, 6), GREY) +
  note(55, 17.5, "Protein\nfalls with age", hjust = 0, lineheight = 0.95) +
  note(55, 6.4, "Protein\nkept", hjust = 0, lineheight = 0.95) +
  annotate("rect", xmin = 74, xmax = 114, ymin = 1, ymax = 29, fill = tint(M6BRAIN, 0.86)) +
  note(94, 26, "m6Brain") +
  note(94, 21.5, "Predicts the age effect of\nevery gene from its transcript", lineheight = 0.95) +
  arrow(94, 94, 17.6, 13.6) +
  note(94, 8.5, "Guides a proteomic clock\ncalibrated on brain tissue", lineheight = 0.95) +
  note(94, 3.2, "Reads age in 25 cohorts") +
  arrow(67.5, 73, 12) +
  arrow(115, 121, 12) +
  note(123, 26, "Brains the clock calls older", hjust = 0) +
  annotate("point", x = 124.5, y = 21.2, shape = 15, size = 1.8, colour = HUMAN) +
  note(127, 21.2, "AD +10.2 years, HD +2.45 months per 100 CAG", hjust = 0) +
  annotate("point", x = 124.5, y = 16.2, shape = 15, size = 1.8, colour = CELL[["Neuron"]]) +
  note(127, 16.2, "Neuronal proteins lost in AD and HD", hjust = 0) +
  annotate("point", x = 124.5, y = 11.2, shape = 15, size = 1.8, colour = CELL[["Astrocyte"]]) +
  note(127, 11.2, "Astrocyte proteins gained in AD", hjust = 0) +
  annotate("point", x = 124.5, y = 6.2, shape = 15, size = 1.8, colour = M6BRAIN) +
  note(127, 6.2, "m6A-marked proteins lost", hjust = 0)

fig <- fig_panels(a, b, c_, d_plot, e, f, g, design = "AAAAAAAAABBBBBCCCCCC\nDDDDDDDDDDDDDDDDDDDD\nEEEEEEEEEEEEEEFFFFFF\nGGGGGGGGGGGGGGGGGGGG",
                  heights = c(46, 53, 39, 30))
save_fig(fig, "fig6", 170)
