# Additional file 1 Fig. S5: the clock in disease
# a, age acceleration of every sample by region, cases and controls; b, the 30 clocks scored on the same donors; c, the
# acceleration against the individual m6A index in aging cohorts
# reads source_data/sfig5a, sfig5a_tests, sfig5b and sfig5c (scripts 53 and 54); writes output/figures/sfig5
source("figures/style.R")

STUDIES <- c("Ping 2020", "Merrihew 2023", "AD combined", "Blumenreich 2024", "Umoh 2018", "Langfelder 2016")
wrap <- function(x, width) vapply(x, function(s) paste(strwrap(s, width), collapse = "\n"), "")
# row label as plotmath: m6A with its superscript, the m6Brain clock in bold
math_label <- function(s, bold) {
  parts <- strsplit(gsub(paste0("(", MARK, ")"), "\r\\1\r", s), "\r")[[1]]
  pieces <- lapply(parts, function(p) {
    if (grepl(paste0("^", MARK, "$"), p)) bquote(m^.(substr(p, 2, 2)) * .(as.name(substr(p, 3, 3)))) else p
  })
  label <- as.call(c(as.name("paste"), pieces))
  if (bold) bquote(bold(.(label))) else label
}

# a: cases (brick) and controls (grey) of each region (box, median and quartiles; every donor), with the coefficient of
# the disease in the region (years) and its P above the bracket; Alzheimer's disease and ALS above, Parkinson's below
s <- read_panel("sfig5a")
tests <- read_panel("sfig5a_tests")
for (d in list(s, tests)) d[, strip := paste0(study, "\n", wrap(region, 24))]
tests[, `:=`(row = fifelse(disease == "PD", 2, 1), text = signed(est, 1))]
s[, group := factor(fifelse(group == "Case", disease, "Control"), c("Control", "AD", "PD", "ALS"))]
BRACKET <- 33
box_row <- function(r) {
  t <- tests[row == r][order(match(study, STUDIES))][, strip := factor(strip, strip)]
  d <- s[strip %in% t$strip][, strip := factor(strip, levels(t$strip))]
  fig_box(d, "group", "acceleration", colours = c(Control = GREY, AD = HUMAN, PD = HUMAN, ALS = HUMAN),
          n_labels = FALSE, ref = 0, ylab = "Age acceleration (years)") +
    geom_path(data = t[, .(x = c(1, 1, 2, 2), y = BRACKET - c(2, 0, 0, 2)), by = strip], aes(x, y),
              inherit.aes = FALSE, linewidth = lw(0.5)) +
    geom_text(data = t, aes(x = 1.5, y = BRACKET + 2.8, label = sprintf("'%s, '*%s", text, fmt_p(p))), parse = TRUE,
              vjust = 0, size = pt_size(6), inherit.aes = FALSE) +
    facet_wrap(~strip, nrow = 1, scales = "free_x") +
    scale_y_continuous(labels = MINUS_LABELS, breaks = c(-20, 0, 20)) +
    coord_cartesian(ylim = c(-23, 55)) +
    theme_box() +
    theme(strip.text = element_text(size = 6, lineheight = 0.9))
}
a <- box_row(1) / box_row(2)

# b: the 30 clocks on the same donors, coefficient of the disease on each score (colour; SD of the score, HD per 100
# CAG repeats) and the rank of the clock in the study (number; bold, P < 0.05); rows ordered by the AD combined
# coefficient, the m6Brain clock in bold
r <- read_panel("sfig5b")[unit == "SD"]
r[, `:=`(column = fifelse(study == "AD combined", "Combined", sub(" ", "\n", study)),
         block = factor(disease, c("AD", "PD", "ALS", "HD")))]
r[, column := factor(column, unique(column[order(match(study, STUDIES))]))]
rows <- r[study == "AD combined"][order(est), label]
r[, row := factor(label, rows)]
strip_fill <- c(rep(tint(HUMAN, 0.45), 3), tint(MOUSE, 0.45))
b <- ggplot(r, aes(column, row, fill = est)) +
  geom_tile(colour = "white", linewidth = lw(0.6)) +
  geom_text(aes(label = rank, fontface = fifelse(p < 0.05, "bold", "plain")), size = pt_size(6)) +
  ggh4x::facet_grid2(. ~ block, scales = "free_x", space = "free_x",
                     strip = ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = strip_fill,
                                                                                      colour = strip_fill))) +
  scale_fill_signed(1.5, "Disease coefficient\n(SD of the score;\nHD per 100 CAG)", light = TRUE) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(labels = as.expression(Map(math_label, rows, rows == "m6Brain clock")), expand = c(0, 0)) +
  theme_tile() +
  theme(legend.position = "inside", legend.position.inside = c(0, 1), legend.justification = c(1, 0),
        legend.margin = margin(r = 3, unit = "mm"), strip.placement = "outside",
        strip.text.x = element_text(size = 6), panel.spacing.x = unit(1, "mm"),
        axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, lineheight = 0.9))

# c: acceleration per SD of the age-adjusted m6A index in every cohort, sorted (bar and 95% CI; faded when the CI
# includes 0) and pooled by random effects (dashed line, with its P)
SETS <- c(mouse = "Held-out mouse\n(months)", external = "External mouse\n(months)", human = "Held-out human\n(years)")
COL <- c(mouse = MOUSE, external = EXTERNAL, human = HUMAN)
c_ <- read_panel("sfig5c")
c_[, `:=`(panel = factor(SETS[set], SETS), label = sprintf("%s, n = %d", cohort, n))]
bars <- c_[kind == "cohort"][order(match(set, names(SETS)), -est)]
bars[, `:=`(row = factor(paste(set, label), unique(paste(set, label))),
            fill = fifelse(ci_low < 0 & ci_high > 0, tint(COL[set], 0.6), COL[set]))]
pooled <- c_[kind == "pooled"][bars[, .(top = .N), by = panel], on = "panel"]
pooled[, text := sprintf("'Pooled %s, '*%s", number(est), fmt_p(p))]
cc <- ggplot(bars, aes(est, row)) +
  geom_vline(xintercept = 0, linewidth = lw(0.5)) +
  geom_col(aes(fill = I(fill)), width = 0.72) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, linewidth = lw(0.5), colour = EDGE) +
  geom_segment(data = pooled, aes(x = est, xend = est, y = 0.4, yend = top + 0.6), linetype = "22",
               linewidth = lw(0.75), inherit.aes = FALSE) +
  geom_text(data = pooled, aes(x = est, y = top + 0.95, label = text), parse = TRUE, hjust = 1, vjust = 0, size = pt_size(6),
            inherit.aes = FALSE) +
  ggh4x::facet_grid2(panel ~ ., scales = "free", space = "free_y", independent = "x") +
  scale_y_discrete(labels = function(x) sub("^[a-z]+ ", "", x), expand = expansion(add = c(0.5, 2.3))) +
  scale_x_continuous(labels = MINUS_LABELS, breaks = scales::breaks_pretty(4)) +
  labs(x = "Acceleration per SD of the age-adjusted m6A index", y = NULL) +
  theme_box() +
  theme(strip.text.y = element_text(size = 6, lineheight = 0.9), axis.ticks.y = element_blank())

fig <- fig_compose(a, patchwork::wrap_elements(full = b), patchwork::wrap_elements(full = cc), design = "AA\nBC",
                   widths = c(93, 90), heights = c(74, 96))
save_fig(fig, "sfig5", 170)
