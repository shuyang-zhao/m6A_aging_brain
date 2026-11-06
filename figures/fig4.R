# Figure 4: reader levels and the m6A index
# reads source_data/fig4a-f (scripts 30, 33, 38 and 54); writes output/figures/fig4
source("figures/style.R")

# a: m6A index of every donor against age, a least-squares line per cohort; pooled Pearson r with its P
a <- read_panel("fig4a_donors")
pooled <- read_panel("fig4a")[kind == "pooled"]
index_panel <- function(set, hue, name, xlab, limits, ylab = NULL) {
  d <- a[species == set]
  fit <- pooled[species == set]
  ggplot(d, aes(age, index)) +
    ref_line(yintercept = 0) +
    geom_point(colour = tint(hue, 0.3), size = 0.75, stroke = 0, position = position_jitter(width = diff(range(d$age)) / 120, height = 0, seed = 1)) +
    geom_smooth(aes(group = cohort), method = "lm", formula = y ~ x, se = FALSE, colour = shade(hue, 0.35), linewidth = lw(0.7)) +
    note_corner(sprintf("%s, %d cohorts", name, uniqueN(d$cohort))) +
    note_corner(lab_stats(fit$r, fit$p), vjust = 2.6, parse = TRUE) +
    scale_y_continuous(labels = MINUS_LABELS, breaks = seq(-4, 4, 2)) +
    coord_cartesian(ylim = limits) +
    labs(x = xlab, y = ylab) +
    theme_pub()
}
a <- index_panel("mouse", MOUSE, "Mouse discovery", "Age (months)", c(-2.6, 4.6), "m6A index (SD)") +
  index_panel("external", EXTERNAL, "Mouse external", "Age (months)", c(-3.3, 5)) +
  index_panel("human", HUMAN, "Human", "Age (years)", c(-4.2, 8.4))

# b: age effect of the readers in mouse brain: combined (bar, with its q) and every cohort (dots)
b <- read_panel("fig4b")
b[, regulator := factor(regulator, c("FMR1", "FXR1", "FXR2", "YTHDF3"))]
total <- b[kind == "combined"]
b <- ggplot(b[kind == "cohort"], aes(regulator, effect)) +
  geom_col(data = total, fill = tint(MOUSE, 0.6), width = 0.68) +
  geom_hline(yintercept = 0, linewidth = lw(0.5)) +
  ggbeeswarm::geom_quasirandom(fill = MOUSE, shape = 21, colour = EDGE, stroke = st(0.3), size = 1.1, width = 0.2) +
  geom_text(data = total, aes(y = -31, label = paste0("q = ", signif(q, 2))), angle = 90, hjust = 0, size = pt_size(6)) +
  scale_y_continuous(labels = MINUS_LABELS, breaks = seq(-30, 10, 10)) +
  labs(x = NULL, y = "Age effect, mouse protein (%)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"))

# c: coefficient of the index on the reader level, readers among the control RNA-binding proteins, ranked
readers <- read_panel("fig4c")
controls <- read_panel("fig4c_controls")
rank_panel <- function(set, hue, xlab, at, top, ylab = NULL, slant = 0) {
  r <- readers[species == set, .(name = regulator, est, control_pct, reader = TRUE)]
  d <- rbind(controls[species == set, .(name = regulator, est, control_pct = NA_real_, reader = FALSE)], r)[order(est)]
  d[, rank := .I]
  shown <- d[reader == TRUE]
  shown[, label := paste0(name, ", ", ordinal(control_pct))]
  ggplot(d, aes(rank, est)) +
    ref_line(yintercept = 0) +
    geom_point(data = d[reader == FALSE], colour = GREY, size = 0.8, stroke = 0) +
    geom_point(data = shown, fill = hue, shape = 21, colour = EDGE, stroke = st(0.4), size = 1.5) +
    label_column(shown, "rank", "est", "label", at, top, 0.11, hjust = 1, slant = slant) +
    scale_x_continuous(breaks = c(0, 25, 50, 75), limits = c(-28, nrow(d) + 2), expand = c(0, 0)) +
    scale_y_continuous(labels = MINUS_LABELS, limits = c(-0.8, 0.8)) +
    labs(x = xlab, y = ylab) +
    theme_pub()
}
c_ <- rank_panel("mouse", MOUSE, "Rank, mouse (88 controls)", 66, 0.74, "Index per SD of the reader (SD)") +
  rank_panel("human", HUMAN, "Rank, human (90 controls)", 64, 0.74, slant = -3)

# d: the same coefficient in every unit (fill; size, evidence; open, P >= 0.05)
d <- read_panel("fig4d")
d[, `:=`(species = factor(SPECIES_NAME[species], SPECIES_NAME),
         regulator = factor(regulator, c("FMR1", "FMRP family", "YTHDF module", "YTHDF3")),
         evidence = pmin(-log10(p), 4))]
d[, unit := factor(unit, rev(unique(unit[order(species, unit)])))]
d <- ggplot(d, aes(regulator, unit)) +
  geom_point(data = ~ .x[p >= 0.05], shape = 1, size = 0.9, colour = GREY, stroke = st(0.5)) +
  geom_point(data = ~ .x[p < 0.05], aes(size = evidence, fill = est), shape = 21, colour = EDGE, stroke = st(0.3)) +
  ggh4x::facet_grid2(species ~ ., scales = "free_y", space = "free_y") +
  scale_fill_signed(1.5, "SD per SD") +
  scale_size_area(max_size = 2.8, name = expression(-log[10] ~ italic(P)), breaks = c(2, 4), labels = c("2", "≥4")) +
  scale_x_discrete(position = "top") +
  theme_tile() +
  theme(panel.grid.major = element_line(linewidth = lw(0.25), colour = "grey93"), legend.position = "bottom",
        legend.box = "vertical", legend.spacing.y = unit(0.5, "mm"), axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5),
        strip.text.y = element_text(size = 6), panel.spacing = unit(1.2, "mm"))

# e: FMR1 coupling with its m6A-marked targets against its unmarked targets in every proteome unit, two CLIP lists;
# filled, |marked coupling| above the 95th percentile of the control proteins of the unit
e <- read_panel("fig4e")[layer == "protein" & !is.na(b_u)]
e[, `:=`(species = factor(SPECIES_NAME[species], SPECIES_NAME), above = abs(b_s) > control_p95)]
counts <- e[, .(text = sprintf("%d of %d %s", sum(b_s > b_u), .N, tolower(species[1]))), by = .(list, species)]
counts[, y := fifelse(species == "Mouse", 0.455, 0.395)]
e <- ggplot(e, aes(b_u, b_s)) +
  ref_line(slope = 1) +
  ref_line(yintercept = 0) +
  ref_line(xintercept = 0) +
  geom_point(aes(colour = species, fill = interaction(species, above)), shape = 21, stroke = st(0.5), size = 1.2) +
  geom_text(data = counts, aes(x = -0.265, y = y, label = text), hjust = 0, size = pt_size(6)) +
  facet_wrap(~list) +
  scale_colour_manual(values = SPECIES, guide = "none") +
  scale_fill_manual(values = c("Mouse.TRUE" = MOUSE, "Human.TRUE" = HUMAN, "Mouse.FALSE" = "white", "Human.FALSE" = "white"), guide = "none") +
  scale_x_continuous(labels = MINUS_LABELS, breaks = c(-0.2, 0, 0.2)) +
  scale_y_continuous(labels = MINUS_LABELS, breaks = c(-0.2, 0, 0.2)) +
  coord_equal(xlim = c(-0.28, 0.5), ylim = c(-0.28, 0.5)) +
  labs(x = "Coupling with unmarked targets", y = "Coupling with m6A-marked targets") +
  theme_box()

# f: percentile of the YTHDF module and YTHDF3 coupling with the load among the control proteins; bold, at least 95th
f <- read_panel("fig4f")
f[, `:=`(row = factor(paste(SPECIES_NAME[species], layer), rev(levels(layer_rows(species, layer)))),
         regulator = factor(c(YTHDF = "YTHDF\nmodule", YTHDF3 = "YTHDF3")[regulator], c("YTHDF\nmodule", "YTHDF3")))]
f <- ggplot(f, aes(regulator, row, fill = pctile_signed)) +
  geom_tile(colour = "white", linewidth = lw(0.6)) +
  geom_text(aes(label = round(pctile_signed), fontface = fifelse(pctile_signed >= 95, "bold", "plain")), size = pt_size(6)) +
  scale_fill_gradientn(colours = c("#FFFFFF", tint(M6BRAIN, 0.75), tint(M6BRAIN, 0.45)), limits = c(0, 100), guide = "none") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_tile()

fig <- fig_panels(a, b, c_, d, e, f, design = paste(strrep("A", 40), paste0(strrep("B", 7), strrep("C", 18), strrep("D", 15)),
                                               paste0(strrep("E", 14), strrep("F", 11), strrep("D", 15)), sep = "\n"), heights = c(48, 60, 54))
save_fig(fig, "fig4", 168)
