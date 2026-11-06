# Figure 3: m6A load marks proteins that decline with age
# reads source_data/fig3a-e (scripts 21 and 50); writes output/figures/fig3
source("figures/style.R")

SETS <- c("Mouse discovery" = MOUSE, "Mouse external" = EXTERNAL, "Human protein" = HUMAN)

# a: distribution of the age z of every protein by number of m6A sites (peaks in human)
a <- read_panel("fig3a")
a[, sites := factor(sites, SITES)]
ecdf_panel <- function(target_name, hue, limit, sites_label) {
  d <- a[target == target_name]
  n <- d[, .N, keyby = sites]
  step <- 0.085
  ggplot(d, aes(observed, colour = sites)) +
    ref_line(yintercept = 0.5) +
    ref_line(xintercept = 0) +
    stat_ecdf(geom = "step", n = 256, linewidth = lw(0.9)) +
    key_rows(paste(n$sites, sprintf("(n = %s)", thousands(n$N))), sites_ramp(hue),
             x = limit * 0.08, y = 0.33, step = step, width = limit * 0.12) +
    note(limit * 0.08, 0.33 + step, sites_label, hjust = 0) +
    scale_colour_manual(values = sites_ramp(hue), guide = "none") +
    scale_x_continuous(labels = MINUS_LABELS) +
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), expand = expansion(mult = 0.01)) +
    coord_cartesian(xlim = c(-limit, limit)) +
    labs(x = paste("Age z,", tolower(target_name)), y = "Cumulative fraction") +
    theme_pub()
}
a <- ecdf_panel("Mouse protein", MOUSE, 8, "m6A sites") +
  (ecdf_panel("Human protein", HUMAN, 5, "m6A peaks") + labs(y = NULL))

# b: per cohort, effect of 1, 2-3 and 4 or more sites against none (%; human % per decade) and load coefficient per SD
# of load; colour, z of each coefficient; last row of each set, the random-effects estimate
sites <- read_panel("fig3b_sites")
load <- read_panel("fig3b_load")
strips <- c("Mouse discovery" = "Mouse discovery", "Mouse external" = "Mouse external", "Human protein" = "Human")
b <- rbind(sites[, .(panel, cohort, kind, column = sites, effect, z, block = "Sites against none (%)")],
           load[, .(panel, cohort, kind, column = "Load", effect, z, block = "Per SD (%)")])
b[, `:=`(row = factor(paste(panel, cohort), rev(load[order(match(panel, names(SETS)), kind == "pooled", effect),
                                                     paste(panel, cohort)])),
         set = factor(strips[panel], strips), block = factor(block, c("Sites against none (%)", "Per SD (%)")),
         column = factor(column, c(SITES[-1], "Load")),
         text = mapply(number, effect, fifelse(panel == "Human protein", 2, 1)))]
b <- ggplot(b, aes(column, row, fill = z)) +
  geom_tile(colour = "white", linewidth = lw(0.6)) +
  geom_text(aes(label = text, fontface = fifelse(kind == "pooled", "bold", "plain")), size = pt_size(6)) +
  ggh4x::facet_grid2(
    set ~ block, scales = "free", space = "free", switch = "y",
    strip = ggh4x::strip_themed(background_y = ggh4x::elem_list_rect(fill = tint(unname(SETS), 0.45),
                                                                     colour = tint(unname(SETS), 0.45)))
  ) +
  scale_fill_signed(6, "z", light = TRUE) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(labels = setNames(load$cohort, paste(load$panel, load$cohort)), expand = c(0, 0)) +
  theme_tile() +
  theme(legend.position = "bottom", legend.margin = margin(t = 1, unit = "mm"), strip.placement = "outside",
        strip.text.x = element_text(size = 6), strip.text.y.left = element_text(angle = 90, size = 6),
        panel.spacing.x = unit(1.2, "mm"), panel.spacing.y = unit(1.4, "mm"))

# c: the m6A load and gene length in the same model (95% CI)
c_ <- read_panel("fig3c")
c_[, `:=`(set = factor(panel, names(SETS)), predictor = factor(predictor, c("m6A load", "Gene length")))]
levels(c_$predictor) <- c("m6A load", "Length")
c_[, fill := fifelse(predictor == "Length", GREY, SETS[as.character(panel)])]
c_[, label := sprintf("'%s: '*%s", fifelse(predictor == "Length", "Length", "Load"), fmt_p(p))]
titles <- c("Mouse discovery" = "Mouse discovery", "Mouse external" = "Mouse external", "Human protein" = "Human, per decade")
c_ <- ggplot(c_, aes(predictor, effect)) +
  geom_hline(yintercept = 0, linewidth = lw(0.5)) +
  geom_col(aes(fill = I(fill)), width = 0.66) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0, linewidth = lw(0.5), colour = EDGE) +
  geom_text(data = ~ .x[predictor == "m6A load"], aes(x = -Inf, y = -Inf, label = label), parse = TRUE, hjust = -0.05,
            vjust = -1.75, size = pt_size(6)) +
  geom_text(data = ~ .x[predictor == "Length"], aes(x = -Inf, y = -Inf, label = label), parse = TRUE, hjust = -0.05,
            vjust = -0.4, size = pt_size(6)) +
  facet_wrap(~set, scales = "free_y", labeller = as_labeller(titles)) +
  scale_y_continuous(labels = MINUS_LABELS, expand = expansion(mult = c(0.34, 0.04))) +
  labs(x = NULL, y = "Coefficient (% per SD)") +
  theme_box()

# d: load coefficient retained as covariates are added (solid line) and in two further models, as a percentage of the
# coefficient without covariates (95% CI); the design of each model under the axis
ladder <- read_panel("fig3d_ladder")
checks <- read_panel("fig3d_checks")[check != "m6A load"]
base <- ladder[step == "none", .(panel, none = effect)]
STEPS <- c("none", "lengths", "expression", "neuron share", "programmes", "DRACH as covariate", "ubiquitous genes")
d <- rbind(ladder[, .(panel, model = step, effect, ci_low, ci_high)],
           checks[, .(panel, model = check, effect, ci_low, ci_high)])
d <- d[base, on = "panel"]
d[, `:=`(x = match(model, STEPS), mean = 100 * effect / none, low = 100 * ci_high / none, high = 100 * ci_low / none,
         species = sub(" .*", "", panel))]
line <- d[x <= 5]
plot_d <- ggplot(d, aes(x, mean, colour = species)) +
  ref_line(yintercept = 0) +
  ref_line(yintercept = 100) +
  geom_vline(xintercept = 5.5, linewidth = lw(0.5), colour = LIGHT) +
  geom_ribbon(data = line, aes(ymin = low, ymax = high, fill = species), colour = NA, alpha = 0.18) +
  geom_line(data = line, linewidth = lw(0.9)) +
  geom_errorbar(data = d[x > 5], aes(ymin = low, ymax = high), width = 0, linewidth = lw(0.5),
                position = position_dodge(width = 0.3)) +
  geom_point(data = line, aes(fill = species), shape = 21, colour = EDGE, stroke = st(0.4), size = 1.5) +
  geom_point(data = d[x > 5], aes(fill = species), shape = 21, colour = EDGE, stroke = st(0.4), size = 1.5,
             position = position_dodge(width = 0.3)) +
  key_rows(c("Mouse discovery", "Human"), SPECIES, x = 2.1, y = 186, step = 15, width = 0.42) +
  scale_colour_manual(values = SPECIES, guide = "none") +
  scale_fill_manual(values = SPECIES, guide = "none") +
  scale_x_continuous(limits = c(0.5, 7.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 50, 100, 150), labels = MINUS_LABELS) +
  coord_cartesian(ylim = c(0, 200)) +
  labs(x = NULL, y = "Coefficient retained (%)") +
  theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
design <- data.frame(
  Lengths = c(FALSE, rep(TRUE, 6)), Expression = c(FALSE, FALSE, rep(TRUE, 5)), "Neuron share" = c(rep(FALSE, 3), rep(TRUE, 4)),
  Programs = c(rep(FALSE, 4), rep(TRUE, 3)), "DRACH motifs" = c(rep(FALSE, 5), TRUE, FALSE),
  "Ubiquitous genes only" = c(rep(FALSE, 6), TRUE), check.names = FALSE
)
d <- plot_d / design_strip(design) + plot_layout(heights = c(1, 0.62))

# e: coefficient of the load after permuting ages within the blocks of each cohort (grey) and the observed coefficient
null <- read_panel("fig3e_null")
observed <- read_panel("fig3e_observed")
null_panel <- function(name, hue, label) {
  o <- observed[panel == name]
  ggplot(null[panel == name], aes(effect)) +
    geom_histogram(bins = 36, fill = GREY, colour = NA) +
    geom_hline(yintercept = 0, linewidth = lw(0.5)) +
    annotate("segment", x = o$effect, xend = o$effect, y = 128, yend = 3, colour = hue, linewidth = lw(1),
             arrow = grid::arrow(length = unit(1.1, "mm"), type = "closed")) +
    note(o$effect, 136, paste0("Observed ", number(o$effect), "%"), hjust = 0) +
    note(o$effect, 151, fmt_p(signif(o$p_permutation, 2)), hjust = 0, parse = TRUE) +
    scale_x_continuous(labels = MINUS_LABELS) +
    scale_y_continuous(breaks = c(0, 50, 100), expand = expansion(mult = c(0, 0.04))) +
    coord_cartesian(ylim = c(0, 160)) +
    labs(x = label, y = "Permuted copies") +
    theme_pub()
}
e <- null_panel("Mouse discovery", MOUSE, "Coefficient, mouse (%)") +
  (null_panel("Human protein", HUMAN, "Coefficient, human (%)") + labs(y = NULL))

fig <- fig_compose(a, b, c_, d, e, design = "AAAAABBBB\nCCCCCBBBB\nDDDDDEEEE", heights = c(58, 46, 58))
save_fig(fig, "fig3", 168)
