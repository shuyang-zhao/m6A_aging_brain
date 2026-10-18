# Figure 3: m6A load marks proteins that decline with age
# reads source_data/fig3a, fig3b, fig3c_ladder, fig3c_checks, fig3c_null and fig3c_null_observed;
# writes output/figures/fig3.png and .pdf
source("figures/style.R")

PANELS <- c("Mouse discovery", "Mouse external", "Human protein")
COL <- c("Mouse discovery" = FOCAL, "Mouse external" = CONTRAST, "Human protein" = RED)
PER <- "%; human per decade"

# a: mean age effect of the genes by number of m6A sites (peaks in human): mean over cohorts, 95% t-interval across cohorts
a <- read_panel("fig3a")
a[, `:=`(panel = factor(panel, PANELS), sites = factor(sites, c("0", "1", "2-3", "≥4"), c("0", "1", "2–3", "≥4")))]
strips <- a[, .(n = n_cohorts[1]), by = panel][, setNames(sprintf("%s (%d cohorts%s)", panel, n,
                                                              fifelse(panel == "Human protein", ", per decade", "")), panel)]
pa <- ggplot(a, aes(sites, effect, group = panel)) +
  ZERO_Y +
  geom_linerange(aes(ymin = ci_low, ymax = ci_high), colour = CI, linewidth = 0.6) +
  geom_line(aes(colour = panel), linewidth = 1.1) +
  geom_point(aes(fill = panel), shape = 21, colour = EDGE, size = 3.8, stroke = 0.5) +
  facet_wrap(~panel, nrow = 1, scales = "free_y", labeller = as_labeller(strips)) +
  scale_colour_manual(values = COL, guide = "none") +
  scale_fill_manual(values = COL, guide = "none") +
  scale_y_continuous(labels = MINUS_LABELS, expand = expansion(mult = c(0.08, 0.1))) +
  labs(x = "m6A sites per gene (human: brain peaks)", y = "Mean age effect (%)") +
  theme_fig() +
  BOXED +
  theme(panel.spacing = unit(7, "mm"))

# b: load coefficient of every cohort (95% CI) and the random-effects estimate (diamond; I2 = heterogeneity); cohorts
# ordered by effect and aligned at the random-effects row, filled where the 95% CI excludes zero; a CI beyond the axis is cut and ends in an arrow
b <- read_panel("fig3b")
XLIM <- list("Mouse discovery" = c(-3.6, 0.4), "Mouse external" = c(-3.8, 0.8), "Human protein" = c(-1.6, 0.9))
M <- b[kind == "cohort", .N, by = panel][, max(N)] # rows of the longest forest
forest <- function(pn) {
  lim <- XLIM[[pn]]
  d <- b[panel == pn]
  co <- d[kind == "cohort"][order(effect)][, y := M - .N + seq_len(.N)]
  pooled <- d[kind == "pooled"][, y := M + 1.9]
  co[, sig := ci_low > 0 | ci_high < 0]
  note <- sprintf("%s%% (%s to %s)\n%s, I² = %.0f%%", signed(pooled$effect), MINUS(sprintf("%.2f", pooled$ci_low)),
                  MINUS(sprintf("%.2f", pooled$ci_high)), fmt_p(pooled$p), pooled$i2)
  cut <- co[ci_low < lim[1]]
  ggplot(co, aes(y = y)) +
    ZERO_X +
    geom_segment(aes(x = pmax(ci_low, lim[1]), xend = pmin(ci_high, lim[2]), yend = y), colour = CI, linewidth = 0.6) +
    geom_segment(data = cut, aes(x = lim[1] + 0.06 * diff(lim), xend = lim[1], yend = y), colour = CI, linewidth = 0.6,
                 arrow = arrow(length = unit(1.8, "mm"), type = "closed")) +
    geom_point(data = co[effect >= lim[1] & effect <= lim[2]],
               aes(effect, fill = I(fifelse(sig, COL[[pn]], "white")), colour = I(fifelse(sig, EDGE, GREY))),
               shape = 21, size = 3.6, stroke = 0.5) +
    geom_segment(data = pooled, aes(x = ci_low, xend = ci_high, yend = y), colour = CI, linewidth = 0.6) +
    geom_point(data = pooled, aes(effect), shape = 23, size = 5.4, fill = COL[[pn]], colour = EDGE, stroke = 0.5) +
    annotate("label", x = mean(lim), y = pooled$y + 0.9, label = note, size = PT, lineheight = 0.95, vjust = 1,
             fill = "white", linewidth = 0, border.colour = NA, label.padding = unit(0.6, "mm")) +
    scale_y_reverse(breaks = c(co$y, pooled$y), labels = c(co$cohort, "Random effects"), limits = c(M + 5, 0.3),
                    expand = expansion(0)) +
    scale_x_continuous(labels = MINUS_LABELS, expand = expansion(0)) +
    coord_cartesian(xlim = lim, clip = "off") +
    labs(x = paste0("Age effect per log unit\nof m6A load (", if (pn == "Human protein") "% per decade" else "%", ")"),
         y = NULL) +
    theme_rows() +
    theme(plot.margin = margin(4.4, 6, 3, 2, "mm"), axis.title.x = element_text(lineheight = 0.95))
}
pb <- wrap_plots(lapply(PANELS, forest), nrow = 1) + FLUSH

# c left: load coefficient with the covariates added in turn (ladder; the label is the share of the coefficient without
# covariates), with DRACH motifs as a covariate and with the ubiquitously expressed genes only; Mouse discovery and
# Human protein; 95% CI, random effects over cohorts
CP <- c("Mouse discovery", "Human protein")
STEP <- c(lengths = "+ Lengths", expression = "+ Expression", "neuron share" = "+ Neuron share",
          programmes = "+ Programmes (all covariates)")
CHECK <- c("DRACH as covariate" = "DRACH motifs as covariate", "ubiquitous genes" = "Ubiquitous genes only")
ladder <- read_panel("fig3c_ladder")[, row := STEP[step]][, y := match(step, names(STEP))]
checks <- read_panel("fig3c_checks")[check %in% names(CHECK)][, row := CHECK[check]][, y := 4.9 + match(check, names(CHECK))]
d <- rbind(ladder, checks, fill = TRUE)
d[, panel := factor(panel, CP)]
ticks <- unique(d[, .(y, row)])
cl <- ggplot(d, aes(effect, y)) +
  ZERO_X +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = CI, linewidth = 0.6, orientation = "y") +
  geom_point(aes(fill = panel), shape = 21, colour = EDGE, size = 3.8, stroke = 0.5) +
  geom_text(data = d[!is.na(retained_pct)], aes(x = Inf, label = paste0(round(retained_pct), "%")), hjust = 1.1,
            size = PT) +
  annotate("text", x = Inf, y = 0.1, label = "Retained", hjust = 1.1, size = PT) +
  facet_wrap(~panel, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = COL, guide = "none") +
  scale_y_reverse(breaks = ticks$y, labels = ticks$row, limits = c(7.4, -0.3), expand = expansion(0)) +
  scale_x_continuous(labels = MINUS_LABELS, expand = expansion(mult = c(0.05, 0.38)),
                     breaks = \(l) Filter(\(x) x <= 0, scales::breaks_pretty(4)(l))) +
  labs(x = paste0("Age effect per log unit of m6A load (", PER, ")"), y = NULL) +
  theme_rows() +
  BOXED +
  theme(panel.spacing = unit(5, "mm"))

# c right: the coefficient (line) against 1000 age-permuted copies of the load models (histogram), with the permutation P
nul <- read_panel("fig3c_null")[panel %in% CP][, panel := factor(panel, CP)]
obs <- read_panel("fig3c_null_observed")[, panel := factor(panel, CP)]
obs[, label := sprintf("%s%%\n%s", signed(effect), fmt_p(p_permutation))]
cn <- ggplot(nul, aes(effect)) +
  geom_histogram(bins = 36, fill = GREY, colour = "white", linewidth = 0.3) +
  geom_vline(data = obs, aes(xintercept = effect, colour = panel), linewidth = 1.3) +
  geom_text(data = obs, aes(x = Inf, y = Inf, label = label), hjust = 1.05, vjust = 1.2, size = PT, lineheight = 0.95) +
  facet_wrap(~panel, nrow = 1, scales = "free") +
  scale_colour_manual(values = COL, guide = "none") +
  scale_x_continuous(labels = MINUS_LABELS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.5))) +
  labs(x = paste0("Coefficient with permuted ages (", PER, ")"), y = "Permuted copies") +
  theme_fig() +
  BOXED +
  theme(panel.spacing = unit(5, "mm"))

pc <- (cl | cn) + plot_layout(widths = c(200, 166)) + FLUSH
fig <- panel_grid(list(a = pa, b = pb, c = pc), design = "A\nB\nC", heights = c(84, 112, 90))
save_fig(fig, "fig3", 286)
