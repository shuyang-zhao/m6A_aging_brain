# Figure 5: the m6Brain clock
# reads source_data/fig5a-e and fig5d_training (script 55); writes output/figures/fig5.png and .pdf
source("figures/style.R")

SETS <- c("held-out mouse" = FOCAL, "external mouse" = CONTRAST, "held-out human" = RED)
SET_LABELS <- c("held-out mouse" = "Held-out mouse", "external mouse" = "External mouse",
                "held-out human" = "Held-out human")
set_of <- function(set, species) factor(paste(set, species), names(SETS))

# a: Pearson r of the predicted with chronological age in every cohort (point, Fisher-z 95% CI; filled, P < 0.05),
# cohorts ordered by r within the set, and the random-effects estimate over the cohorts of the set (diamond)
a <- read_panel("fig5a")
a[, `:=`(group = set_of(set, species), pooled = kind == "pooled")]
setorder(a, group, pooled, -r)
a[, x := factor(fifelse(pooled, paste0("Pooled ", group), cohort), unique(fifelse(pooled, paste0("Pooled ", group),
                                                                               cohort)))]
a[, `:=`(fill = fifelse(pooled | p < 0.05, SETS[as.character(group)], "white"))]
pooled_text <- a[pooled == TRUE, .(group, x, text = sprintf("r = %.2f", r))]
pa <- ggplot(a, aes(x, r)) +
  geom_linerange(aes(ymin = ci_low, ymax = ci_high), colour = CI, linewidth = 0.6) +
  geom_point(data = ~ .x[pooled == FALSE], aes(fill = I(fill)), colour = EDGE, shape = 21, size = 3.4, stroke = 0.5) +
  geom_point(data = ~ .x[pooled == TRUE], aes(fill = I(fill)), colour = EDGE, shape = 23, size = 5, stroke = 0.5) +
  geom_text(data = pooled_text, aes(x, 0.36, label = text), hjust = 0, angle = 90, size = PT, colour = INK) +
  facet_grid(. ~ group, scales = "free_x", space = "free_x", labeller = labeller(group = SET_LABELS)) +
  scale_x_discrete(labels = function(x) sub("^Pooled .*", "Pooled", x)) +
  scale_y_continuous(breaks = seq(0.4, 1, 0.2), limits = c(0.34, 1)) +
  labs(x = NULL, y = "Pearson r of predicted\nwith chronological age") +
  theme_fig() +
  BOXED +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.spacing.x = unit(3, "mm"))

# b: predicted against chronological age (donor means of the clock's raw prediction) in four cohorts; dashed line,
# identity; r of the cohort
b <- read_panel("fig5b")
b[, `:=`(key = paste(set, species), cohort = factor(cohort, unique(cohort)))]
cohort_scatter <- function(d, y_title) {
  unit <- d$age_unit[1]
  d[, strip := sprintf("%s, %s", cohort, tolower(SET_LABELS[key]))]
  note <- unique(d[, .(strip, text = sprintf("r = %.2f\nn = %d", r, n))])
  lim <- range(c(d$age, d$predicted_age))
  ggplot(d, aes(age, predicted_age)) +
    geom_abline(slope = 1, intercept = 0, colour = EDGE, linewidth = 0.6, linetype = "22") +
    geom_point(fill = SETS[d$key[1]], colour = EDGE, shape = 21, size = 2.6, stroke = 0.4, alpha = 0.9,
               position = position_jitter(width = 0.02 * diff(lim), height = 0, seed = 1)) +
    geom_text(data = note, aes(x = -Inf, y = Inf, label = text), hjust = -0.08, vjust = 1.15, size = PT,
              lineheight = 0.95, colour = INK, inherit.aes = FALSE) +
    facet_wrap(~strip) +
    scale_x_continuous(breaks = scales::breaks_pretty(4)) +
    scale_y_continuous(breaks = scales::breaks_pretty(4), expand = expansion(mult = c(0.05, 0.3))) +
    coord_cartesian(xlim = lim) +
    labs(x = paste0("Age (", unit, ")"), y = y_title) +
    theme_fig() +
    BOXED
}
scatters <- Map(cohort_scatter, split(b, b$cohort), y_title = list("Predicted age", NULL, "Predicted age", NULL))
pb <- wrap_plots(scatters, nrow = 2) + FLUSH

# c: benchmark, pooled r (95% CI) of every clock per set; rows grouped (the m6Brain clock and the two reference nets,
# then the published clocks by kind) and ordered by the mean r over the sets; dashed line, the m6Brain clock
KINDS <- c("m6Brain clock and references", "Proteomic clocks", "Transcriptome clocks", "Signatures and half-lives")
KIND_LABELS <- c("m6Brain\nand references", "Proteomic\nclocks", "Transcriptome\nclocks", "Signatures,\nhalf-lives")
MINE <- c(m6brain_clock = FOCAL, elastic_net = GREY, m6a_load_prior = FOCAL_LIGHT)
cc <- read_panel("fig5c")
cc[, `:=`(group = set_of(set, species), kind = factor(kind, KINDS, KIND_LABELS))]
cc[, mean_r := mean(r), by = weights]
row_order <- unique(cc[order(kind, -(weights == "m6brain_clock"), -mean_r), label])
cc[, `:=`(row = factor(label, rev(row_order)), fill = fifelse(weights %in% names(MINE), MINE[weights], "white"))]
mine <- cc[weights == "m6brain_clock", .(group, r)]
pc <- ggplot(cc, aes(r, row)) +
  geom_vline(data = mine, aes(xintercept = r), colour = FOCAL, linewidth = 0.5, linetype = "22") +
  ZERO_X +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = CI, linewidth = 0.5) +
  geom_point(aes(fill = I(fill), shape = weights == "m6brain_clock", size = weights == "m6brain_clock"),
             colour = EDGE, stroke = 0.5) +
  facet_grid(kind ~ group, scales = "free_y", space = "free_y", labeller = labeller(group = SET_LABELS)) +
  scale_shape_manual(values = c("TRUE" = 23, "FALSE" = 21), guide = "none") +
  scale_size_manual(values = c("TRUE" = 4.4, "FALSE" = 3), guide = "none") +
  scale_x_continuous(labels = MINUS_LABELS, breaks = c(-0.5, 0, 0.5, 1), limits = c(-0.75, 1)) +
  labs(x = "Pooled Pearson r with age (95% CI)", y = NULL) +
  theme_rows() +
  BOXED +
  theme(panel.spacing.y = unit(2.4, "mm"), panel.spacing.x = unit(3, "mm"),
        strip.text.y = element_text(angle = 0, hjust = 0, margin = margin(l = 2, unit = "mm")))

# d: within-cohort error, MAE of 10 x 10 nested cross-validation (mean ± SD over the seeds) of the three nets and of the
# training-mean null, with the paired Wilcoxon P against the unguided net; right, MAE by number of training donors in
# BrainGVEX (mean ± SD over 5 seeds) with the predictions intact, permuted across genes or not used (unguided)
VARIANT_COL <- c(m6brain_clock = FOCAL, elastic_net = GREY, m6a_load_prior = FOCAL_LIGHT, null = "white",
                 shuffled_prior = RED_LIGHT)
w <- read_panel("fig5d")
w[, `:=`(strip = sprintf("%s (%s)", cohort, age_unit), row = factor(variant, rev(names(VARIANT_COL)[1:4]),
                                                                     rev(unique(label[order(match(variant,
                                                                     names(VARIANT_COL)))]))))]
w[, strip := factor(strip, unique(strip))]
w[, note := fifelse(is.na(p_wilcoxon), "", fmt_p(p_wilcoxon))]
pd1 <- ggplot(w, aes(mae_mean, row)) +
  geom_col(aes(fill = I(VARIANT_COL[variant])), colour = EDGE, linewidth = 0.5, width = 0.66) +
  geom_errorbar(aes(xmin = mae_mean - mae_sd, xmax = mae_mean + mae_sd), width = 0.3, colour = EDGE,
                linewidth = 0.5, orientation = "y") +
  geom_text(aes(x = mae_mean + mae_sd, label = note), hjust = -0.12, size = PT, colour = INK) +
  facet_wrap(~strip, ncol = 1, scales = "free_x") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.75)), breaks = scales::breaks_pretty(3)) +
  labs(x = "MAE", y = NULL) +
  theme_rows() +
  BOXED
t <- read_panel("fig5d_training")
t[, `:=`(x = factor(n_train, c("25", "50", "100", "all"), c("25", "50", "100", "All")),
         variant = factor(variant, c("m6brain_clock", "elastic_net", "shuffled_prior")))]
dodge <- position_dodge(0.4)
pd2 <- ggplot(t, aes(x, mae_mean, group = variant)) +
  geom_line(aes(colour = variant), linewidth = 0.8, position = dodge) +
  geom_errorbar(aes(ymin = mae_mean - mae_sd, ymax = mae_mean + mae_sd), width = 0, colour = CI, linewidth = 0.5,
                position = dodge) +
  geom_point(aes(fill = variant), shape = 21, colour = EDGE, size = 3, stroke = 0.5, position = dodge) +
  scale_colour_manual(values = c(m6brain_clock = FOCAL, elastic_net = "#8C8C8C", shuffled_prior = RED),
                      guide = "none") +
  scale_fill_manual(values = VARIANT_COL[levels(t$variant)], labels = c("m6Brain clock", "Unguided",
                                                                         "Permuted"), name = NULL) +
  facet_wrap(~"Luo 2024 (years)") +
  scale_y_continuous(limits = c(3.4, NA), breaks = 5:10) +
  labs(x = "Training donors", y = "MAE") +
  theme_fig() +
  BOXED +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.02), legend.justification = c(0, 0))
pd <- (pd1 + pd2 + plot_layout(widths = c(1, 1))) + FLUSH

# e: the same donors, r of the m6Brain clock on the protein and of the transcriptome clocks on the RNA (95% CI for the
# m6Brain clock and the best transcriptome clock, with the P of their difference; open circles, the other clocks)
e <- read_panel("fig5e")
e[, mine := weights == "m6brain_clock"]
e[, rank := frank(-r, ties.method = "first"), by = .(cohort, mine)]
e[, row := fifelse(mine, "m6Brain clock, protein", fifelse(rank == 1, paste0(label, ", RNA"), "Other RNA clocks"))]
e[, strip := sprintf("%s, n = %d", cohort, n)]
e[, row := factor(row, rev(unique(row[order(-mine, rank)])))]
shown <- e[mine | rank == 1]
pe <- ggplot(e, aes(r, row)) +
  ZERO_X +
  geom_point(data = e[!mine & rank > 1], shape = 21, fill = "white", colour = "#8C8C8C", size = 3, stroke = 0.5) +
  geom_linerange(data = shown, aes(xmin = ci_low, xmax = ci_high), colour = CI, linewidth = 0.6) +
  geom_point(data = shown, aes(fill = mine, shape = mine, size = mine), colour = EDGE, stroke = 0.5) +
  geom_text(data = shown[!is.na(p_difference)], aes(x = Inf, label = fmt_p(p_difference)), hjust = -0.1, size = PT) +
  facet_wrap(~strip, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("TRUE" = FOCAL, "FALSE" = GREY), guide = "none") +
  scale_shape_manual(values = c("TRUE" = 23, "FALSE" = 21), guide = "none") +
  scale_size_manual(values = c("TRUE" = 4.4, "FALSE" = 3.4), guide = "none") +
  scale_x_continuous(labels = MINUS_LABELS, breaks = c(0, 0.5, 1), limits = c(-0.2, 1)) +
  coord_cartesian(clip = "off") +
  labs(x = "Pearson r with age on the same donors", y = NULL) +
  theme_rows() +
  BOXED +
  theme(plot.margin = margin(4.4, 24, 3, 2, "mm"))

fig <- panel_grid(list(a = pa, b = pb, c = pc, d = pd, e = pe), design = "ABB\nCCD\nCCE",
                  widths = c(150, 70, 146), heights = c(120, 105, 95))
save_fig(fig, "fig5", 320)
