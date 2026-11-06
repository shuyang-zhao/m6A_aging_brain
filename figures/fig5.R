# Figure 5: the m6Brain clock
# reads source_data/fig5a-g (scripts 52 and 55); writes output/figures/fig5
source("figures/style.R")

SETS <- c("held-out mouse" = "Held-out mouse", "external mouse" = "External mouse", "held-out human" = "Held-out human")
set_of <- function(set, species) factor(SETS[paste(set, species)], SETS)
METHOD <- c("m6Brain clock" = M6BRAIN, "Unguided elastic net" = GREY, "m6A-load-guided net" = LOAD)

# a: proteins chosen by each net of all cohorts by decile of the predicted |age z| that guides the clock
a <- read_panel("fig5a")
a[, `:=`(species = factor(SPECIES_NAME[species], SPECIES_NAME), label = factor(label, names(METHOD)[1:2]))]
key <- data.table(species = factor("Mouse", levels(a$species)), y = c(150, 130), text = c("m6Brain clock", "Unguided net"),
                  label = factor(names(METHOD)[1:2], names(METHOD)[1:2]))
a <- ggplot(a, aes(decile, chosen, fill = label)) +
  geom_col(position = position_dodge(width = 0.82), width = 0.78) +
  geom_point(data = key, aes(1, y, colour = label), shape = 15, size = 1.8, inherit.aes = FALSE) +
  geom_text(data = key, aes(1.7, y, label = text), hjust = 0, size = pt_size(6), inherit.aes = FALSE) +
  facet_wrap(~species) +
  scale_fill_manual(values = METHOD, guide = "none") +
  scale_colour_manual(values = METHOD, guide = "none") +
  scale_x_continuous(breaks = c(1, 5, 10)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Decile of the predicted |age z|", y = "Proteins in the clock") +
  theme_box()

# b: predicted against chronological age in four cohorts with many ages (donor means); dashed, identity
b <- read_panel("fig5b")
b[, group := set_of(set, species)]
b[, strip := factor(cohort, unique(cohort))]
labels <- b[, .(text = lab_stats(r[1], n = n[1])), by = strip]
b <- ggplot(b, aes(age, predicted_age)) +
  geom_abline(slope = 1, linetype = "22", linewidth = lw(0.5), colour = "grey55") +
  geom_point(aes(fill = group), shape = 21, colour = EDGE, stroke = st(0.25), size = 1, alpha = 0.9) +
  geom_text(data = labels, aes(x = -Inf, y = Inf, label = text), parse = TRUE, hjust = -0.06, vjust = 1.3, size = pt_size(6)) +
  facet_wrap(~strip, nrow = 1, scales = "free") +
  scale_fill_manual(values = SET, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(x = "Chronological age (mouse, months; human, years)", y = "Predicted age") +
  theme_box()

# c: Pearson r of the m6Brain clock against the best published clock of the same cohort; dashed, identity
c_ <- read_panel("fig5c")[kind == "cohort"]
c_[, group := set_of(set, species)]
above <- c_[, sum(r > r_best)]
c_ <- ggplot(c_, aes(r_best, r)) +
  ref_line(slope = 1) +
  geom_point(aes(fill = group), shape = 21, colour = EDGE, stroke = st(0.3), size = 1.5) +
  note(0.03, 1, sprintf("%d of %d cohorts\nabove identity", above, nrow(c_)), hjust = 0, vjust = 1, lineheight = 0.95) +
  key_rows(names(SET), SET, x = 0.3, y = 0.24, step = 0.09, width = 0.05, square = TRUE) +
  scale_fill_manual(values = SET, guide = "none") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Best published clock, r", y = "m6Brain clock, r") +
  theme_pub()

# d: pooled Pearson r per set of the m6Brain clock, its reference nets and the 27 published clocks and signatures
d <- read_panel("fig5d")
KINDS <- c("m6Brain clock and references" = "m6Brain", "Proteomic clocks" = "Proteomic clocks",
           "Transcriptome clocks" = "Transcriptome clocks", "Signatures and half-lives" = "Signatures")
d[, `:=`(group = factor(sub(" ", "\n", SETS[paste(set, species)]), sub(" ", "\n", SETS)), block = factor(KINDS[kind], KINDS))]
d[, mean_r := mean(r), by = label]
rows <- unique(d[order(block, -(label == "m6Brain clock"), -mean_r), label])
d[, row := factor(label, rev(rows))]
d <- ggplot(d, aes(group, row, fill = r)) +
  geom_tile(colour = "white", linewidth = lw(0.5)) +
  geom_text(aes(label = number(r), fontface = fifelse(label == "m6Brain clock", "bold", "plain")), size = pt_size(6)) +
  ggh4x::facet_grid2(block ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_signed(1, "Pooled r with age", light = TRUE) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_tile() +
  theme(legend.position = "bottom", legend.margin = margin(t = 1, unit = "mm"), panel.spacing = unit(1.2, "mm"),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 90, size = 6), strip.clip = "off",
        strip.background = element_rect(fill = LIGHT, colour = LIGHT))

# e: mean absolute error of 10 x 10 nested cross-validation within two cohorts (bar, mean; dots, seeds), with the
# paired Wilcoxon P against the unguided net
e <- read_panel("fig5e")
seeds <- read_panel("fig5e_seeds")
ROWS <- c("m6Brain clock" = "m6Brain\nclock", "Unguided elastic net" = "Unguided\nnet", "m6A-load-guided net" = "m6A load",
          "Training-mean null" = "Training\nmean")
FILL <- setNames(c(M6BRAIN, GREY, LOAD, "white"), ROWS)
strips <- e[, sprintf("%s (%s)", cohort, age_unit)]
e[, `:=`(y = factor(ROWS[label], rev(ROWS)), strip = factor(strips, unique(strips)))]
seeds[, `:=`(y = factor(ROWS[label], rev(ROWS)), strip = factor(e$strip[match(cohort, e$cohort)], levels(e$strip)))]
tested <- e[label == "m6Brain clock", .(strip, y, x = mae_mean + mae_sd * 2, text = fmt_p(p_wilcoxon, digits = 1))]
e <- ggplot(e, aes(mae_mean, y)) +
  geom_col(aes(fill = y), colour = EDGE, linewidth = lw(0.4), width = 0.72, orientation = "y") +
  ggbeeswarm::geom_quasirandom(data = seeds, aes(x = mae), size = 0.55, stroke = 0, width = 0.2, colour = "black", orientation = "y") +
  geom_text(data = tested, aes(x, y, label = text), parse = TRUE, hjust = -0.08, size = pt_size(6), inherit.aes = FALSE) +
  facet_wrap(~strip, scales = "free_x") +
  scale_fill_manual(values = FILL, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.62)), breaks = scales::breaks_pretty(3)) +
  labs(x = "MAE", y = NULL) +
  theme_box() +
  theme(axis.ticks.y = element_blank())

# f: error by number of training donors in Luo 2024 (mean and SD over 5 seeds), with the predictions, with the
# predictions permuted across genes and without them
f <- read_panel("fig5f")
LINES <- c("Permuted predictions" = "#6E6E6E", "Unguided elastic net" = GREY, "m6Brain clock" = M6BRAIN)
f[, `:=`(x = match(n_train, c("25", "50", "100", "all")), label = factor(label, names(LINES)))]
f <- ggplot(f, aes(x, mae_mean, colour = label, fill = label)) +
  geom_ribbon(aes(ymin = mae_mean - mae_sd, ymax = mae_mean + mae_sd), colour = NA, alpha = 0.2) +
  geom_line(aes(linetype = label), linewidth = lw(0.9)) +
  geom_point(size = 1) +
  key_rows(c("Permuted", "Unguided", "m6Brain clock"), LINES, x = 2, y = 9.9, step = 0.55, width = 0.28) +
  scale_colour_manual(values = LINES, guide = "none") +
  scale_fill_manual(values = LINES, guide = "none") +
  scale_linetype_manual(values = c("22", "solid", "solid"), guide = "none") +
  scale_x_continuous(breaks = 1:4, labels = c("25", "50", "100", "All"), expand = expansion(add = 0.2)) +
  labs(x = "Training donors, Luo 2024", y = "MAE (years)") +
  theme_pub()

# g: the same donors: r of the m6Brain clock on the protein and of each transcriptome clock on the RNA, sorted; P of the
# difference from the best transcriptome clock
g <- read_panel("fig5g")
g[, `:=`(mine = weights == "m6brain_clock", strip = sprintf("%s, n = %d", cohort, n))]
g[, strip := factor(strip, rev(unique(strip)))]
setorder(g, strip, -r)
g[, x := seq_len(.N), by = strip]
best <- g[mine == FALSE, .SD[which.max(r)], by = strip]
tested <- best[, .(strip, text = fmt_p(p_difference))]
mine <- g[mine == TRUE]
g <- ggplot(g, aes(x, r)) +
  geom_hline(yintercept = 0, linewidth = lw(0.5)) +
  geom_col(aes(fill = mine), width = 0.76) +
  geom_text(data = best, aes(x = 3, y = -0.17, label = paste("Best RNA", number(r))), hjust = 0, size = pt_size(6)) +
  geom_text(data = mine, aes(x = 3, y = -0.32, label = paste("m6Brain", number(r))), hjust = 0, size = pt_size(6)) +
  geom_text(data = tested, aes(x = 3, y = -0.47, label = text), hjust = 0, parse = TRUE, size = pt_size(6), inherit.aes = FALSE) +
  facet_wrap(~strip) +
  scale_fill_manual(values = c("TRUE" = M6BRAIN, "FALSE" = GREY), guide = "none") +
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), expand = expansion(mult = 0.02)) +
  coord_cartesian(ylim = c(-0.55, 1)) +
  labs(x = "Clocks, ranked", y = "r with age, same donors") +
  theme_box() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.text = element_text(size = 6))

fig <- fig_panels(a, b, c_, d, e, f, g, design = "AAAAAAABBBBBBBBBBB\nCCCCCDDDDDDDEEEEEE\nFFFFFDDDDDDDGGGGGG",
                  heights = c(46, 62, 56))
save_fig(fig, "fig5", 170)
