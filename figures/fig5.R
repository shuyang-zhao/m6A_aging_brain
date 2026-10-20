# Figure 5: the zero-shot proteome clock
# reads source_data/fig5a, fig5b, fig5c, fig5d and fig5d_null (script 52); writes output/figures/fig5.png and .pdf
source("figures/style.R")

SETS <- c("held-out mouse" = FOCAL, "external mouse" = CONTRAST, "held-out human" = RED)
SET_LABELS <- c("held-out mouse" = "Held-out mouse", "external mouse" = "External mouse",
                "held-out human" = "Held-out human")

# a: Spearman rho of the clock score with age in every cohort (point; filled where P < 0.05) with its Fisher-z 95%
# interval, and the random-effects estimate over the cohorts of the set (diamond; I2, the heterogeneity)
a <- read_panel("fig5a")
a[, `:=`(group = factor(paste(set, species), names(SETS)), pooled = kind == "pooled")]
a[, label := fifelse(pooled, "Pooled", sprintf("%s, n = %d", cohort, n))]
a[, label := factor(label, rev(unique(label)))]
a[, `:=`(fill = fifelse(pooled | p < 0.05, SETS[as.character(group)], "white"),
         edge = fifelse(pooled | p < 0.05, EDGE, GREY))]
pooled_text <- a[pooled == TRUE, .(group, label, x = -0.86,
                                   text = sprintf("ρ = %.2f, I² = %.0f%%", rho, i2))]
a <- ggplot(a, aes(rho, label)) +
  ZERO_X +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = CI, linewidth = 0.6) +
  geom_point(data = ~ .x[pooled == FALSE], aes(fill = I(fill), colour = I(edge)), shape = 21, size = 3.6, stroke = 0.5) +
  geom_point(data = ~ .x[pooled == TRUE], aes(fill = I(fill)), colour = EDGE, shape = 23, size = 5.2, stroke = 0.5) +
  geom_label(data = pooled_text, aes(x, label, label = text), hjust = 0, size = PT, colour = INK, fill = "white",
             border.colour = NA, label.padding = unit(0.6, "mm"), inherit.aes = FALSE) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", labeller = labeller(group = function(x) SET_LABELS[x])) +
  scale_x_continuous(labels = MINUS_LABELS, breaks = seq(-0.5, 1, 0.5), limits = c(-0.9, 1.02)) +
  labs(x = "Spearman ρ of the score with age", y = NULL) +
  theme_rows() +
  BOXED +
  theme(panel.spacing.y = unit(4, "mm"), strip.text.y = element_text(angle = -90))

# b: clock score against age (points; line, least squares) in four cohorts, as in the source data; Spearman rho and P
b <- read_panel("fig5b")
b[, `:=`(key = paste(set, species), cohort = factor(cohort, unique(cohort)))]
cohort_scatter <- function(d, y_title) {
  ages <- uniqueN(d$age)
  d[, strip := sprintf("%s, %s", cohort, tolower(SET_LABELS[key]))]
  # two-age cohorts: the text sits between the two groups of points, which is empty
  note <- unique(d[, .(strip, text = paste0("ρ = ", sprintf("%.2f", rho), "\n", fmt_p(p)),
                       x = if (ages == 2) mean(range(age)) else -Inf, hjust = if (ages == 2) 0.5 else -0.1)])
  ggplot(d, aes(age, score)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = EDGE, linewidth = 0.8) +
    geom_point(fill = SETS[d$key[1]], colour = EDGE, shape = 21, size = 2.6, stroke = 0.4, alpha = 0.9,
               position = position_jitter(width = 0.04 * diff(range(d$age)), height = 0, seed = 1)) +
    geom_text(data = note, aes(x = x, y = Inf, label = text, hjust = hjust), vjust = 1.15, size = PT,
              lineheight = 0.95, colour = INK, inherit.aes = FALSE) +
    facet_wrap(~strip) +
    scale_x_continuous(breaks = if (ages < 8) sort(unique(d$age)) else scales::breaks_pretty(4)) +
    scale_y_continuous(labels = MINUS_LABELS, breaks = scales::breaks_pretty(4), expand = expansion(mult = c(0.05, 0.4))) +
    labs(x = paste0("Age (", d$age_unit[1], ")"), y = y_title) +
    theme_fig() +
    BOXED
}
scatters <- Map(cohort_scatter, split(b, b$cohort), y_title = list("Clock score", NULL))
b <- wrap_plots(scatters, nrow = 2) + FLUSH

# c: pooled rho of every weight set in each set of cohorts (bold, m6Brain); rows grouped by whether the weights were
# derived from age effects of the genes (supervised: observed effects, aging gene sets, published proteome clocks) or
# not (zero-shot); random weights, the median of 1,000 draws
c_ <- read_panel("fig5c")
c_[, `:=`(column = factor(SET_LABELS[paste(set, species)], SET_LABELS[names(SETS)]),
          row = factor(label, rev(unique(label))),
          group = factor(group, c("supervised", "zero-shot"), c("Supervised", "Zero-shot")))]
cc <- ggplot(c_, aes(column, row, fill = rho)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = MINUS(sprintf("%.2f", rho)), fontface = fifelse(weights == "model", "bold", "plain")),
            size = PT, colour = INK) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_div(0.9, "Pooled Spearman ρ of the score with age") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_heat() +
  theme(legend.position = "bottom", panel.spacing.y = unit(4, "mm"),
        strip.text.y = element_text(angle = -90, margin = margin(l = 2, unit = "mm")))

# d: the clock from the weights of the unseen genes alone (no protein age effect in any discovery cohort) in the
# external cohorts with at least 200 of them quantified (point and 95% interval; filled where P < 0.05), the
# random-effects estimate (diamond) and the pooled rho of 1,000 permutations of the weights among the same genes
# (violin) with the share of them that reach the observed pooled rho (P)
u <- read_panel("fig5d")
u[, pooled := kind == "pooled"]
u[, label := fifelse(pooled, "Pooled", sprintf("%s, n = %d, %d genes", cohort, n, n_genes))]
NULL_ROW <- "Permuted weights"
u[, `:=`(row = factor(label, rev(c(label, NULL_ROW))), fill = fifelse(pooled | p < 0.05, CONTRAST, "white"),
         edge = fifelse(pooled | p < 0.05, EDGE, GREY))]
perm <- read_panel("fig5d_null")[, row := factor(NULL_ROW, levels(u$row))]
notes <- u[pooled == TRUE, .(row = factor(c(label, NULL_ROW), levels(u$row)), x = c(-0.43, 0.9), hjust = c(0, 1),
                             text = c(sprintf("ρ = %.2f, I² = %.0f%%", rho, i2), fmt_p(p_random)))]
d <- ggplot(u, aes(rho, row)) +
  ZERO_X +
  geom_violin(data = perm, orientation = "y", fill = LIGHT, colour = GREY, linewidth = 0.5, width = 0.8) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = CI, linewidth = 0.6) +
  geom_point(data = ~ .x[pooled == FALSE], aes(fill = I(fill), colour = I(edge)), shape = 21, size = 3.6, stroke = 0.5) +
  geom_point(data = ~ .x[pooled == TRUE], aes(fill = I(fill)), colour = EDGE, shape = 23, size = 5.2, stroke = 0.5) +
  geom_label(data = notes, aes(x, row, label = text, hjust = hjust), size = PT, colour = INK, fill = "white",
             border.colour = NA, label.padding = unit(0.6, "mm"), inherit.aes = FALSE) +
  scale_x_continuous(labels = MINUS_LABELS, breaks = seq(-0.25, 0.75, 0.25), limits = c(-0.45, 0.9)) +
  labs(x = "Spearman ρ of the unseen-gene score with age", y = NULL) +
  theme_rows()

fig <- panel_grid(list(a = a, b = b, c = cc, d = d), design = "AB\nAC\nDC", widths = c(156, 210),
                  heights = c(125, 95, 70))
save_fig(fig, "fig5", 290)
