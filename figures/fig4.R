# Figure 4: the zero-shot proteome clock
# reads source_data/fig4a, fig4b, fig4c and fig4c_null (script 52); writes output/figures/fig4.png and .pdf
source("figures/style.R")

SETS <- c("held-out mouse" = FOCAL, "external mouse" = CONTRAST, "held-out human" = RED)
SET_LABELS <- c("held-out mouse" = "Held-out mouse", "external mouse" = "External mouse",
                "held-out human" = "Held-out human")

# a: Spearman rho of the clock score with age in every cohort (point; filled where P < 0.05) with its Fisher-z 95%
# interval, and the random-effects estimate over the cohorts of the set (diamond; I2, the heterogeneity)
a <- read_panel("fig4a")
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
b <- read_panel("fig4b")
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

# c: pooled rho of the model weights and of five control weights (point and 95% interval; model coloured) per set,
# above the 1,000 sets of random weights (violin) and the share of them that reach the model (P)
c_pool <- read_panel("fig4c")
c_pool[, `:=`(group = factor(paste(set, species), names(SETS)))]
CONTROL_ROWS <- c("Model", "Length only", "Linear m6A load", "MSigDB DEMAGALHAES_AGING", "MSigDB LU_AGING_BRAIN",
                  "Observed z of other cohorts")
ROW_LEVELS <- rev(c("Random weights", CONTROL_ROWS))
c_pool[, row := factor(weights, ROW_LEVELS)]
c_pool[, fill := fifelse(weights == "Model", SETS[as.character(group)], GREY)]
null <- read_panel("fig4c_null")[, `:=`(group = factor(paste(set, species), names(SETS)), row = factor("Random weights", ROW_LEVELS))]
null_p <- c_pool[weights == "Model", .(group, row = factor("Random weights", ROW_LEVELS), x = 1.01,
                                       text = fifelse(p_random < 0.0011, "P < 0.001", fmt_p(p_random)))]
cc <- ggplot(c_pool, aes(rho, row)) +
  ZERO_X +
  geom_violin(data = null, aes(x = rho), orientation = "y", fill = LIGHT, colour = GREY, linewidth = 0.5,
              scale = "width", width = 0.8) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = CI, linewidth = 0.6) +
  geom_point(aes(fill = I(fill)), colour = EDGE, shape = 21, size = 3.8, stroke = 0.5) +
  geom_label(data = null_p, aes(x, row, label = text), hjust = 1, size = PT, colour = INK, fill = "white",
             border.colour = NA, label.padding = unit(0.6, "mm"), inherit.aes = FALSE) +
  facet_grid(group ~ ., labeller = labeller(group = function(x) SET_LABELS[x])) +
  scale_x_continuous(labels = MINUS_LABELS, breaks = seq(-0.5, 1, 0.5), limits = c(-0.7, 1.02)) +
  labs(x = "Pooled Spearman ρ of the score with age", y = NULL) +
  theme_rows() +
  BOXED +
  theme(panel.spacing.y = unit(4, "mm"), strip.text.y = element_text(angle = -90))

fig <- panel_grid(list(a = a, b = b, c = cc), design = "AB\nAC", widths = c(156, 210), heights = c(125, 165))
save_fig(fig, "fig4", 290)
