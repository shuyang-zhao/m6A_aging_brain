# Additional file 1 Fig. S5: the clock in disease
# a, age acceleration of every sample by region, cases and controls; b, the 30 clocks scored on the same donors; c, the
# acceleration against the individual m6A index in aging cohorts
# reads source_data/sfig5a, sfig5a_tests, sfig5b and sfig5c (scripts 53 and 54); writes output/figures/sfig5.png and .pdf
source("figures/style.R")

STUDIES <- c("Ping 2020", "Merrihew 2023", "AD combined", "Blumenreich 2024", "Umoh 2018", "Langfelder 2016")

# a: cases (red) and controls (grey) of each region (box, median and quartiles), with the coefficient of the disease in
# the region (years) and its P
s <- read_panel("sfig5a")
tests <- read_panel("sfig5a_tests")
strip_of <- function(study, region) paste0(study, "\n", vapply(region, function(r) paste(strwrap(r, 22), collapse = "\n"),
                                                               ""))
s[, `:=`(strip = strip_of(study, region),
         group = factor(fifelse(group == "Case", disease, "Control"), c("Control", "AD", "PD", "ALS")))]
s[, strip := factor(strip, unique(strip))]
tests[, `:=`(strip = factor(strip_of(study, region), levels(s$strip)),
             text = paste0(MINUS(sprintf("%+.1f", est)), " y\n", fmt_p(p)))]
a <- ggplot(s, aes(group, acceleration)) +
  ZERO_Y +
  geom_boxplot(aes(colour = group), fill = NA, width = 0.5, linewidth = 0.6, outliers = FALSE, show.legend = FALSE) +
  geom_point(aes(fill = group), colour = EDGE, shape = 21, size = 2.2, stroke = 0.4, alpha = 0.9,
             position = position_jitter(width = 0.14, height = 0, seed = 1)) +
  geom_text(data = tests, aes(x = -Inf, y = Inf, label = text), hjust = -0.06, vjust = 1.15, size = PT,
            lineheight = 0.95, inherit.aes = FALSE) +
  facet_wrap(~strip, nrow = 2, scales = "free_x") +
  scale_fill_manual(values = c(Control = GREY, AD = RED, PD = RED, ALS = RED), guide = "none") +
  scale_colour_manual(values = c(Control = EDGE, AD = EDGE, PD = EDGE, ALS = EDGE)) +
  scale_y_continuous(labels = MINUS_LABELS, expand = expansion(mult = c(0.05, 0.35))) +
  labs(x = NULL, y = "Age acceleration (years)") +
  theme_fig() +
  BOXED

# b: the 30 clocks on the same donors, coefficient of the disease on each score (SD of the score; HD per 100 CAG
# repeats; colour) and the rank of the clock in the study (number; bold, P < 0.05); rows ordered by the AD combined
# coefficient
r <- read_panel("sfig5b")[unit == "SD"]
r[, column := factor(fifelse(study == "Langfelder 2016", "Langfelder 2016\n(HD, per 100 CAG)",
                             fifelse(study == "AD combined", "AD combined", sprintf("%s\n(%s)", study, disease))))]
r[, column := factor(column, unique(column[order(match(study, STUDIES))]))]
r[, row := factor(label, r[study == "AD combined"][order(est), label])]
b <- ggplot(r, aes(column, row, fill = est)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = rank, fontface = fifelse(p < 0.05, "bold", "plain")), size = PT, colour = INK) +
  scale_fill_div(1.5, "Disease coefficient (SD of the score)") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_heat() +
  theme(legend.position = "bottom", axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5))

# c: acceleration per SD of the age-adjusted m6A index per cohort (95% CI; filled, P < 0.05) and the random-effects
# estimate (diamond), held-out mouse and external mouse (months), human (years)
SETS <- c(mouse = "Held-out mouse (months)", external = "External mouse (months)", human = "Held-out human (years)")
COL <- c(mouse = FOCAL, external = CONTRAST, human = RED)
c_ <- read_panel("sfig5c")
c_[, pooled := kind == "pooled"]
c_[, label := fifelse(pooled, "Random effects", sprintf("%s, n = %d", cohort, n))]
c_[, `:=`(panel = factor(SETS[set], SETS), fill = fifelse(pooled | p < 0.05, COL[set], "white"))]
c_[, row := factor(paste(set, label), rev(paste(set, label)[order(match(set, names(SETS)), pooled, est)]))]
cc <- ggplot(c_, aes(est, row)) +
  ZERO_X +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = CI, linewidth = 0.6) +
  geom_point(aes(fill = I(fill), shape = pooled, size = pooled), colour = EDGE, stroke = 0.5) +
  geom_text(data = c_[pooled == TRUE], aes(x = Inf, label = fmt_p(p)), hjust = 1.08, size = PT) +
  facet_grid(panel ~ ., scales = "free", space = "free_y") +
  scale_shape_manual(values = c("TRUE" = 23, "FALSE" = 21), guide = "none") +
  scale_size_manual(values = c("TRUE" = 4.6, "FALSE" = 3.2), guide = "none") +
  scale_y_discrete(labels = function(x) sub("^[a-z]+ ", "", x)) +
  scale_x_continuous(labels = MINUS_LABELS, breaks = scales::breaks_pretty(4)) +
  coord_cartesian(xlim = c(-8, 5)) +
  labs(x = "Acceleration per SD of the age-adjusted m6A index", y = NULL) +
  theme_rows() +
  BOXED +
  theme(strip.text.y = element_text(angle = -90))

fig <- panel_grid(list(a = a, b = b, c = cc), design = "AA\nBC", widths = c(196, 170), heights = c(110, 215))
save_fig(fig, "sfig5", 325)
