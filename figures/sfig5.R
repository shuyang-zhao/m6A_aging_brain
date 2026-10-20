# Additional file 1 Fig. S5: the clock in disease
# a, clock score of every sample by region, cases and controls; b, coefficients of the m6A part of the score (model minus
# the model without the m6A features) in disease
# reads source_data/sfig5a and sfig5b (script 53); writes output/figures/sfig5.png and .pdf
source("figures/style.R")

# a: cases (red) and controls (grey) of each region (box, median and quartiles); cases and controls of a region were
# standardized together
s <- read_panel("sfig5a")
s[, `:=`(strip = paste0(study, "\n", vapply(region, function(r) paste(strwrap(r, 22), collapse = "\n"), "")),
         group = factor(fifelse(group == "Case", disease, "Control"), c("Control", "AD", "PD", "ALS")))]
s[, strip := factor(strip, unique(strip))]
a <- ggplot(s, aes(group, score)) +
  geom_boxplot(aes(colour = group), fill = NA, width = 0.5, linewidth = 0.6, outliers = FALSE, show.legend = FALSE) +
  geom_point(aes(fill = group), colour = EDGE, shape = 21, size = 2.2, stroke = 0.4, alpha = 0.9,
             position = position_jitter(width = 0.14, height = 0, seed = 1)) +
  facet_wrap(~strip, nrow = 2, scales = "free_x") +
  scale_fill_manual(values = c(Control = GREY, AD = RED, PD = RED, ALS = RED), guide = "none") +
  scale_colour_manual(values = c(Control = EDGE, AD = EDGE, PD = EDGE, ALS = EDGE)) +
  scale_y_continuous(labels = MINUS_LABELS) +
  labs(x = NULL, y = "Clock score") +
  theme_fig() +
  BOXED

b <- read_panel("sfig5b")
STUDIES <- c("Ping 2020", "Merrihew 2023", "AD combined", "Blumenreich 2024", "Umoh 2018", "Langfelder 2016")
b[, `:=`(key = fifelse(p < 0.05, RED, "white"), combined = study == "AD combined",
         label = sprintf("%s (n = %d)", study, n),
         disease = factor(disease, c("AD", "PD", "ALS", "HD"), c("AD", "PD", "ALS", "HD, per\n100 CAG")))]
b[, label := factor(label, rev(label[match(STUDIES, study)]))]
b <- ggplot(b, aes(est, label)) +
  ZERO_X +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), orientation = "y", colour = CI, linewidth = 0.6) +
  geom_point(aes(fill = I(key), shape = combined, size = combined), colour = EDGE, stroke = 0.5) +
  geom_text(aes(x = Inf, label = paste0("  ", fmt_p(p))), hjust = 0, size = PT, colour = INK) +
  facet_grid(disease ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 23), guide = "none") +
  scale_size_manual(values = c("FALSE" = 3.6, "TRUE" = 4.6), guide = "none") +
  scale_x_continuous(labels = MINUS_LABELS, expand = expansion(mult = 0.08)) +
  labs(x = "Disease coefficient of the m6A part of the score (SD of the score)", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_rows() +
  BOXED +
  theme(plot.margin = margin(4.4, 28, 3, 2, "mm"), strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, margin = margin(r = 2, unit = "mm")))

save_fig(panel_grid(list(a = a, b = b), ncol = 1, heights = c(95, 130)), "sfig5", 235)
