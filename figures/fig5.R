# Figure 5: the proteome clock in neurodegeneration
# reads source_data/fig5a, fig5b and fig5c (script 53); writes output/figures/fig5.png and .pdf
source("figures/style.R")

tests <- read_panel("fig5c")
CLOCK_Y <- "Clock score"

# a: AD, donor scores of cases (red) and controls (grey) in two studies (box, median and quartiles); above, the
# coefficient of AD on the score (SD; adjusted for age, sex and batch) with its P
a <- read_panel("fig5a")
note <- tests[disease == "AD" & study %in% a$study, .(strip = factor(study, unique(a$study)),
                                                      text = paste0(signed(est), " SD\n", fmt_p(p)))]
a[, `:=`(group = factor(group, c("Control", "AD")), strip = factor(study, unique(study)))]
a <- ggplot(a, aes(group, score)) +
  geom_boxplot(aes(colour = group), fill = NA, width = 0.5, linewidth = 0.6, outliers = FALSE, show.legend = FALSE) +
  geom_point(aes(fill = group), colour = EDGE, shape = 21, size = 2.8, stroke = 0.4, alpha = 0.9,
             position = position_jitter(width = 0.14, height = 0, seed = 1)) +
  geom_text(data = note, aes(x = -Inf, y = Inf, label = text), hjust = -0.08, vjust = 1.15, size = PT, lineheight = 0.95,
            colour = INK, inherit.aes = FALSE) +
  facet_wrap(~strip) +
  scale_fill_manual(values = c(Control = GREY, AD = RED), guide = "none") +
  scale_colour_manual(values = c(Control = EDGE, AD = EDGE)) +
  scale_y_continuous(labels = MINUS_LABELS, expand = expansion(mult = c(0.05, 0.3))) +
  labs(x = NULL, y = CLOCK_Y) +
  theme_fig() +
  BOXED

# b: HD allelic series at 6 months, score of each mouse against the CAG length of its genotype (red, expanded alleles;
# grey, 7 and 20 repeats); line, least squares
b <- read_panel("fig5b")
b[, expanded := cag > 20]
b <- ggplot(b, aes(cag, score)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = EDGE, linewidth = 0.8) +
  geom_point(aes(fill = expanded), colour = EDGE, shape = 21, size = 2.8, stroke = 0.4, alpha = 0.9,
             position = position_jitter(width = 4, height = 0, seed = 1)) +
  scale_fill_manual(values = c("TRUE" = RED, "FALSE" = GREY), guide = "none") +
  scale_x_continuous(breaks = c(0, 100, 200), expand = expansion(add = 12)) +
  scale_y_continuous(labels = MINUS_LABELS) +
  labs(x = sprintf("CAG repeats of the genotype (%d-month-old mice)", b$age[1]), y = CLOCK_Y) +
  theme_fig() +
  BOXED

# c: coefficient of the disease on the score (SD; point and 95% CI; red, P < 0.05) for each study, AD studies combined
# (diamond), HD per 100 CAG repeats; P on the right
c_rows <- c("Ping 2020", "Merrihew 2023", "AD combined", "Blumenreich 2024", "Umoh 2018", "Langfelder 2016")
d <- copy(tests)
d[, label := fifelse(study == "AD combined", "Combined", study)]
d[, `:=`(row = factor(label, rev(label[match(c_rows, study)])), sig = p < 0.05,
         disease = factor(disease, c("AD", "PD", "ALS", "HD"), c("AD", "PD", "ALS", "HD, per\n100 CAG")))]
cc <- ggplot(d, aes(est, row)) +
  ZERO_X +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = CI, linewidth = 0.6) +
  geom_point(data = ~ .x[study != "AD combined"], aes(fill = sig), colour = EDGE, shape = 21, size = 3.8, stroke = 0.5) +
  geom_point(data = ~ .x[study == "AD combined"], aes(fill = sig), colour = EDGE, shape = 23, size = 5.2, stroke = 0.5) +
  geom_text(aes(x = Inf, label = fmt_p(p)), hjust = -0.12, size = PT, colour = INK) +
  facet_grid(disease ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = c("TRUE" = RED, "FALSE" = "white"), guide = "none") +
  scale_x_continuous(labels = MINUS_LABELS, breaks = scales::breaks_pretty(4)) +
  coord_cartesian(clip = "off") +
  labs(x = "Coefficient on the clock score (SD)", y = NULL) +
  theme_rows() +
  BOXED +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0), panel.spacing.y = unit(3, "mm"),
        plot.margin = margin(10.4, 28, 3, 2, "mm"))

fig <- panel_grid(list(a = a, b = b, c = cc), nrow = 1, widths = c(100, 105, 161))
save_fig(fig, "fig5", 100)
