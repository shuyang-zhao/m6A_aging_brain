# Figure 6: the m6Brain clock in neurodegeneration
# reads source_data/fig6a-e (script 53); writes output/figures/fig6.png and .pdf
source("figures/style.R")

tests <- read_panel("fig6c")
ACC_Y <- "Age acceleration (years)"

# a: AD, age acceleration of cases (red) and controls (grey) in two studies (box, median and quartiles); above, the
# coefficient of AD (years; adjusted for age, sex and batch) with its P
a <- read_panel("fig6a")
note <- tests[disease == "AD" & contrast == "disease vs control" & study %in% a$study,
              .(strip = factor(study, unique(a$study)), text = paste0(MINUS(sprintf("%+.1f", est)), " years\n", fmt_p(p)))]
a[, `:=`(group = factor(group, c("Control", "AD")), strip = factor(study, unique(study)))]
a <- ggplot(a, aes(group, acceleration)) +
  ZERO_Y +
  geom_boxplot(aes(colour = group), fill = NA, width = 0.5, linewidth = 0.6, outliers = FALSE, show.legend = FALSE) +
  geom_point(aes(fill = group), colour = EDGE, shape = 21, size = 2.8, stroke = 0.4, alpha = 0.9,
             position = position_jitter(width = 0.14, height = 0, seed = 1)) +
  geom_text(data = note, aes(x = -Inf, y = Inf, label = text), hjust = -0.08, vjust = 1.15, size = PT, lineheight = 0.95,
            colour = INK, inherit.aes = FALSE) +
  facet_wrap(~strip) +
  scale_fill_manual(values = c(Control = GREY, AD = RED), guide = "none") +
  scale_colour_manual(values = c(Control = EDGE, AD = EDGE)) +
  scale_y_continuous(labels = MINUS_LABELS, expand = expansion(mult = c(0.05, 0.3))) +
  labs(x = NULL, y = ACC_Y) +
  theme_fig() +
  BOXED

# b: HD allelic series at 6 months, acceleration of each mouse against the CAG length of its genotype (red, expanded
# alleles; grey, 7 and 20 repeats, below the 36 of the HD range); line, least squares
b <- read_panel("fig6b")
b[, expanded := cag >= 36]
hd <- tests[disease == "HD"]
b <- ggplot(b, aes(cag, acceleration)) +
  ZERO_Y +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = EDGE, linewidth = 0.8) +
  geom_point(aes(fill = expanded), colour = EDGE, shape = 21, size = 2.8, stroke = 0.4, alpha = 0.9,
             position = position_jitter(width = 4, height = 0, seed = 1)) +
  annotate("text", x = -Inf, y = Inf, label = sprintf("%s months per 100 CAG\n%s", signed(hd$est), fmt_p(hd$p)),
           hjust = -0.06, vjust = 1.15, size = PT, lineheight = 0.95) +
  scale_fill_manual(values = c("TRUE" = RED, "FALSE" = GREY), guide = "none") +
  scale_x_continuous(breaks = c(0, 100, 200), expand = expansion(add = 12)) +
  scale_y_continuous(labels = MINUS_LABELS, expand = expansion(mult = c(0.05, 0.3))) +
  labs(x = sprintf("CAG repeats of the genotype (%d-month-old mice)", b$age[1]), y = "Age acceleration (months)") +
  theme_fig() +
  BOXED

# c and e: coefficient of the disease (point and 95% CI; filled, P < 0.05) for each study, AD studies combined
# (diamond), HD per 100 CAG repeats; P on the right; c on the age acceleration (years; HD months) with the slopes per
# Braak stage and CERAD score of the AD studies, e on the individual m6A index (SD)
STUDY_ROWS <- c("Ping 2020", "Merrihew 2023", "AD combined", "Blumenreich 2024", "Umoh 2018", "Langfelder 2016")
GROUP_ROWS <- c(AD = "AD", PD = "PD", ALS = "ALS", HD = "HD, months\nper 100 CAG", "per Braak stage" = "Per Braak\nstage",
                "per CERAD score" = "Per CERAD\nscore")
forest <- function(d, key, x) {
  d <- copy(d)
  d[, `:=`(label = fifelse(study == "AD combined", "Combined", study),
           group = fifelse(contrast %in% names(GROUP_ROWS), contrast, disease))]
  d[, `:=`(row = factor(paste(group, label), rev(paste(group, label)[order(match(group, names(GROUP_ROWS)),
                                                                           match(study, STUDY_ROWS))])),
           sig = p < 0.05, group = factor(GROUP_ROWS[group], GROUP_ROWS))]
  ggplot(d, aes(est, row)) +
    ZERO_X +
    geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = CI, linewidth = 0.6) +
    geom_point(data = ~ .x[study != "AD combined"], aes(fill = sig), colour = EDGE, shape = 21, size = 3.8, stroke = 0.5) +
    geom_point(data = ~ .x[study == "AD combined"], aes(fill = sig), colour = EDGE, shape = 23, size = 5.2, stroke = 0.5) +
    geom_text(aes(x = Inf, label = fmt_p(p)), hjust = -0.12, size = PT, colour = INK) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = c("TRUE" = key, "FALSE" = "white"), guide = "none") +
    scale_x_continuous(labels = MINUS_LABELS, breaks = scales::breaks_pretty(4)) +
    scale_y_discrete(labels = function(x) sub("^.* (Combined|[A-Z][a-z]+ [0-9]{4})$", "\\1", x)) +
    coord_cartesian(clip = "off") +
    labs(x = x, y = NULL) +
    theme_rows() +
    BOXED +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0), panel.spacing.y = unit(3, "mm"),
          plot.margin = margin(10.4, 28, 3, 2, "mm"))
}
cc <- forest(tests, RED, "Age acceleration (years)")

# d: the clock is linear, so the acceleration splits exactly over groups of proteins; contribution of each group (bar,
# 95% CI) to the acceleration of AD (studies combined, years) and of HD (months per 100 CAG repeats), by cell source and
# direction of the weight (declining, negative weight; rising, positive weight); tick, the contribution expected from
# the group's share of the summed absolute weights (total x weight share); right, share of the total
SOURCES <- c(neuron = "Neuron", astrocyte = "Astrocyte", oligodendrocyte = "Oligodendrocyte", microglia = "Microglia",
             ubiquitous = "Ubiquitous", "no cell-type data" = "No cell-type data")
DIRECTION <- c(declining = "Declining (weight < 0)", rising = "Rising (weight > 0)")
EXPECTED <- "Expected from the weight share"
dd <- read_panel("fig6d")
heads <- dd[partition == "total", setNames(fifelse(disease == "HD", sprintf("HD: %s months per 100 CAG", signed(est)),
                                                   sprintf("AD combined: %s years", MINUS(sprintf("%+.1f", est)))), study)]
dd <- dd[partition != "total"]
dd[, c("source", "direction") := tstrsplit(group, " | ", fixed = TRUE)]
dd[, `:=`(row = factor(SOURCES[source], rev(SOURCES)), direction = factor(DIRECTION[direction], DIRECTION),
          head = factor(heads[study], heads[c("AD combined", "Langfelder 2016")]))]
dodge <- position_dodge(0.84, preserve = "single", reverse = TRUE)
d <- ggplot(dd, aes(est, row, fill = direction)) +
  ZERO_X +
  geom_col(width = 0.84, colour = EDGE, linewidth = 0.4, position = dodge) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = EDGE, linewidth = 0.5, position = dodge) +
  geom_errorbar(aes(xmin = expected, xmax = expected, linetype = EXPECTED), width = 0.94, colour = INK, linewidth = 0.9,
                position = dodge, key_glyph = "vline") +
  geom_text(aes(x = Inf, label = MINUS(sprintf("%.0f%%", 100 * share_of_total))), hjust = 1.1, size = PT, colour = INK,
            position = dodge) +
  facet_grid(. ~ head, scales = "free_x") +
  scale_fill_manual(values = setNames(c(FOCAL, RED), DIRECTION), name = NULL) +
  scale_linetype_manual(values = setNames("solid", EXPECTED), name = NULL) +
  scale_x_continuous(labels = MINUS_LABELS, breaks = scales::breaks_pretty(4), expand = expansion(mult = c(0.05, 0.4))) +
  guides(fill = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  labs(x = "Contribution to the acceleration", y = NULL) +
  theme_rows() +
  BOXED +
  theme(legend.position = "bottom", panel.spacing.x = unit(5, "mm"))

e <- forest(read_panel("fig6e"), FOCAL, "Coefficient on the m6A index (SD)")

fig <- panel_grid(list(a = a, b = b, c = cc, d = d, e = e), design = "ABC\nDDE", widths = c(100, 105, 161),
                  heights = c(110, 135))
save_fig(fig, "fig6", 245)
