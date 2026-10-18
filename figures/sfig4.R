# Additional file 1 Fig. S4: details of the clock, rho of every weight set in every cohort, the scores of the model
# without the m6A features against the model's and the random-weight null per set
# reads source_data/sfig4a-c; writes output/figures/sfig4.png and .pdf
source("figures/style.R")

SETS <- c("held-out|mouse" = "Held-out mouse", "held-out|human" = "Held-out human", "external|mouse" = "External mouse")
SET_COL <- c("Held-out mouse" = FOCAL, "Held-out human" = RED, "External mouse" = CONTRAST)
WEIGHTS <- c("Model", "Without m6A features", "Length only", "Linear m6A load", "MSigDB DEMAGALHAES_AGING",
             "MSigDB LU_AGING_BRAIN", "Observed z of other cohorts")
tag_set <- function(d) d[, set := factor(SETS[paste(set, species, sep = "|")], SETS)]

# a: Spearman rho of the score of every weight set with age in every cohort (cohorts ordered by the model within the
# set); b: correlation of the scores without the m6A features with the scores of the model
d <- tag_set(read_panel("sfig4a"))
order_of <- d[weights == "Model"][order(set, -rho), cohort]
d[, `:=`(cohort = factor(cohort, order_of), weights = factor(weights, WEIGHTS))]
a <- ggplot(d, aes(weights, cohort, fill = rho)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = MINUS(sprintf("%.2f", rho)), fontface = fifelse(weights == "Model", "bold", "plain")),
            size = PT, colour = INK) +
  facet_grid(set ~ ., scales = "free_y", space = "free_y", switch = "y", labeller = label_wrap_gen(8)) +
  scale_fill_gradientn(
    colours = DIVERGING, limits = c(-1, 1), name = "Spearman ρ of the score with age", breaks = c(-1, 0, 1),
    labels = MINUS_LABELS, guide = guide_colourbar(
      theme = theme(legend.key.width = unit(26, "mm"), legend.key.height = unit(3.2, "mm")),
      title.position = "top", title.hjust = 0.5, ticks = FALSE, frame.colour = NA)
  ) +
  scale_x_discrete(position = "top", labels = function(x) sub("_", " ", x), expand = c(0, 0)) +
  scale_y_discrete(limits = rev, expand = c(0, 0)) +
  theme_heat() +
  theme(legend.position = "bottom", strip.placement = "outside",
        axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, margin = margin(b = 2, unit = "mm")),
        strip.text.y.left = element_text(angle = 0, hjust = 1, lineheight = 0.95, margin = margin(r = 2, unit = "mm")))

r <- tag_set(read_panel("sfig4b"))
r[, cohort := factor(cohort, order_of)]
b <- ggplot(r, aes(r_model, cohort)) +
  geom_point(aes(fill = set), shape = 21, colour = EDGE, size = 3.6, stroke = 0.5) +
  facet_grid(set ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = SET_COL, guide = "none") +
  scale_y_discrete(limits = rev, expand = c(0, 0.6)) +
  scale_x_continuous(limits = c(0.98, 1.0005), breaks = c(0.98, 0.99, 1)) +
  labs(x = "r of the scores without m6A vs the model's", y = NULL) +
  theme_fig() +
  BOXED +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), strip.text.y = element_blank(),
        panel.spacing = unit(1.6, "mm"))

# c: pooled rho of the random weights per set (density of the draws; line, the model)
n <- tag_set(read_panel("sfig4c"))
m <- unique(n[, .(set, model_rho, p_random)])
m[, `:=`(label = paste0("ρ = ", sprintf("%.2f", model_rho), "\n", fmt_p(p_random)),
         right = model_rho > 0.5)]
c_ <- ggplot(n, aes(rho)) +
  geom_density(fill = LIGHT, colour = GREY, linewidth = 0.5, adjust = 0.8) +
  geom_vline(data = m, aes(xintercept = model_rho, colour = set), linewidth = 1.2) +
  geom_text(data = m, aes(x = model_rho, y = Inf, label = label, hjust = fifelse(right, 1.12, -0.12)), vjust = 1.3,
            size = PT, lineheight = 0.95, colour = INK) +
  facet_wrap(~set, nrow = 1) +
  scale_colour_manual(values = SET_COL, guide = "none") +
  scale_x_continuous(labels = MINUS_LABELS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Pooled ρ of random weights", y = "Density") +
  theme_fig() +
  BOXED +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

fig <- (a + b + c_ + plot_layout(design = "AB\nCC", widths = c(115, 68), heights = c(132, 32))) +
  plot_annotation(tag_levels = "a", theme = PANEL_TAG)
save_fig(fig, "sfig4", 2 * 164)
