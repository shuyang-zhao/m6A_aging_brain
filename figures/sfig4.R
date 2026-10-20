# Additional file 1 Fig. S4: details of the clock, rho of every weight set in every cohort and pooled per set, the
# scores of the model without the m6A features against the model's and the random-weight null per set
# reads source_data/sfig4a-c; writes output/figures/sfig4.png and .pdf
source("figures/style.R")

SETS <- c("held-out|mouse" = "Held-out mouse", "external|mouse" = "External mouse", "held-out|human" = "Held-out human")
SET_COL <- c("Held-out mouse" = FOCAL, "External mouse" = CONTRAST, "Held-out human" = RED)
GROUPS <- c(supervised = "Supervised", "zero-shot" = "Zero-shot", ablation = "Ablation")
SUPERVISED <- c("Observed z of other cohorts", "MSigDB DEMAGALHAES_AGING", "MSigDB LU_AGING_BRAIN")
ZERO_SHOT <- c("m6Brain", "Length only", "Linear m6A load", "Random weights")
tag_set <- function(d) d[, set := factor(SETS[paste(set, species, sep = "|")], SETS)]

# a: Spearman rho of the score of every weight set with age in every cohort (ordered by m6Brain within the set) and
# pooled over the cohorts of each set (bold, m6Brain); rows grouped as in Fig. 5c, published clocks after the gene
# sets; random weights, the median of 1,000 draws
d <- tag_set(read_panel("sfig4a"))
d[, pooled := kind == "pooled"]
order_of <- d[weights == "m6Brain" & !pooled][order(set, -rho), cohort]
clocks <- sort(setdiff(d[group == "supervised", unique(weights)], SUPERVISED))
d[, `:=`(block = factor(fifelse(pooled, "Pooled", as.character(set)), c(SETS, "Pooled")),
         column = factor(fifelse(pooled, as.character(set), cohort), c(order_of, SETS)),
         row = factor(weights, rev(c(SUPERVISED, clocks, ZERO_SHOT, "Without m6A features"))),
         group = factor(GROUPS[group], GROUPS))]
a <- ggplot(d, aes(column, row, fill = rho)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = cell_number(rho), fontface = fifelse(weights == "m6Brain", "bold", "plain")),
            size = PT, colour = INK) +
  facet_grid(group ~ block, scales = "free", space = "free", switch = "y") +
  scale_fill_gradientn(
    colours = DIVERGING, limits = c(-1, 1), name = "Spearman ρ of the score with age", breaks = c(-1, 0, 1),
    labels = MINUS_LABELS, guide = guide_colourbar(
      theme = theme(legend.key.width = unit(26, "mm"), legend.key.height = unit(3.2, "mm")),
      title.position = "top", title.hjust = 0.5, ticks = FALSE, frame.colour = NA)
  ) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_heat() +
  theme(legend.position = "bottom", strip.placement = "outside", panel.spacing = unit(3, "mm"),
        axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, margin = margin(b = 2, unit = "mm")),
        strip.text.y.left = element_text(angle = 0, hjust = 1, margin = margin(r = 3, unit = "mm")))

# b: correlation of the scores without the m6A features with the scores of the model in every cohort
r <- tag_set(read_panel("sfig4b"))
r[, cohort := factor(cohort, order_of)]
b <- ggplot(r, aes(cohort, r_model)) +
  geom_point(aes(fill = set), shape = 21, colour = EDGE, size = 3.6, stroke = 0.5) +
  facet_grid(. ~ set, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = SET_COL, guide = "none") +
  scale_x_discrete(expand = c(0, 0.6)) +
  scale_y_continuous(limits = c(0.98, 1.0005), breaks = c(0.98, 0.99, 1)) +
  labs(x = NULL, y = "r of scores without m6A vs m6Brain") +
  theme_fig() +
  BOXED +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.spacing = unit(3, "mm"))

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
  scale_x_continuous(labels = MINUS_LABELS, breaks = seq(-0.4, 0.8, 0.4), limits = c(NA, 1.15)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Pooled ρ of random weights", y = "Density") +
  theme_fig() +
  BOXED +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

fig <- panel_grid(list(a = a, b = b, c = c_), design = "AA\nBC", widths = c(190, 176), heights = c(175, 125))
save_fig(fig, "sfig4", 300)
