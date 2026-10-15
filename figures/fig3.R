# Figure 3: the m6A-protein association in human brain
# reads source_data/fig3a, fig3a_given_mrna, fig3b, fig3b_units, fig3c and fig3d and the main model of fig2d; writes
# output/figures/fig3.png and .pdf
source("figures/style.R")

# a: age effect of human genes by number of m6A peaks of the brain map: protein, mean over units with a 95% interval
# across units (number, the load coefficient over cohorts as cited in the text); protein given mRNA, the paired cohort
a <- sites_panel("fig3a", "Human", shapes = c(21, 21), open = "Luo 2024", per = "/decade",
                 axes = labs(x = "m6A peaks per gene (brain map)", y = "Age effect (% per decade)"))

# b: m6A load in the human brain maps: one dot per unit, the random-effects estimate over cohorts with its 95%
# interval (number of cohorts in parentheses); units beyond the axis are arrows at its edge
HUMAN <- c(human_protein = "Protein", human_striatum = "Striatum mRNA", human_cortex = "Cortex mRNA",
           human_hippocampus = "Hippocampus mRNA", human_cerebellum = "Cerebellum mRNA",
           human_neurons = "Neuronal mRNA")
re <- read_panel("fig3b")
re[, `:=`(label = sprintf("%s (%d)", HUMAN[context], n_cohorts), y = 7 - match(context, names(HUMAN)),
          key = fcase(context == "human_protein", "Protein", context == "human_neurons", "Single-nucleus RNA",
                      default = "Bulk RNA"))]
units <- read_panel("fig3b_units")
setorder(units, unit_id, context)
XLIM <- c(-1.2, 0.6)
set.seed(1)
units[, `:=`(y = 7 - match(context, names(HUMAN)) + runif(.N, -0.17, 0.17), key = re$key[match(context, re$context)])]
b <- ggplot(re, aes(effect, y)) +
  ZERO_X +
  geom_point(data = units[effect >= XLIM[1]], aes(fill = I(LAYER[key])), colour = EDGE, shape = 21, size = 2.4,
             stroke = 0.5) +
  geom_segment(data = units[effect < XLIM[1]], aes(x = XLIM[1] + 0.1, xend = XLIM[1] + 0.01, yend = y,
                                                       colour = I(LAYER[key])),
               linewidth = 0.8, arrow = arrow(length = unit(1.4, "mm"), type = "closed")) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = CI, linewidth = 0.6) +
  geom_point(aes(fill = I(LAYER[key])), colour = EDGE, shape = 23, size = 5.2, stroke = 0.5) +
  scale_x_continuous(labels = MINUS_LABELS, limits = XLIM, expand = expansion(add = 0.02)) +
  scale_y_continuous(breaks = re$y, labels = re$label, expand = expansion(add = 0.5)) +
  labs(x = "Age effect of the m6A load per log unit (% per decade)", y = NULL) +
  theme_rows()

# c: validation tiers of the m6A tests; dot size and colour, z of each validation set; filled where the one-sided test
# passes with the discovery sign, open where it does not
TESTS <- c("m6A load (brain map)", "m6A signature")
VALID <- c(human_protein = "Human protein", human_striatum = "Human striatum", human_cortex = "Human cortex",
           human_hippocampus = "Human hippocampus", human_cerebellum = "Human cerebellum",
           human_neurons = "Human neurons", "other cohorts" = "Other cohorts")
TIERS <- c(mouse_to_human = "Mouse discovery", within_human = "BrainGVEX\ndiscovery")
d <- read_panel("fig3c")[test %in% TESTS]
d[, `:=`(test = factor(test, TESTS, c("Load", "Signature")),
         row = factor(VALID[validation], rev(VALID)), tier = factor(TIERS[tier], TIERS))]
cc <- ggplot(d, aes(test, row)) +
  geom_point(aes(size = abs(validation_z), colour = validation_z, fill = fifelse(passed, validation_z, NA_real_)),
             shape = 21, stroke = 0.5) +
  facet_grid(tier ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_div(4, "Validation z") +
  scale_colour_gradientn(colours = DIVERGING, limits = c(-4, 4), oob = scales::squish, guide = "none") +
  scale_size(range = c(3.6, 13), limits = c(0, 5), guide = "none") +
  scale_x_discrete(position = "top") +
  theme_matrix() +
  BOXED +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90, lineheight = 0.95),
        legend.position = "bottom")

# d: human m6A load counted from the brain map or from the cell-line maps (bars, 95% CI)
MAP <- c("Brain map" = FOCAL, "Cell-line maps" = CONTRAST)
m <- read_panel("fig3d")
m[, `:=`(row = factor(HUMAN[context], rev(HUMAN)),
         def = factor(definition, c("brain peaks", "cell line"), names(MAP)))]
dd <- ggplot(m, aes(effect_pct, row, fill = def)) +
  ZERO_X +
  geom_col(position = position_dodge(0.8, orientation = "y", reverse = TRUE), width = 0.75, colour = EDGE,
           linewidth = 0.5) +
  geom_linerange(aes(xmin = ci_low_pct, xmax = ci_high_pct),
                 position = position_dodge(0.8, orientation = "y", reverse = TRUE), linewidth = 0.6, colour = EDGE) +
  scale_fill_manual(values = MAP, name = NULL) +
  scale_x_continuous(labels = MINUS_LABELS) +
  labs(x = "Age effect of the m6A load per log unit (% per decade)", y = NULL) +
  theme_rows() +
  theme(legend.position = "bottom")

fig <- panel_grid(list(a = a, b = b, c = cc, d = dd), design = "AC\nBD", widths = c(100, 83), heights = c(78, 80))
save_fig(fig, "fig3", 316)
