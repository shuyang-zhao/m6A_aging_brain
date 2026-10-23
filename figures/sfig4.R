# Additional file 1 Fig. S4: details of the clock; a, r of every weight set in every cohort and pooled per set; b, the
# label-free mode; c, unseen proteins; d, the error across cohorts against the training-mean null
# reads source_data/sfig4a-d, sfig4b_null and sfig4c_null (scripts 52 and 55); writes output/figures/sfig4.png and .pdf
source("figures/style.R")

SETS <- c("held-out|mouse" = "Held-out mouse", "external|mouse" = "External mouse", "held-out|human" = "Held-out human")
SET_COL <- c("Held-out mouse" = FOCAL, "External mouse" = CONTRAST, "Held-out human" = RED)
GROUPS <- c("m6Brain clock and references" = "m6Brain", "Proteomic clocks" = "Proteomic",
            "Transcriptome clocks" = "Transcriptome", "Signatures and half-lives" = "Signatures",
            "Label-free and gene sets" = "Label-free")
tag_set <- function(d) d[, set := factor(SETS[paste(set, species, sep = "|")], SETS)]

# a: Pearson r of every weight set with age in every cohort (ordered by the m6Brain clock within the set) and pooled
# over the cohorts of each set; rows grouped as in Fig. 5c, then the label-free weight sets (random weights, the median
# of 1,000 draws), ordered by the mean pooled r
d <- tag_set(read_panel("sfig4a"))
d[, pooled := kind == "pooled"]
order_of <- d[weights == "m6brain_clock" & !pooled][order(set, -r), cohort]
d[, mean_r := mean(r[pooled]), by = label]
rows <- unique(d[order(match(group, names(GROUPS)), -(weights == "m6brain_clock"), -mean_r), label])
d[, `:=`(block = factor(fifelse(pooled, "Pooled", as.character(set)), c(SETS, "Pooled")),
         column = factor(fifelse(pooled, as.character(set), cohort), c(unique(order_of), SETS)),
         row = factor(label, rev(rows)), group = factor(GROUPS[group], GROUPS))]
a <- ggplot(d, aes(column, row, fill = r)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = cell_number(r), fontface = fifelse(weights == "m6brain_clock", "bold", "plain")),
            size = PT, colour = INK) +
  facet_grid(group ~ block, scales = "free", space = "free") +
  scale_fill_div(1, "Pearson r of the score with age") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_heat() +
  theme(legend.position = "bottom", panel.spacing = unit(2.4, "mm"),
        axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, margin = margin(b = 2, unit = "mm")),
        strip.text.y = element_text(angle = -90, margin = margin(l = 2, unit = "mm")))

# b: label-free mode, pooled r (95% CI) of the m6Brain predictions used as weights and of the comparison weights per
# set; violin, the pooled r of 1,000 random weights, with the P of the m6Brain predictions against them
NULL_ROW <- "Random weights"
lf <- tag_set(read_panel("sfig4b"))
lf[, row := factor(label, rev(c("m6Brain label-free", "Label-free without m6A features", "Linear m6A load",
                                "Length only", NULL_ROW)))]
random <- tag_set(read_panel("sfig4b_null"))[, row := factor(NULL_ROW, levels(lf$row))]
note <- lf[weights == "model", .(set, row = factor(NULL_ROW, levels(lf$row)), text = fmt_p(p_random))]
b <- ggplot(lf, aes(r, row)) +
  ZERO_X +
  geom_violin(data = random, orientation = "y", fill = LIGHT, colour = GREY, linewidth = 0.5, width = 0.8) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = CI, linewidth = 0.6) +
  geom_point(aes(fill = fifelse(weights == "model", SET_COL[as.character(set)], "white")), shape = 21, colour = EDGE,
             size = 3.4, stroke = 0.5) +
  geom_text(data = note, aes(x = -0.72, label = text), hjust = 0, vjust = 2.6, size = PT) +
  facet_wrap(~set, nrow = 1) +
  scale_fill_identity() +
  scale_x_continuous(labels = MINUS_LABELS, breaks = c(-0.5, 0, 0.5, 1), limits = c(-0.75, 1)) +
  scale_y_discrete(expand = expansion(add = c(1.1, 0.6))) +
  labs(x = "Pooled Pearson r with age", y = NULL) +
  theme_rows() +
  BOXED +
  theme(panel.spacing.x = unit(3, "mm"))

# c: unseen proteins, the label-free score from the weights of the genes without a protein age effect in any discovery
# cohort in the external deep-DIA cohorts (95% CI; filled, P < 0.05), the random-effects estimate (diamond) and the
# pooled r of 1,000 permutations of the weights among the same genes (violin), with their P
u <- read_panel("sfig4c")
u[, pooled := kind == "pooled"]
u[, label := fifelse(pooled, "Pooled", sprintf("%s, n = %d,\n%d genes", cohort, n, n_genes))]
PERMUTED <- "Permuted weights"
u[, `:=`(row = factor(label, rev(c(label, PERMUTED))), fill = fifelse(pooled | p < 0.05, CONTRAST, "white"))]
perm <- read_panel("sfig4c_null")[, row := factor(PERMUTED, levels(u$row))]
notes <- u[pooled == TRUE, .(row = factor(PERMUTED, levels(u$row)), text = fmt_p(p_random))]
cc <- ggplot(u, aes(r, row)) +
  ZERO_X +
  geom_violin(data = perm, orientation = "y", fill = LIGHT, colour = GREY, linewidth = 0.5, width = 0.8) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = CI, linewidth = 0.6) +
  geom_point(aes(fill = I(fill), shape = pooled, size = pooled), colour = EDGE, stroke = 0.5) +
  geom_text(data = notes, aes(x = -0.28, label = text), hjust = 0, vjust = 2.6, size = PT) +
  scale_shape_manual(values = c("TRUE" = 23, "FALSE" = 21), guide = "none") +
  scale_size_manual(values = c("TRUE" = 4.6, "FALSE" = 3.4), guide = "none") +
  scale_x_continuous(labels = MINUS_LABELS, breaks = c(-0.2, 0, 0.2, 0.4, 0.6), limits = c(-0.3, 0.62)) +
  scale_y_discrete(expand = expansion(add = c(1.1, 0.6))) +
  labs(x = "r of the unseen-gene score", y = NULL) +
  theme_rows()

# d: MAE of the m6Brain clock in every cohort (filled) and of the training-mean null (open), against the offset of the
# cohort's mean age from the mean age of its training samples; mouse in months, human in years
m <- tag_set(read_panel("sfig4d"))
m[, panel := factor(fifelse(species == "mouse", "Mouse (months)", "Human (years)"), c("Mouse (months)",
                                                                                      "Human (years)"))]
dd <- ggplot(m, aes(age_offset)) +
  geom_segment(aes(y = mae_null, yend = mae), colour = CI, linewidth = 0.5) +
  geom_point(aes(y = mae_null, fill = "Training-mean null"), shape = 21, colour = EDGE, size = 2.8, stroke = 0.5) +
  geom_point(aes(y = mae, fill = set), shape = 21, colour = EDGE, size = 3.2, stroke = 0.5) +
  facet_wrap(~panel, scales = "free") +
  scale_fill_manual(values = c(SET_COL, "Training-mean null" = "white"), breaks = c(unname(SETS), "Training-mean null"),
                    name = NULL, guide = guide_legend(nrow = 2)) +
  scale_x_continuous(labels = MINUS_LABELS, breaks = scales::breaks_pretty(4)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Cohort mean age − training mean age", y = "MAE") +
  theme_fig() +
  BOXED +
  theme(legend.position = "bottom")

fig <- panel_grid(list(a = a, b = b, c = cc, d = dd), design = "AAA\nBCD", widths = c(158, 98, 110),
                  heights = c(258, 82))
save_fig(fig, "sfig4", 340)
