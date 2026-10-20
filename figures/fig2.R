# Figure 2: gene-level model of the age effect from transcript features
# reads source_data/fig2a, fig2b, fig2c, fig2c_permuted and fig2d; writes output/figures/fig2.png and .pdf
source("figures/style.R")
suppressPackageStartupMessages(library(ggbeeswarm))

model_scales <- function(levels) {
  i <- match(levels, MODELS)
  list(
    scale_fill_manual(values = setNames(MODEL_FILL[i], levels), name = NULL),
    scale_colour_manual(values = setNames(MODEL_EDGE[i], levels), name = NULL),
    scale_shape_manual(values = setNames(MODEL_SHAPE[i], levels), name = NULL),
    guides(fill = guide_legend(nrow = 1), colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
  )
}
rho_panel <- function(d, levels, targets) {
  n <- d[, .(n = n[1]), by = target]
  d[, `:=`(target = factor(target, rev(targets)), model = factor(model, levels))]
  labels <- n[, setNames(sprintf("%s\nn = %s", target, format(n, big.mark = ",", trim = TRUE)), target)]
  dodge <- position_dodge(0.78, reverse = TRUE)
  ggplot(d, aes(rho, target, group = model)) +
    ZERO_X +
    geom_linerange(aes(xmin = rho_ci_low, xmax = rho_ci_high), colour = CI, linewidth = 0.6, position = dodge) +
    geom_point(aes(fill = model, colour = model, shape = model), size = 3.4, stroke = 0.6, position = dodge) +
    model_scales(levels) +
    scale_y_discrete(labels = labels) +
    scale_x_continuous(labels = MINUS_LABELS) +
    labs(x = "Held-out Spearman ρ", y = NULL) +
    theme_rows() +
    theme(legend.position = "bottom", legend.key.spacing.x = unit(3, "mm"),
          axis.text.y = element_text(lineheight = 0.95))
}

# a: held-out Spearman rho of the four models on the five targets (95% bootstrap CI), genes of the held-out chromosomes
a <- rho_panel(read_panel("fig2a"), MODELS[1:4], TARGETS)

# b: loss of held-out R2 when a feature group is removed (bars), mouse and human protein; for the tested group the
# permutation P against 100 models with its features permuted across genes (the 95% range of these models, below
# 0.001, is not drawn); m6A groups on a grey band
b <- read_panel("fig2b")
b[, target := factor(target, c("Mouse protein", "Human protein"))]
order <- b[, .(m = mean(delta_r2)), by = block][order(m), block]
b[, `:=`(row = factor(block, order), tested = !is.na(p))]
bands <- data.table(block = order, i = seq_along(order))[grepl("m6A", block)]
TARGET_COL <- c("Mouse protein" = FOCAL, "Human protein" = RED)
dodge_b <- position_dodge(0.78, reverse = TRUE)
bb <- ggplot(b, aes(delta_r2, row, group = target)) +
  geom_rect(data = bands, aes(xmin = -Inf, xmax = Inf, ymin = i - 0.5, ymax = i + 0.5), fill = "#F1F1F1",
            inherit.aes = FALSE) +
  ZERO_X +
  geom_col(aes(fill = target), position = dodge_b, width = 0.7, colour = EDGE, linewidth = 0.5) +
  geom_text(data = b[tested == TRUE], aes(x = Inf, label = paste0(sub(" protein", "", target), ", ", fmt_p(p))),
            hjust = 1.05, size = PT, position = dodge_b) +
  scale_fill_manual(values = TARGET_COL, name = NULL) +
  scale_y_discrete(limits = order) +
  scale_x_continuous(labels = MINUS_LABELS, expand = expansion(mult = c(0.1, 0.25))) +
  labs(x = "Loss of held-out R² on removing the group", y = NULL) +
  theme_rows() +
  theme(legend.position = "bottom")

# c: in silico deletion of the m6A sites: mean change of the prediction of the genes with sites (observed, diamond)
# against 100 copies of the sites permuted across genes (dots, line: 95% of the copies), all sites and per region
r <- read_panel("fig2c")
perm <- read_panel("fig2c_permuted")
REGIONS <- c("All sites", "5' UTR", "CDS", "3' UTR")
REGION_LABEL <- c("All sites", "5′ UTR", "CDS", "3′ UTR")
prep <- function(x) {
  x[, `:=`(target = factor(target, names(TARGET_COL)), region = factor(region, rev(REGIONS), rev(REGION_LABEL)))]
}
prep(r)
prep(perm)
cc <- ggplot(r, aes(y = region)) +
  ZERO_X +
  geom_quasirandom(data = perm, aes(effect), orientation = "y", width = 0.3, size = 1.1, colour = "#8C8C8C",
                   alpha = 0.7, stroke = 0) +
  geom_linerange(aes(xmin = null_low, xmax = null_high), colour = INK, linewidth = 0.6) +
  geom_point(aes(effect, fill = target), shape = 23, size = 5, colour = EDGE, stroke = 0.6) +
  geom_text(aes(x = Inf, label = fmt_p(p)), hjust = 1.05, size = PT) +
  facet_wrap(~target, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = TARGET_COL, guide = "none") +
  scale_x_continuous(labels = MINUS_LABELS, expand = expansion(mult = c(0.06, 0.5)),
                     breaks = \(l) Filter(\(x) x <= 0.001, scales::breaks_pretty(4)(l))) +
  labs(x = "Prediction with minus without sites (age z)", y = NULL) +
  theme_rows() +
  BOXED +
  theme(panel.spacing = unit(6, "mm"))

# d: mouse protein model applied to human genes (rho with the human targets, 95% bootstrap CI) against the shuffled target
d <- rho_panel(read_panel("fig2d"), MODELS[c(5, 4)], c("Human protein", "Human mRNA"))

top <- max(112, 40 + 10.5 * length(order)) # panel b grows with the number of feature groups
fig <- panel_grid(list(a = a, b = bb, c = cc, d = d),
                  design = c(area(1, 1, 1, 180), area(1, 181, 1, 366), area(2, 1, 2, 250), area(2, 251, 2, 366)),
                  heights = c(top, 84))
save_fig(fig, "fig2", top + 84)
