# Additional file 1 Fig. S3: model details, accuracy of every model on every target, the deletion effect by number of
# sites and the feature groups
# reads source_data/sfig3a-c; writes output/figures/sfig3.png and .pdf
source("figures/style.R")

TRANSFER <- c("Human protein (mouse model)", "Human mRNA (mouse model)")
SPECIES <- function(target) fifelse(grepl("^Mouse", target), FOCAL, RED)

# a: held-out Spearman rho of each model (point, 95% CI); last rows, the mouse protein model applied to the human
# genes
a <- read_panel("sfig3a")
ROW_LEVELS <- c(TARGETS, TRANSFER)
a[, `:=`(row = match(target, rev(ROW_LEVELS)), model = factor(model, MODELS[1:4]))]
a[, y := row + c(0.27, 0.09, -0.09, -0.27)[as.integer(model)]]
a <- ggplot(a, aes(rho, y)) +
  ZERO_X +
  geom_linerange(aes(xmin = rho_ci_low, xmax = rho_ci_high), orientation = "y", colour = CI, linewidth = 0.6) +
  geom_point(aes(colour = model, fill = model, shape = model), size = 3.4, stroke = 0.5) +
  scale_colour_manual(values = setNames(MODEL_EDGE[1:4], MODELS[1:4]), name = NULL) +
  scale_fill_manual(values = setNames(MODEL_FILL[1:4], MODELS[1:4]), name = NULL) +
  scale_shape_manual(values = setNames(MODEL_SHAPE[1:4], MODELS[1:4]), name = NULL) +
  scale_y_continuous(breaks = rev(seq_along(ROW_LEVELS)), labels = ROW_LEVELS, expand = expansion(add = 0.5)) +
  scale_x_continuous(labels = MINUS_LABELS) +
  labs(x = "Held-out Spearman ρ", y = NULL) +
  theme_rows() +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0)) +
  guides(colour = guide_legend(nrow = 1), fill = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))

# b: prediction with minus without the sites of a gene (negative: the sites make the predicted decline stronger), by
# its number of m6A sites (mean, 95% CI)
b <- read_panel("sfig3b")
b[, `:=`(sites = factor(sites, c("1", "2-3", "4-7", "≥8"), c("1", "2–3", "4–7", "≥8")),
         target = factor(target, c("Mouse protein", "Human protein")))]
dodge <- position_dodge(0.3)
b <- ggplot(b, aes(sites, mean, group = target)) +
  ZERO_Y +
  geom_linerange(aes(ymin = ci_low, ymax = ci_high), colour = CI, linewidth = 0.6, position = dodge) +
  geom_line(aes(colour = target), linewidth = 1, position = dodge) +
  geom_point(aes(fill = target), colour = EDGE, shape = 21, size = 3.4, stroke = 0.5, position = dodge) +
  scale_colour_manual(values = c("Mouse protein" = FOCAL, "Human protein" = RED), name = NULL) +
  scale_fill_manual(values = c("Mouse protein" = FOCAL, "Human protein" = RED), name = NULL) +
  scale_y_continuous(labels = MINUS_LABELS) +
  labs(x = "m6A sites (human: brain peaks) in the gene", y = "Prediction with minus\nwithout sites (z)") +
  theme_fig() +
  theme(legend.position = "bottom")

# c: held-out R2 lost when a feature group is removed (x 0.001), one panel per target (point; grey bar, 95% of the
# permuted copies of the tested group "m6A sites and binding"; its P above the panel)
d <- read_panel("sfig3c")
M6A <- c("m6A sites and binding", "m6A sites", "m6A reader binding")
others <- d[!block %in% M6A, .(m = mean(delta_r2)), by = block][order(-m), block]
GROUPS <- c(M6A, others)
strips <- d[!is.na(p), setNames(paste0(target, "\n", fmt_p(p, "m6A P")), target)]
d[, `:=`(block = factor(block, GROUPS), delta_r2 = 1000 * delta_r2, null_low = 1000 * null_low,
         null_high = 1000 * null_high, target = factor(target, TARGETS))]
c_ <- ggplot(d, aes(delta_r2, block)) +
  ZERO_X +
  geom_linerange(data = d[!is.na(null_low)], aes(xmin = null_low, xmax = null_high), orientation = "y",
                 linewidth = 4, colour = NS) +
  geom_point(aes(fill = I(SPECIES(target))), shape = 21, colour = EDGE, size = 3.6, stroke = 0.5) +
  facet_grid(. ~ target, scales = "free_x", labeller = labeller(target = strips)) +
  scale_y_discrete(limits = rev(GROUPS), expand = expansion(add = 0.6)) +
  scale_x_continuous(labels = MINUS_LABELS, n.breaks = 4, expand = expansion(mult = 0.1)) +
  labs(x = "Held-out R² lost when the group is removed (× 0.001)", y = NULL) +
  theme_rows() +
  BOXED +
  theme(strip.text = element_text(lineheight = 0.95))

h_c <- 24 + 6.5 * length(GROUPS)
fig <- panel_grid(list(a = a, b = b, c = c_), design = "AB\nCC", widths = c(100, 83), heights = c(75, h_c))
save_fig(fig, "sfig3", 2 * (75 + h_c))
