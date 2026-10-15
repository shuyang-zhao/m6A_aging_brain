# Additional file 1 Fig. S3: the models of gene-level ageing (script 50): design, accuracy and the deletion effect by
# number of sites
# reads source_data/sfig3b, c; writes output/figures/sfig3.png and .pdf
source("figures/style.R")

TRANSFER <- "Human protein (mouse model)"

# a: schematic of the model: seven feature groups per gene are the input of TabPFN, which predicts the combined age z
# of each gene; the genes of each of five chromosome groups are held out in turn
GROUPS_A <- data.table(
  y = 7:1, label = c("m6A sites", "m6A reader binding", "Transcript structure", "Expression and turnover",
                     "Cell-type origin", "Functional programmes", "Sequence representation"),
  fill = c(FOCAL_LIGHT, RED_LIGHT, "#E3E3E3", "#E3E3E3", "#E3E3E3", "#E3E3E3", "#E3E3E3"),
  edge = c(FOCAL, RED, GREY, GREY, GREY, GREY, GREY)
)
folds <- data.table(x = 75 + 5 * (0:4), held = c(FALSE, FALSE, TRUE, FALSE, FALSE))
a <- ggplot() +
  geom_tile(data = GROUPS_A, aes(21, y * 11 + 4, fill = I(fill), colour = I(edge)), width = 40, height = 8.5,
            linewidth = 0.7) +
  geom_text(data = GROUPS_A, aes(21, y * 11 + 4, label = label), size = PT, colour = INK) +
  geom_segment(aes(x = 42.5, xend = 47.5, y = 48, yend = 48), colour = EDGE, linewidth = 0.8,
               arrow = arrow(length = unit(2.2, "mm"), type = "closed")) +
  annotate("tile", x = 57.5, y = 48, width = 19, height = 16, fill = CONTRAST_LIGHT, colour = CONTRAST,
           linewidth = 0.7) +
  annotate("text", x = 57.5, y = 48, label = "TabPFN", size = PT, fontface = "bold", colour = INK) +
  geom_segment(aes(x = 67.5, xend = 71, y = 48, yend = 48), colour = EDGE, linewidth = 0.8,
               arrow = arrow(length = unit(2.2, "mm"), type = "closed")) +
  geom_tile(data = folds, aes(x, 48, fill = I(fifelse(held, BLUE, BLUE_LIGHT))), width = 4, height = 16,
            colour = "white") +
  annotate("text", x = 85, y = 35, label = "Age z of each gene;\nchromosome groups\nheld out in turn",
           size = PT, colour = INK, lineheight = 0.95, vjust = 1) +
  coord_cartesian(xlim = c(0, 100), ylim = c(5, 86.5), expand = FALSE) +
  theme_void(base_size = TEXT, base_family = FONT) +
  theme(plot.margin = margin(4.4, 5, 3, 2, "mm"))

# b: held-out Spearman rho of each model (point, 95% CI); last row, the mouse protein model applied to human proteins
b <- read_panel("sfig3b")
b[, `:=`(target = factor(target, rev(c(TARGETS, TRANSFER))),
         model = factor(c(tabpfn = "TabPFN", lightgbm = "LightGBM", elastic_net = "Elastic net",
                          shuffled_target = "Shuffled target")[model],
                        c("TabPFN", "LightGBM", "Elastic net", "Shuffled target")))]
b <- ggplot(b, aes(rho, target)) +
  ZERO_X +
  geom_linerange(aes(xmin = rho_ci_low, xmax = rho_ci_high), colour = CI, linewidth = 0.6) +
  geom_point(aes(colour = model, fill = model), shape = 21, size = 3.8, stroke = 0.5) +
  scale_colour_manual(values = c(TabPFN = EDGE, LightGBM = EDGE, "Elastic net" = EDGE, "Shuffled target" = GREY),
                      name = NULL) +
  scale_fill_manual(values = c(TabPFN = FOCAL, LightGBM = CONTRAST, "Elastic net" = BLUE, "Shuffled target" = "white"),
                    name = NULL) +
  scale_x_continuous(labels = MINUS_LABELS) +
  labs(x = "Held-out Spearman ρ", y = NULL) +
  theme_rows() +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))

# c: prediction with minus without the sites of a gene (negative: the sites make the predicted decline stronger), by
# its number of m6A sites (mean, 95% CI)
cc <- read_panel("sfig3c")
cc[, `:=`(sites = factor(sites, c("1", "2-3", "4-7", ">=8"), c("1", "2–3", "4–7", "≥8")),
         target = factor(target, c("Mouse protein", "Human protein")))]
cc <- ggplot(cc, aes(sites, mean, group = target)) +
  ZERO_Y +
  geom_linerange(aes(ymin = ci_low, ymax = ci_high), colour = CI, linewidth = 0.6) +
  geom_line(aes(colour = target), linewidth = 1) +
  geom_point(aes(fill = target), colour = EDGE, shape = 21, size = 3.6, stroke = 0.5) +
  scale_colour_manual(values = c("Mouse protein" = FOCAL, "Human protein" = CONTRAST), name = NULL) +
  scale_fill_manual(values = c("Mouse protein" = FOCAL, "Human protein" = CONTRAST), name = NULL) +
  scale_y_continuous(labels = MINUS_LABELS) +
  labs(x = "m6A sites (human: brain peaks) in the gene", y = "Prediction with minus\nwithout sites (z)") +
  theme_fig() +
  theme(legend.position = "bottom")

fig <- panel_grid(list(a = a, b = b, c = cc), design = "AB\nCB", widths = c(80, 103), heights = c(62, 62))
save_fig(fig, "sfig3", 248)
