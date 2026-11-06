# Additional file 1 Fig. S3: model details, accuracy of every model on every target, the deletion effect by number of
# sites and the feature groups
# reads source_data/sfig3a-c; writes output/figures/sfig3.png and .pdf
source("figures/style.R")

MODELS <- c("TabPFN", "LightGBM", "Elastic net", "Shuffled target")
TRANSFER <- c("Human protein (mouse model)", "Human mRNA (mouse model)")

# a: held-out Spearman rho of each model on the five targets (95% bootstrap CI below); last rows, the mouse protein
# model applied to the human genes
a <- read_panel("sfig3a")
a[, `:=`(row = factor(target, rev(c(TARGETS, TRANSFER))), model = factor(model, MODELS),
         part = factor(fifelse(target %in% TRANSFER, "transfer", "own"), c("own", "transfer")),
         text = sprintf("%s\n(%s, %s)", number(rho), number(rho_ci_low), number(rho_ci_high)))]
a <- ggplot(a, aes(model, row, fill = rho)) +
  geom_tile(data = a[model != "Shuffled target"], colour = "white", linewidth = lw(0.6)) +
  geom_tile(data = a[model == "Shuffled target"], fill = LIGHT, colour = "white", linewidth = lw(0.6)) +
  geom_text(aes(label = text), size = pt_size(6), lineheight = 0.95) +
  facet_grid(part ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradientn(colours = c("#FFFFFF", tint(M6BRAIN, 0.7), tint(M6BRAIN, 0.35)), limits = c(0, 0.6),
                       oob = scales::squish, guide = "none") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Held-out Spearman ρ (95% CI)") +
  theme_tile() +
  theme(axis.title.x = element_text(), strip.text = element_blank(), panel.spacing.y = unit(2, "mm"))

# b: prediction with minus without the sites of a gene (negative: the sites make the predicted decline stronger) by its
# number of m6A sites (mean, 95% CI as the band)
b <- read_panel("sfig3b")
LEVELS <- c("1", "2-3", "4-7", "≥8")
b[, x := match(sites, LEVELS)]
b <- ggplot(b, aes(x, mean, colour = target, fill = target)) +
  ref_line(yintercept = 0) +
  geom_ribbon(aes(ymin = ci_low, ymax = ci_high), colour = NA, alpha = 0.18) +
  geom_line(linewidth = lw(0.9)) +
  geom_point(shape = 21, colour = EDGE, stroke = st(0.4), size = 1.5) +
  geom_text(data = b[x == 4], aes(label = target), colour = "black", hjust = 0, nudge_x = 0.12, size = pt_size(6)) +
  scale_colour_manual(values = c("Mouse protein" = MOUSE, "Human protein" = HUMAN), guide = "none") +
  scale_fill_manual(values = c("Mouse protein" = MOUSE, "Human protein" = HUMAN), guide = "none") +
  scale_x_continuous(breaks = 1:4, labels = c("1", "2–3", "4–7", "≥8"), expand = expansion(mult = c(0.04, 0.45))) +
  scale_y_continuous(labels = MINUS_LABELS) +
  labs(x = "m6A sites (human: brain peaks) in the gene", y = "Prediction with minus\nwithout sites (age z)") +
  theme_pub()

# c: held-out R2 lost when a feature group is removed, for the five targets (dot area and shade, x 0.001; number
# beside the dot); the m6A groups on a light band, and for the tested group the permutation P against 100 models with
# its features permuted across genes
d <- read_panel("sfig3c")
M6A <- c("m6A sites and binding", "m6A sites", "m6A reader binding")
GROUPS <- c(M6A, d[!block %in% M6A, .(m = mean(delta_r2)), by = block][order(-m), block])
d[, `:=`(row = factor(block, rev(GROUPS)), column = factor(target, TARGETS), value = 1000 * delta_r2)]
d[, `:=`(text = number(value, 1), line = as.numeric(row) + fifelse(is.na(p), 0, 0.2))]
band <- data.table(i = 1:3)
c_ <- ggplot(d, aes(column, row)) +
  geom_rect(data = band, aes(ymin = length(GROUPS) - i + 0.5, ymax = length(GROUPS) - i + 1.5), xmin = -Inf, xmax = Inf,
            fill = LIGHT, alpha = 0.6, inherit.aes = FALSE) +
  geom_point(aes(size = pmax(value, 0), fill = pmax(value, 0)), shape = 21, colour = EDGE, stroke = st(0.4)) +
  geom_text(aes(x = as.numeric(column) + 0.2, y = line, label = text), hjust = 0, size = pt_size(6)) +
  geom_text(data = d[!is.na(p)], aes(x = as.numeric(column) + 0.2, y = as.numeric(row) - 0.2, label = fmt_p(p)),
            hjust = 0, size = pt_size(6), parse = TRUE) +
  scale_size_area(max_size = 8, name = "Held-out R² lost (× 0.001)", breaks = c(2, 20, 100),
                  guide = guide_legend(nrow = 1, title.position = "left")) +
  scale_fill_gradientn(colours = c(tint(M6BRAIN, 0.8), tint(M6BRAIN, 0.4), M6BRAIN), trans = "sqrt", limits = c(0, 125),
                       name = "Held-out R² lost (× 0.001)", breaks = c(2, 20, 100),
                       guide = guide_legend(nrow = 1, title.position = "left")) +
  scale_x_discrete(position = "top", labels = scales::label_wrap(14), expand = expansion(add = c(0.5, 1))) +
  scale_y_discrete(labels = scales::label_wrap(24)) +
  theme_tile() +
  theme(panel.grid.major = element_line(linewidth = lw(0.25), colour = LIGHT), legend.position = "bottom",
        legend.margin = margin(t = 1, unit = "mm"))

fig <- fig_compose(a, b, c_, design = "AAAAAABBBB\nCCCCCCCCCC", heights = c(56, 74))
save_fig(fig, "sfig3", 134)
