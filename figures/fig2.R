# Figure 2: gene-level model of the age effect from transcript features
# reads source_data/fig2a-f (script 50); writes output/figures/fig2
source("figures/style.R")

performance <- read_panel("fig2b")

# a: held-out prediction against the observed age z of every gene, mouse and human protein
a <- read_panel("fig2a")
scatter <- function(name, hue, ylab) {
  fit <- performance[target == name & model == "TabPFN"]
  d <- by_density(a[target == name], "observed", "predicted")
  ggplot(d, aes(observed, predicted, colour = .density)) +
    ref_line(slope = 1) +
    rasterise_pts(geom_point(size = 0.35, stroke = 0)) +
    note_corner(name) +
    note_corner(lab_stats(fit$rho, stat = "rho"), vjust = 2.6, parse = TRUE) +
    note_corner(lab_stats(n = fit$n), vjust = 3.95, parse = TRUE) +
    scale_density(hue) +
    scale_x_continuous(labels = MINUS_LABELS) +
    scale_y_continuous(labels = MINUS_LABELS) +
    coord_cartesian(xlim = quantile(d$observed, c(0.0005, 0.9995)), ylim = range(d$predicted) * 1.05) +
    labs(x = "Observed age z", y = ylab) +
    theme_pub()
}
a <- scatter("Mouse protein", MOUSE, "Predicted age z") + scatter("Human protein", HUMAN, NULL)

# b: held-out Spearman rho of the four models on the five targets, a bar in each cell
b <- copy(performance)
HEADERS <- c("Mouse\nprotein", "Mouse mRNA-\nadjusted protein", "Human\nprotein", "Human\nmRNA", "Human neuronal\nmRNA")
MODELS <- c("TabPFN", "LightGBM", "Elastic net", "Shuffled target")
b[, `:=`(column = match(target, TARGETS), row = match(model, rev(MODELS)))]
b[, `:=`(left = column - 0.45, right = column - 0.45 + pmax(rho, 0) * 1.05, fill = fifelse(model == "TabPFN", M6BRAIN, GREY))]
head <- b[model == "TabPFN", .(column, text = paste0(HEADERS[column], "\nn = ", thousands(n)))]
b <- ggplot(b) +
  geom_rect(aes(xmin = left, xmax = right, ymin = row - 0.3, ymax = row + 0.3, fill = I(fill))) +
  geom_segment(aes(x = left, xend = left, y = row - 0.42, yend = row + 0.42), linewidth = lw(0.5)) +
  geom_text(aes(right + 0.04, row, label = number(rho), fontface = fifelse(model == "TabPFN", "bold", "plain")),
            hjust = 0, size = pt_size(6)) +
  geom_text(data = head, aes(column - 0.45, 4.7, label = text), hjust = 0, vjust = 0, lineheight = 0.95,
            size = pt_size(6)) +
  scale_x_continuous(limits = c(0.55, 5.65), expand = c(0, 0)) +
  scale_y_continuous(breaks = 1:4, labels = rev(MODELS), limits = c(0.5, 6.1), expand = c(0, 0)) +
  labs(x = "Held-out Spearman ρ") +
  theme_tile() +
  theme(axis.title.x = element_text(size = 7, margin = margin(t = 1, unit = "mm")), axis.text.x = element_blank())

# c: held-out R2 lost when a feature group is removed, mouse left and human right, group names in the gutter; band, the
# m6A groups; P, permutation of the tested group against 100 models
c_ <- read_panel("fig2c")
GUTTER <- 0.04
c_ <- dcast(c_, block ~ target, value.var = "delta_r2")[, .(block, mouse = `Mouse protein`, human = `Human protein`)]
c_[, row := rank((mouse + human) / 2)]
tested <- read_panel("fig2c")[!is.na(p)]
band <- c_[grepl("m6A", block), .(row)]
tested_row <- c_[block == tested$block[1], row]
p_of <- tested[, setNames(p, target)]
c_ <- ggplot(c_) +
  geom_rect(data = band, aes(ymin = row - 0.5, ymax = row + 0.5), xmin = -Inf, xmax = Inf, fill = tint(M6BRAIN, 0.88)) +
  geom_rect(aes(xmin = -GUTTER - pmax(mouse, 0), xmax = -GUTTER, ymin = row - 0.34, ymax = row + 0.34), fill = MOUSE) +
  geom_rect(aes(xmin = GUTTER, xmax = GUTTER + pmax(human, 0), ymin = row - 0.34, ymax = row + 0.34), fill = HUMAN) +
  geom_text(aes(0, row, label = block), size = pt_size(6)) +
  geom_text(aes(-GUTTER - pmax(mouse, 0) - 0.002, row, label = number(mouse, 3)), hjust = 1, size = pt_size(6)) +
  geom_text(aes(GUTTER + pmax(human, 0) + 0.002, row, label = number(human, 3)), hjust = 0, size = pt_size(6)) +
  note(c(-GUTTER, GUTTER), 8.95, c("Mouse protein", "Human protein"), hjust = c(1, 0)) +
  note(-GUTTER - 0.022, tested_row, paste0("Mouse P = ", signif(p_of["Mouse protein"], 2), "\nHuman P = ", signif(p_of["Human protein"], 2)),
       hjust = 1, lineheight = 0.95) +
  scale_x_continuous(breaks = c(-GUTTER - c(0.1, 0.05, 0), GUTTER + c(0, 0.05)), labels = c("0.1", "0.05", "0", "0", "0.05"),
                     expand = expansion(mult = c(0.06, 0.1))) +
  scale_y_continuous(limits = c(0.5, 9.3), expand = c(0, 0)) +
  labs(x = expression("Held-out" ~ italic(R)^2 ~ "lost on removing the group"), y = NULL) +
  theme_pub() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())

# d: in silico deletion of the sites of each region: mean change of the prediction of marked genes (bar), the 95% range
# of 100 models with the m6A features permuted (grey band) and the permutation P, along the transcript
d <- read_panel("fig2d")
REGIONS <- c("5' UTR", "CDS", "3' UTR")
regions <- d[region %in% REGIONS]
regions[, x := match(region, REGIONS) - 0.5]
all_sites <- d[region == "All sites"]
deletion <- function(name, hue) {
  r <- regions[target == name]
  o <- all_sites[target == name]
  floor <- min(r$effect) * 1.5
  ggplot(r, aes(x, effect)) +
    geom_rect(aes(xmin = x - 0.36, xmax = x + 0.36, ymin = null_low, ymax = null_high), fill = LIGHT) +
    geom_col(fill = hue, width = 0.46) +
    geom_hline(yintercept = 0, linewidth = lw(0.5)) +
    geom_text(aes(y = effect, label = fmt_p(signif(p, 2))), parse = TRUE, vjust = 1.35, size = pt_size(6)) +
    note(0.02, floor * 0.98, sprintf("%s, all sites %s", sub(" protein", "", name), number(o$effect)), hjust = 0, vjust = 0) +
    scale_x_continuous(limits = c(0, 3), expand = c(0, 0)) +
    scale_y_continuous(labels = MINUS_LABELS, breaks = scales::breaks_pretty(3), expand = expansion(mult = c(0, 0.04))) +
    coord_cartesian(ylim = c(floor, NA)) +
    labs(x = NULL, y = "Age z change") +
    theme_pub() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
}
site_flags <- data.frame(x = c(0.35, 1.15, 1.6, 1.85, 2.2, 2.45, 2.6, 2.8))
molecule <- tx_cartoon(c("5′ UTR", "CDS", "3′ UTR"), c(FALSE, TRUE, FALSE), fill = tint(M6BRAIN, 0.7)) +
  geom_segment(data = site_flags, aes(x, 0.2, xend = x, yend = 0.75), linewidth = lw(0.5), colour = M6BRAIN) +
  geom_point(data = site_flags, aes(x, 0.75), size = 0.7, colour = M6BRAIN) +
  scale_y_continuous(limits = c(-1, 0.95))
d <- deletion("Mouse protein", MOUSE) / deletion("Human protein", HUMAN) / molecule +
  plot_layout(heights = c(1, 1, 0.42))

# e: deletion effect of every gene with sites, ranked; the genes at both ends named
e <- read_panel("fig2e")
ranked <- function(name, hue, low, high, ylab, names_low, names_high) {
  g <- e[target == name][order(effect_all)]
  g[, `:=`(rank = .I, tail = .I <= low | .I > .N - high)]
  named <- g[tail == TRUE]
  ggplot(g, aes(rank, effect_all)) +
    ref_line(yintercept = 0) +
    rasterise_pts(geom_point(aes(colour = tail, size = tail), stroke = 0)) +
    label_column(named[effect_all < 0], "rank", "effect_all", "symbol", nrow(g) * 0.14, names_low[1], names_low[2],
                 fontface = "italic") +
    label_column(named[effect_all > 0], "rank", "effect_all", "symbol", nrow(g) * 0.84, names_high[1], names_high[2],
                 hjust = 1, fontface = "italic") +
    scale_colour_manual(values = c("TRUE" = hue, "FALSE" = GREY), guide = "none") +
    scale_size_manual(values = c("TRUE" = 1.1, "FALSE" = 0.4), guide = "none") +
    scale_x_continuous(labels = thousands, breaks = scales::breaks_pretty(3), expand = expansion(mult = 0.02)) +
    scale_y_continuous(labels = MINUS_LABELS) +
    labs(x = paste("Rank,", tolower(name)), y = ylab) +
    theme_pub()
}
e <- ranked("Mouse protein", MOUSE, 7, 3, "Prediction with minus without sites", c(-0.8, 0.34), c(0.8, 0.2)) +
  ranked("Human protein", HUMAN, 5, 3, NULL, c(-0.1, 0.045), c(0.24, 0.05))

# f: the mouse protein model applied to human orthologs: observed human age z by decile of the prediction (mean, SE)
f <- read_panel("fig2f_deciles")
stats <- read_panel("fig2f")[model == "Mouse protein model"]
TRANSFER <- c("Human protein" = HUMAN, "Human mRNA" = shade(HUMAN, 0.45))
f[, label := paste0(target, "\nρ = ", number(stats[match(f$target, target), rho]))]
ends <- f[decile == 10][, y := observed + fifelse(target == "Human protein", -0.3, 0.3)]
f <- ggplot(f, aes(decile, observed, colour = target, fill = target)) +
  ref_line(yintercept = 0) +
  geom_ribbon(aes(ymin = observed - se, ymax = observed + se), colour = NA, alpha = 0.2) +
  geom_line(linewidth = lw(0.9)) +
  geom_point(size = 0.9) +
  geom_text(data = ends, aes(x = 10.3, y = y, label = label), hjust = 0, lineheight = 0.95, size = pt_size(6), colour = "black") +
  scale_colour_manual(values = TRANSFER, guide = "none") +
  scale_fill_manual(values = TRANSFER, guide = "none") +
  scale_x_continuous(breaks = c(1, 5, 10), expand = expansion(add = c(0.4, 0))) +
  scale_y_continuous(labels = MINUS_LABELS) +
  coord_cartesian(xlim = c(1, 14), clip = "off") +
  labs(x = "Decile of the mouse model prediction", y = "Observed human age z") +
  theme_pub()

fig <- fig_compose(a, b, c_, d, e, f, design = "AAAAAAAABBBBBBBBBB\nCCCCCCCCCCDDDDDDDD\nEEEEEEEEEEEEFFFFFF",
                   heights = c(52, 58, 50))
save_fig(fig, "fig2", 168)
