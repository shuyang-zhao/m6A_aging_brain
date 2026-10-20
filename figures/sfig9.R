# Additional file 1 Fig. S9: IGF2BP2 falls in the striatum but its targets do not follow (cell-type specificity, neural
# targets, coupling, protein score, perturbations and the neuron-versus-glia recheck)
# reads source_data/sfig9a-f; writes output/figures/sfig9.png and .pdf
source("figures/style.R")

COMP <- c(unbound_m6A = "vs unbound m6A genes", bound_unmarked = "vs bound unmarked genes")
COMP_COL <- c("vs unbound m6A genes" = CONTRAST, "vs bound unmarked genes" = BLUE)
MRNA_COL <- "#4D4D4D"

# a: IGF2BP2 mRNA and the target scores of its primary tiers, difference in age slope (% per decade, 95% CI) between
# the striatal and the other regions of GTEx and between the expressing and non-expressing populations of the
# single-nucleus cohorts
PANELS_A <- c(GTEx = "GTEx: striatum minus\nother regions",
              "Single nucleus" = "Single nucleus: expressing minus\nnon-expressing populations")
d <- read_panel("sfig9a")[, panel := factor(PANELS_A[panel], PANELS_A)]
TIERS_A <- c(gene = "IGF2BP2 mRNA", T1 = "T1 neural binding", T2 = "T2 m6A-dependent", T3 = "T3 falls on knockdown")
COMP_A <- c("-" = "mRNA", unbound_m6A = "vs unbound m6A genes", bound_unmarked = "vs bound unmarked genes")
d[, `:=`(row = factor(TIERS_A[tier], rev(TIERS_A)), comp = factor(COMP_A[comparison], COMP_A))]
dodge <- position_dodge(0.6)
a <- ggplot(d, aes(estimate_pct, row, fill = comp, shape = comp)) +
  ZERO_X +
  geom_linerange(aes(xmin = ci_low_pct, xmax = ci_high_pct), colour = CI, linewidth = 0.6, position = dodge,
                 show.legend = FALSE) +
  geom_point(size = 3.2, stroke = 0.5, colour = EDGE, position = dodge) +
  scale_fill_manual(values = c(setNames(MRNA_COL, COMP_A[1]), COMP_COL), name = NULL) +
  scale_shape_manual(values = setNames(c(23, 21, 24), COMP_A), name = NULL) +
  facet_wrap(~panel) +
  scale_x_continuous(labels = MINUS_LABELS) +
  labs(x = "Difference in age slope (% per decade)", y = NULL) +
  theme_rows() +
  BOXED +
  theme(legend.position = "bottom", strip.text = element_text(lineheight = 0.95))

# b: genes bound by IGF2BP2 with m6A in each combination of the four neural sources (UpSet); T1, at least two sources
SOURCES <- c(conway2016_H9_eCLIP = "H9 stem cells eCLIP", degrauwe2016_MGG8_PARCLIP = "MGG8 glioma PAR-CLIP",
             han2022_A172_eCLIP = "A172 glioma eCLIP", preitner2016_E14brain_CLIP = "E14 brain CLIP")
u <- read_panel("sfig9b")[order(-genes)]
u[, `:=`(combo = factor(seq_len(.N)), t1 = rowSums(.SD) >= 2), .SDcols = names(SOURCES)]
dots <- melt(u, id.vars = c("combo", "genes", "t1"), measure.vars = names(SOURCES), variable.name = "source",
             value.name = "member")
dots[, source := factor(SOURCES[as.character(source)], rev(SOURCES))]
bars <- ggplot(u, aes(combo, genes, fill = t1)) +
  geom_col(width = 0.65, colour = EDGE, linewidth = 0.5) +
  geom_text(data = u[genes >= 30], aes(label = format(genes, big.mark = ",")), angle = 90, hjust = -0.15, size = PT,
            colour = INK) +
  scale_fill_manual(values = c("TRUE" = FOCAL, "FALSE" = GREY),
                    labels = c("TRUE" = "T1 (two or more sources)", "FALSE" = "One source"),
                    breaks = c("TRUE", "FALSE"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3)), labels = scales::label_comma()) +
  labs(x = NULL, y = "IGF2BP2 target\ngenes with m6A") +
  theme_fig() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(4.4, 5, 0, 2, "mm"),
        legend.position = "inside", legend.position.inside = c(0.98, 0.98), legend.justification = c(1, 1))
matrix <- ggplot(dots, aes(combo, source)) +
  geom_line(data = dots[member == TRUE], aes(group = combo), colour = EDGE, linewidth = 0.8) +
  geom_point(aes(colour = I(fifelse(member, fifelse(t1, FOCAL, GREY), "#E3E3E3"))), size = 3) +
  labs(x = NULL, y = NULL) +
  theme_fig() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
        plot.margin = margin(0, 5, 3, 2, "mm"))
b <- (bars / matrix + plot_layout(heights = c(1.5, 1))) + FLUSH

# c: coupling of IGF2BP2 mRNA with the T1 target score within striatal single-nucleus populations (Stability arm);
# bars, edge b with 95% CI; numbers, percentile of |b| among the control RBPs
d <- read_panel("sfig9c")
d[, label := paste(dataset, sub("_", " ", population))]
d[, label := factor(label, rev(unique(label)))]
cc <- ggplot(d, aes(b_s, label)) +
  ZERO_X +
  geom_col(width = 0.65, fill = RED, colour = EDGE, linewidth = 0.5) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), linewidth = 0.6, colour = EDGE) +
  geom_text(aes(x = fifelse(b_s < 0, ci_low, ci_high), label = sprintf("%.0f", pctile_abs),
                hjust = fifelse(b_s < 0, 1.3, -0.3)), size = PT, colour = INK) +
  scale_x_continuous(labels = MINUS_LABELS, breaks = scales::breaks_pretty(3), expand = expansion(mult = 0.25)) +
  labs(x = "Coupling edge (b)", y = NULL) +
  theme_fig() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

# d: age effect on the protein score of the primary targets (mouse, late versus early; human, per decade), pooled
# over units (95% CI)
p <- read_panel("sfig9d")
p[, `:=`(comparison = factor(COMP[comparison], COMP), species = factor(SPECIES_NAME[species], SPECIES_NAME))]
dodge <- position_dodge(0.6)
dd <- ggplot(p, aes(tier, estimate_pct, fill = comparison)) +
  ZERO_Y +
  geom_linerange(aes(ymin = ci_low_pct, ymax = ci_high_pct), colour = CI, linewidth = 0.6, position = dodge) +
  geom_point(shape = 21, size = 4, stroke = 0.5, colour = EDGE, position = dodge) +
  facet_wrap(~species, scales = "free_y") +
  scale_fill_manual(values = COMP_COL, name = NULL) +
  scale_y_continuous(labels = MINUS_LABELS) +
  labs(x = "Target tier", y = "Age effect on protein score (%)") +
  theme_fig() +
  BOXED +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 1))

# e: change of the target score after IGF2BP2 loss in every perturbation dataset (T1 to T3, both comparisons; bold,
# q < 0.05)
LABELS <- c(
  degrauwe2016_IMP2s_KD = "Degrauwe 2016 KD", degrauwe2016_brain_IMP2_KO = "Degrauwe 2016 KO",
  encode_hepg2_shRNA = "ENCODE HepG2", encode_k562_shRNA = "ENCODE K562",
  gse146726_J82_siIGF2BP2 = "Okholm 2020 J82", gse146726_UMUC3_siIGF2BP2 = "Okholm 2020 UM-UC-3",
  gse292558_bmm_ko = "Schymik 2026 KO", huang2018_hepg2_shIGF2BP2 = "Huang 2018 HepG2",
  regue2019_muscle_ko = "Regué 2019 KO"
)
ROWS_LOSS <- paste(rep(c("T1", "T2", "T3"), each = 2), c("vs unbound m6A", "vs bound unmarked"))
x <- read_panel("sfig9e")
x[, `:=`(row = factor(paste(tier, c(unbound_m6A = "vs unbound m6A", bound_unmarked = "vs bound unmarked")[comparison]),
                      rev(ROWS_LOSS)),
         label = factor(LABELS[dataset], LABELS))]
e <- ggplot(x, aes(label, row, fill = effect_pct)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = cell_number(effect_pct), fontface = fifelse(q < 0.05, "bold", "plain")), size = PT,
            colour = INK) +
  scale_fill_div(20, "Change of target score after IGF2BP2 loss (%)") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_heat() +
  theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5), legend.position = "bottom")

# f: neuron-versus-glia recheck, percentile of the IGF2BP2 target sets among 39 control RBPs (all genes or load-matched
# genes as the comparison); dashed lines, 2.5th and 97.5th percentile
n <- read_panel("sfig9f")
n[, `:=`(set = sub("bound-unmarked", "bound unmarked", sub("IGF2BP2 (primary )?", "", set)),
         comparison = factor(c(all = "All genes", matched = "Load-matched genes")[comparison],
                             c("All genes", "Load-matched genes")))]
n[, set := sub("^\\(CLIP within 50 nt\\)", "CLIP within 50 nt", set)]
n[, set := factor(set, rev(unique(set)))]
f <- ggplot(n, aes(pctile_z, set)) +
  geom_vline(xintercept = c(2.5, 97.5), linewidth = 0.5, colour = EDGE, linetype = "22") +
  geom_line(aes(group = set), colour = EDGE, linewidth = 1.2) +
  geom_point(aes(fill = comparison), shape = 21, size = 4, stroke = 0.5, colour = EDGE) +
  scale_fill_manual(values = c("All genes" = GREY, "Load-matched genes" = BLUE), name = NULL) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  labs(x = "Percentile among\n39 control RBPs\n(neurons versus glia)", y = NULL) +
  theme_rows() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 1))

fig <- panel_grid(list(a = a, b = b, c = cc, d = dd, e = e, f = f),
                  design = c(area(1, 1, 1, 240), area(1, 241, 1, 366), area(2, 1, 2, 184), area(2, 185, 2, 366),
                             area(3, 1, 3, 265), area(3, 266, 3, 366)),
                  heights = c(52, 52, 62))
save_fig(fig, "sfig9", 2 * 166)
