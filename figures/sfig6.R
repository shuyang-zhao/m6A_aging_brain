# Additional file 1 Fig. S6: IGF2BP2 falls in the striatum but its targets do not follow (region and cell-type
# specificity, neural targets, coupling, mediation, protein score, perturbations and the neuron-versus-glia recheck)
# reads source_data/sfig6a-g; writes output/figures/sfig6.png and .pdf
source("figures/style.R")
suppressPackageStartupMessages(library(ggbeeswarm))

COMP <- c(unbound_m6A = "vs unbound m6A genes", bound_unmarked = "vs bound unmarked genes")
COMP_COL <- c("vs unbound m6A genes" = FOCAL, "vs bound unmarked genes" = RED)

# a: IGF2BP2 mRNA and the target scores of its primary tiers, region x age in GTEx (393 donors; striatum minus the other
# regions) and expressing minus non-expressing populations of the single-nucleus cohorts; % per decade, 95% CI
d <- read_panel("sfig6a")
d[, panel := c(GTEx = "GTEx: striatum minus\nother regions",
               "Single nucleus" = "Single nucleus: expressing minus\nnon-expressing populations")[panel]]
TIERS_A <- c(gene = "IGF2BP2 mRNA", T1 = "T1 neural binding", T2 = "T2 m6A-dependent", T3 = "T3 falls on knockdown")
COMP_A <- c("-" = "mRNA", unbound_m6A = "vs unbound m6A genes", bound_unmarked = "vs bound unmarked genes")
d[, `:=`(row = factor(TIERS_A[tier], rev(TIERS_A)), comp = factor(COMP_A[comparison], COMP_A))]
dodge <- position_dodge(0.6)
a <- ggplot(d, aes(estimate_pct, row, fill = comp, shape = comp)) +
  ZERO_X +
  geom_linerange(aes(xmin = ci_low_pct, xmax = ci_high_pct), colour = CI, linewidth = 0.6, position = dodge,
                 show.legend = FALSE) +
  geom_point(size = 3.2, stroke = 0.5, colour = EDGE, position = dodge) +
  facet_grid(. ~ panel, scales = "free_x") +
  scale_fill_manual(values = c(setNames(CONTRAST, COMP_A[1]), COMP_COL), name = NULL) +
  scale_shape_manual(values = setNames(c(23, 21, 24), COMP_A), name = NULL) +
  scale_x_continuous(labels = MINUS_LABELS) +
  labs(x = "Difference in age slope (% per decade)", y = NULL) +
  theme_rows() +
  BOXED +
  theme(legend.position = "bottom", strip.text = element_text(lineheight = 0.95))

# b: genes bound by IGF2BP2 with m6A in each combination of the four neural sources (UpSet); T1, at least two sources
SOURCES <- c(conway2016_H9_eCLIP = "H9 stem cells eCLIP", degrauwe2016_MGG8_PARCLIP = "MGG8 glioma PAR-CLIP",
             han2022_A172_eCLIP = "A172 glioma eCLIP", preitner2016_E14brain_CLIP = "E14 brain CLIP")
u <- read_panel("sfig6b")[order(-genes)]
u[, `:=`(combo = factor(seq_len(.N)), t1 = rowSums(.SD) >= 2), .SDcols = names(SOURCES)]
b <- upset(
  u, SOURCES,
  aes(fill = t1),
  geom_col(width = 0.65, colour = EDGE, linewidth = 0.5),
  geom_text(data = u[genes >= 30], aes(label = format(genes, big.mark = ",")), angle = 90, hjust = -0.15, size = PT,
            colour = INK),
  scale_fill_manual(values = c("TRUE" = FOCAL, "FALSE" = GREY),
                    labels = c("TRUE" = "T1 (two or more sources)", "FALSE" = "One source"),
                    breaks = c("TRUE", "FALSE"), name = NULL),
  scale_y_continuous(expand = expansion(mult = c(0, 0.3)), labels = scales::label_comma()),
  labs(y = "IGF2BP2 target genes with m6A"),
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.98), legend.justification = c(1, 1)),
  dot_colour = \(d) fifelse(d$t1, FOCAL, GREY), dot_size = 3, heights = c(1.5, 1)
)

# c: coupling of IGF2BP2 mRNA with the T1 target score within striatal single-nucleus populations (Stability arm);
# bars, edge b with 95% CI; numbers, percentile of |b| among the 90 control RBPs
d <- read_panel("sfig6c")
d[, label := paste(dataset, sub("_", " ", population))]
d[, label := factor(label, rev(unique(label)))]
cc <- ggplot(d, aes(b_s, label)) +
  ZERO_X +
  geom_col(width = 0.65, fill = CONTRAST, colour = EDGE, linewidth = 0.5) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), linewidth = 0.6, colour = EDGE) +
  geom_text(aes(x = fifelse(b_s < 0, ci_low, ci_high), label = sprintf("%.0f", pctile_abs),
                hjust = fifelse(b_s < 0, 1.3, -0.3)), size = PT, colour = INK) +
  scale_x_continuous(labels = MINUS_LABELS, breaks = scales::breaks_pretty(3), expand = expansion(mult = 0.25)) +
  labs(x = "Coupling edge (b)", y = NULL) +
  theme_fig() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank())

# d: fraction of the age slope of the target score taken over by IGF2BP2 mRNA; one dot per striatal population; bar,
# median; axis cut at ±150%
m <- read_panel("sfig6d")
m[, `:=`(comparison = factor(COMP[comparison], COMP), f = pmax(pmin(fraction_pct, 150), -150))]
dd <- ggplot(m, aes(tier, f, fill = comparison)) +
  ZERO_Y +
  stat_summary(aes(colour = comparison), fun = median, geom = "crossbar", width = 0.6, linewidth = 0.7,
               position = position_dodge(0.75), show.legend = FALSE) +
  geom_quasirandom(shape = 21, size = 2.8, stroke = 0.5, colour = EDGE, dodge.width = 0.75, width = 0.12) +
  scale_colour_manual(values = COMP_COL, guide = "none") +
  scale_fill_manual(values = COMP_COL, name = NULL) +
  scale_y_continuous(labels = MINUS_LABELS, limits = c(-150, 150)) +
  labs(x = "Target tier", y = "Mediated by IGF2BP2 mRNA (%)") +
  theme_fig() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 1))

# e: age effect on the protein score of the primary targets (mouse, late versus early; human, per decade), pooled
# over units (95% CI)
p <- read_panel("sfig6e")
p[, `:=`(comparison = factor(COMP[comparison], COMP),
         species = factor(c(mouse = "Mouse", human = "Human")[species], c("Mouse", "Human")))]
dodge <- position_dodge(0.6)
e <- ggplot(p, aes(tier, estimate_pct, fill = comparison)) +
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

# f: change of the target score after IGF2BP2 loss in every perturbation dataset (T1 to T3, both comparisons; bold,
# q < 0.05)
LABELS <- c(
  degrauwe2016_IMP2s_KD = "Degrauwe 2016 KD", degrauwe2016_brain_IMP2_KO = "Degrauwe 2016 KO",
  encode_hepg2_shRNA = "ENCODE HepG2", encode_k562_shRNA = "ENCODE K562",
  gse146726_J82_siIGF2BP2 = "Okholm 2020 J82", gse146726_UMUC3_siIGF2BP2 = "Okholm 2020 UM-UC-3",
  gse292558_bmm_ko = "Schymik 2026 KO", huang2018_hepg2_shIGF2BP2 = "Huang 2018 HepG2",
  regue2019_muscle_ko = "Regué 2019 KO"
)
ROWS_LOSS <- paste(rep(c("T1", "T2", "T3"), each = 2), c("vs unbound m6A", "vs bound unmarked"))
x <- read_panel("sfig6f")
x[, `:=`(row = factor(paste(tier, c(unbound_m6A = "vs unbound m6A", bound_unmarked = "vs bound unmarked")[comparison]),
                      rev(ROWS_LOSS)),
         label = factor(LABELS[dataset], LABELS))]
f <- ggplot(x, aes(label, row, fill = effect_pct)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = cell_number(effect_pct), fontface = fifelse(q < 0.05, "bold", "plain")), size = PT,
            colour = INK) +
  scale_fill_div(20, "Change of target score after IGF2BP2 loss (%)") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_heat() +
  theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5), legend.position = "bottom")

# g: neuron-versus-glia recheck, percentile of the IGF2BP2 target sets among 39 control RBPs (all genes or load-matched
# genes as the comparison); dashed lines, 2.5th and 97.5th percentile
n <- read_panel("sfig6g")
n[, `:=`(set = sub("IGF2BP2 (primary )?", "", set),
         comparison = factor(c(all = "All genes", matched = "Load-matched genes")[comparison],
                             c("All genes", "Load-matched genes")))]
n[, set := sub("^\\(CLIP within 50 nt\\)", "CLIP within 50 nt", set)]
n[, set := factor(set, rev(unique(set)))]
g <- ggplot(n, aes(pctile_z, set)) +
  geom_vline(xintercept = c(2.5, 97.5), linewidth = 0.5, colour = EDGE, linetype = "22") +
  geom_line(aes(group = set), colour = EDGE, linewidth = 1.2) +
  geom_point(aes(fill = comparison), shape = 21, size = 4, stroke = 0.5, colour = EDGE) +
  scale_fill_manual(values = c("All genes" = GREY, "Load-matched genes" = BLUE), name = NULL) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  labs(x = "Percentile among 39 control RBPs\n(neurons versus glia)", y = NULL) +
  theme_rows() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 1))

fig <- panel_grid(list(a = a, b = b, c = cc, d = dd, e = e, f = f, g = g), design = "AAAAA\nBBCCD\nEFFFG",
                  widths = c(44, 18, 40, 20, 62), heights = c(50, 58, 62))
save_fig(fig, "sfig6", 340)
