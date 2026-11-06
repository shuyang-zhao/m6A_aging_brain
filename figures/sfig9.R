# Fig. S9: IGF2BP2 falls in the striatum but its targets do not follow (cell-type specificity, neural targets,
# coupling, protein score, perturbations and the neuron-versus-glia recheck)
# reads source_data/sfig9a-f; writes output/figures/sfig9
source("figures/style.R")

COMP <- c(unbound_m6A = "vs unbound m6A", bound_unmarked = "vs bound unmarked")
TIERS <- c(T1 = "T1 neural\nbinding", T2 = "T2 m6A-dependent", T3 = "T3 falls on\nknockdown")
COMP_ROWS <- function(d) d[, row := factor(COMP[comparison], rev(COMP))]
BAR <- 0.62
P <- function(p) patchwork::wrap_elements(full = p) # panels keep their own widths, so the grid columns stay equal
tier_strips <- theme(strip.placement = "outside", strip.text.x = element_text(size = 6), strip.text.y.left = element_text(angle = 0, hjust = 1, lineheight = 0.95),
                     panel.spacing.y = unit(1.2, "mm"), panel.spacing.x = unit(2, "mm"))

# a: IGF2BP2 mRNA (dark) and the target scores of its tiers (light), difference in age slope (% per decade, 95% CI) between
# the striatal and the other regions of GTEx and between the expressing and non-expressing populations of the
# single-nucleus cohorts
PANELS_A <- c(GTEx = "GTEx: striatum\nminus other regions",
              "Single nucleus" = "Single nucleus: expressing\nminus non-expressing")
a <- read_panel("sfig9a")
a[, `:=`(panel = factor(PANELS_A[panel], PANELS_A), tier = factor(c(gene = "IGF2BP2 mRNA", TIERS)[tier], c("IGF2BP2 mRNA", TIERS)),
         row = factor(fifelse(comparison == "-", " ", COMP[comparison]), rev(c(" ", COMP))),
         fill = fifelse(comparison == "-", HUMAN, tint(HUMAN, 0.55)))]
a <- ggplot(a, aes(estimate_pct, row)) +
  geom_vline(xintercept = 0, linewidth = lw(0.5)) +
  geom_col(aes(fill = I(fill)), width = BAR) +
  geom_errorbar(aes(xmin = ci_low_pct, xmax = ci_high_pct), width = 0, linewidth = lw(0.5), colour = EDGE) +
  facet_grid(tier ~ panel, scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(labels = MINUS_LABELS, breaks = seq(-12, 0, 4)) +
  coord_cartesian(xlim = c(-14, 3)) +
  labs(x = "Difference in age slope (% per decade)", y = NULL) +
  theme_box() +
  tier_strips +
  theme(axis.ticks.y = element_blank())

# b: genes bound by IGF2BP2 with m6A in each combination of the four neural sources (UpSet); T1, at least two sources
SOURCES <- c(conway2016_H9_eCLIP = "H9 stem cells eCLIP", degrauwe2016_MGG8_PARCLIP = "MGG8 glioma PAR-CLIP",
             han2022_A172_eCLIP = "A172 glioma eCLIP", preitner2016_E14brain_CLIP = "E14 brain CLIP")
u <- read_panel("sfig9b")[order(-genes)]
u[, `:=`(combo = seq_len(.N), t1 = rowSums(.SD) >= 2), .SDcols = names(SOURCES)]
dots <- melt(u, id.vars = c("combo", "genes", "t1"), measure.vars = names(SOURCES), variable.name = "source",
             value.name = "member")
dots[, `:=`(source = factor(SOURCES[as.character(source)], rev(SOURCES)), fill = fifelse(member, fifelse(t1, EDGE, GREY), LIGHT))]
u[, fill := fifelse(t1, EDGE, GREY)]
b_bars <- ggplot(u, aes(combo, genes)) +
  geom_col(aes(fill = I(fill)), width = 0.7) +
  geom_text(aes(label = thousands(genes)), angle = 90, hjust = -0.15, size = pt_size(6)) +
  key_rows(c("T1 (2 or more sources)", "One source"), c(EDGE, GREY), x = 4.6, y = 3000, step = 600, width = 0.7, square = TRUE) +
  scale_x_continuous(limits = c(0.4, nrow(u) + 0.6), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.4)), labels = scales::label_comma()) +
  labs(x = NULL, y = "IGF2BP2 target genes\nwith m6A") +
  theme_pub() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
b_matrix <- ggplot(dots, aes(combo, source)) +
  geom_line(data = dots[member == TRUE], aes(group = combo), colour = EDGE, linewidth = lw(0.9)) +
  geom_point(aes(colour = I(fill)), size = 1.6) +
  scale_x_continuous(limits = c(0.4, nrow(u) + 0.6), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_tile() +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_blank(), plot.margin = margin(0, 1.5, 0, 1.5, "mm"))
b <- b_bars / b_matrix + plot_layout(heights = c(1.55, 1))

# c: coupling of IGF2BP2 mRNA with the T1 target score within striatal single-nucleus populations (stability arm), edge b
# (bar, 95% CI) sorted; number, percentile of |b| among the control proteins
CLASS <- c(D1_MSN = "Neuron", D2_MSN = "Neuron", MSN_unassigned = "Neuron", OPC = "Oligodendrocyte")
c_ <- read_panel("sfig9c")
c_[, `:=`(label = paste(dataset, sub("_", " ", population)), fill = CELL[CLASS[population]])]
c_ <- c_[order(b_s)][, label := factor(label, rev(label))]
c_[, `:=`(end = fifelse(b_s < 0, ci_low, ci_high), side = sign(b_s))]
c_ <- ggplot(c_, aes(b_s, label)) +
  geom_vline(xintercept = 0, linewidth = lw(0.5)) +
  geom_col(aes(fill = I(fill)), width = BAR) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, linewidth = lw(0.5), colour = EDGE) +
  geom_text(aes(x = end, label = ordinal(pctile_abs), hjust = fifelse(side < 0, 1.25, -0.25)), size = pt_size(6)) +
  key_rows(c("MSN", "OPC"), CELL[c("Neuron", "Oligodendrocyte")], x = -0.33, y = 2.9, step = 0.75, width = 0.03, square = TRUE) +
  scale_x_continuous(labels = MINUS_LABELS, breaks = c(-0.2, 0, 0.2)) +
  coord_cartesian(xlim = c(-0.34, 0.26)) +
  labs(x = "Coupling edge (b)\nnumber, percentile of |b| among controls", y = NULL) +
  theme_rows()

# d: age effect on the protein score of the tier targets (mouse, late versus early; human, per decade), pooled over
# units (95% CI)
d <- read_panel("sfig9d")
d[, `:=`(species = factor(SPECIES_NAME[species], SPECIES_NAME), tier = factor(TIERS[tier], TIERS), fill = tint(SPECIES[SPECIES_NAME[species]], 0.55))]
COMP_ROWS(d)
d <- ggplot(d, aes(estimate_pct, row)) +
  geom_vline(xintercept = 0, linewidth = lw(0.5)) +
  geom_col(aes(fill = I(fill)), width = BAR) +
  geom_errorbar(aes(xmin = ci_low_pct, xmax = ci_high_pct), width = 0, linewidth = lw(0.5), colour = EDGE) +
  geom_text(data = data.table(species = factor("Human", SPECIES_NAME), tier = factor(TIERS["T2"], TIERS)),
            aes(0.15, 1.5, label = "not tested"), hjust = 0, size = pt_size(6), inherit.aes = FALSE) +
  facet_grid(tier ~ species, scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(labels = MINUS_LABELS) +
  labs(x = "Age effect on protein score (%)", y = NULL) +
  theme_box() +
  tier_strips +
  theme(axis.ticks.y = element_blank())

# e: change of the target score after IGF2BP2 loss in every perturbation dataset (T1 to T3, both comparisons; bold, q < 0.05)
LABELS <- c(
  degrauwe2016_IMP2s_KD = "Degrauwe 2016 KD", degrauwe2016_brain_IMP2_KO = "Degrauwe 2016 KO",
  encode_hepg2_shRNA = "ENCODE HepG2", encode_k562_shRNA = "ENCODE K562",
  gse146726_J82_siIGF2BP2 = "Okholm 2020 J82", gse146726_UMUC3_siIGF2BP2 = "Okholm 2020 UM-UC-3",
  gse292558_bmm_ko = "Schymik 2026 KO", huang2018_hepg2_shIGF2BP2 = "Huang 2018 HepG2",
  regue2019_muscle_ko = "Regué 2019 KO"
)
e <- read_panel("sfig9e")
e[, `:=`(label = factor(LABELS[dataset], rev(LABELS)), comp = factor(COMP[comparison], COMP),
         text = fifelse(abs(effect_pct) < 10, number(effect_pct, 1), number(effect_pct, 0)))]
e <- ggplot(e, aes(tier, label, fill = effect_pct)) +
  geom_tile(colour = "white", linewidth = lw(0.6)) +
  geom_text(aes(label = text, fontface = fifelse(q < 0.05, "bold", "plain")), size = pt_size(6)) +
  facet_grid(. ~ comp) +
  scale_fill_signed(20, "Change of target score after IGF2BP2 loss (%)", light = TRUE) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_tile() +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 6), panel.spacing.x = unit(1.5, "mm"),
        strip.text = element_text(size = 6), strip.placement = "outside")

# f: neuron-versus-glia recheck, percentile of the IGF2BP2 target sets among 39 control RBPs (all genes, grey, or
# load-matched genes, brick, as the comparison); dashed lines, 2.5th and 97.5th percentile
f <- read_panel("sfig9f")
f[, set := sub("bound-unmarked", "bound unmarked", sub("IGF2BP2 (primary )?", "", set))]
f[, set := sub("^\\(CLIP within 50 nt\\)", "CLIP within 50 nt", set)]
ROWS_F <- c("CLIP within 50 nt", paste(rep(c("T1", "T2", "T3"), each = 2), c("vs unbound m6A", "vs bound unmarked")))
f[, set := factor(set, rev(ROWS_F))]
f[, comparison := factor(c(all = "All genes", matched = "Load-matched genes")[comparison], c("All genes", "Load-matched genes"))]
f <- ggplot(f, aes(pctile_z, set)) +
  ref_line(xintercept = c(2.5, 97.5), colour = "black") +
  geom_line(aes(group = set), colour = EDGE, linewidth = lw(0.9)) +
  geom_point(aes(colour = comparison), size = 2) +
  scale_colour_manual(values = c("All genes" = GREY, "Load-matched genes" = HUMAN), name = NULL) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  labs(x = "Percentile among 39 control RBPs\n(neurons versus other nuclei)", y = NULL) +
  theme_rows() +
  theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, -1, 0))

fig <- fig_compose(P(a), P(b), P(c_), P(d), P(e), P(f), design = "AAAAAAABBBBB\nCCCCCDDDDDDD\nEEEEEEFFFFFF", heights = c(50, 46, 54))
save_fig(fig, "sfig9", 152)
