# Additional file 1 Fig. S8: FMRP family and PRRC2A, sign check against knockouts, the values of the network edges and
# edges shared across units
# reads source_data/sfig8a, sfig8b and sfig8c; writes output/figures/sfig8.png and .pdf
source("figures/style.R")

# a: change of the m6A-marked target score after loss of FMR1 or PRRC2A in the knockout data; the edge in the proteome
# predicts a fall (filled, the observed direction agrees with the prediction of the protein edge); * q < 0.05
s <- read_panel("sfig8a")
TISSUES <- c(cortex = "Cortex", hippocampus = "Hippocampus", cerebellum = "Cerebellum",
             "cortex Pvalb neurons" = "Cortex Pvalb neurons", "cortex Camk2a neurons" = "Cortex Camk2a neurons",
             "brain Nestin-Cre" = "Brain, Nestin-Cre", "brain Olig2-Cre" = "Brain, Olig2-Cre")
s[, label := paste0(TISSUES[tissue], ", ", readout, " (", toupper(sub("^.*_", "", perturbation)), ")")]
s[, `:=`(regulator = factor(paste(regulator, "loss"), c("FMR1 loss", "PRRC2A loss")),
         agree = factor(fifelse(agrees == 1, "Agrees", "Disagrees"), c("Agrees", "Disagrees")))]
s[, label := factor(label, rev(unique(label[order(regulator, -effect)])))]
a <- ggplot(s, aes(effect, label, fill = agree)) +
  ZERO_X +
  geom_col(width = 0.66, colour = EDGE, linewidth = 0.5) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), linewidth = 0.6, colour = EDGE) +
  geom_text(data = s[q < 0.05],
            aes(x = fifelse(effect < 0, ci_low, ci_high), label = "*", hjust = fifelse(effect < 0, 1.4, -0.4)),
            size = PT, vjust = 0.75, colour = INK) +
  facet_grid(regulator ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(Agrees = FOCAL, Disagrees = "white"), name = NULL) +
  scale_x_continuous(labels = MINUS_LABELS, expand = expansion(mult = 0.15)) +
  scale_y_discrete(labels = scales::label_wrap(36)) +
  labs(x = "Change of m6A-marked target score", y = NULL) +
  theme_rows() +
  BOXED +
  theme(legend.position = "bottom",
        strip.text.y = element_text(angle = 0, hjust = 0, margin = margin(0, 0, 0, 2, "mm")))

# b: the edges as values: median edge over the units of the primary source in units of the 95th percentile of the
# control RBPs; *, m6A-specific edge in at least one unit; **, also shared across units; the YTHDF load edge counts as
# specific above the 95th percentile of the control RBPs
m <- read_panel("sfig8b")
ROWS <- c("YTHDF1", "YTHDF2", "YTHDF3", "YTHDC1", "YTHDC2", "FMR1", "FXR1", "FXR2", "PRRC2A", "HNRNPC", "RBMX",
          "IGF2BP2", "YTHDF module", "FMRP family", "Writer complex", "Erasers",
          "Module, mouse", "Module, human", "YTHDF3, mouse", "YTHDF3, human")
m[, row := sub(" load,", ",", sub("^YTHDF module load", "Module", row))]
m[, `:=`(row = factor(row, rev(ROWS)), layer = factor(layer, LAYERS),
         group = factor(group, c("Readers", "Modules", "Writers, erasers", "YTHDF load edge")),
         mark = fcase(shared, "**", m6A_specific, "*", default = ""))]
edges <- ggplot(m, aes(layer, row, fill = value)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = mark), size = PT, colour = INK, vjust = 0.75) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y", labeller = label_wrap_gen(9)) +
  scale_fill_div(2, "Edge / control\n95th percentile") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  theme_heat() +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 0), legend.position = "bottom")

# c: posterior mean of the edge of every regulator in every Output unit after sharing across units (mashr); large dots,
# local false sign rate below 0.05
e <- read_panel("sfig8c")
REGS <- c("FMR1", "FMRP family", "FXR1", "FXR2", "PRRC2A", "HNRNPC", "RBMX", "YTHDC1", "YTHDF1", "YTHDF2", "YTHDF3",
          "YTHDF module", "Writer complex", "Erasers")
e[, `:=`(regulator = factor(regulator, REGS), species = factor(SPECIES_NAME[species], SPECIES_NAME))]
e[, unit := factor(unit, unique(unit[order(species, unit)]))]
LIM <- 0.1
mashr <- ggplot(e, aes(regulator, unit)) +
  geom_point(aes(fill = posterior_mean, size = lfsr < 0.05), shape = 21, colour = MATRIX_EDGE, stroke = 0.5) +
  facet_grid(species ~ ., scales = "free_y", space = "free_y") +
  scale_fill_div(LIM, "Posterior mean edge, protein layer") +
  scale_size_manual(values = c("FALSE" = 2.2, "TRUE" = 5.6),
                    labels = c("FALSE" = "lfsr ≥ 0.05", "TRUE" = "lfsr < 0.05"), name = NULL,
                    guide = guide_legend(order = 2, override.aes = list(fill = "#BDBDBD", colour = MATRIX_EDGE))) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev) +
  theme_matrix() +
  BOXED +
  theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5), legend.position = "bottom", legend.box = "vertical",
        legend.spacing.y = unit(1, "mm"), strip.text.y = element_text(angle = 270))

fig <- panel_grid(list(a = a, b = edges, c = mashr), design = "AC\nBC", widths = c(92, 91), heights = c(66, 80))
save_fig(fig, "sfig8", 2 * 146)
