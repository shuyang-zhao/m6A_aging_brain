# Fig. S8: FMRP family and PRRC2A, sign check against knockouts, the values of the network edges and edges shared
# across units
# reads source_data/sfig8a, sfig8b and sfig8c; writes output/figures/sfig8
source("figures/style.R")

# a: change of the m6A-marked target score after loss of FMR1 or PRRC2A in the knockout data (95% CI, exact q beside
# the bar); filled, the observed direction agrees with the fall predicted by the protein coupling
s <- read_panel("sfig8a")
TISSUES <- c(cortex = "Cortex", hippocampus = "Hippocampus", cerebellum = "Cerebellum",
             "cortex Pvalb neurons" = "Cortex Pvalb neurons", "cortex Camk2a neurons" = "Cortex Camk2a neurons",
             "brain Nestin-Cre" = "Brain, Nestin-Cre", "brain Olig2-Cre" = "Brain, Olig2-Cre")
s[, label := paste0(TISSUES[tissue], ", ", readout, " (", toupper(sub("^.*_", "", perturbation)), ")")]
s[, `:=`(regulator = factor(paste(regulator, "loss"), c("FMR1 loss", "PRRC2A loss")),
         agree = factor(fifelse(agrees == 1, "Agrees", "Disagrees"), c("Agrees", "Disagrees")))]
s <- s[order(regulator, -effect)]
s[, row := factor(paste(regulator, label), rev(paste(regulator, label)))]
LIMITS <- c(-0.17, 0.12)
block <- function(name, last) {
  d <- s[regulator == name][, row := droplevels(row)]
  p <- ggplot(d, aes(effect, row)) +
    geom_vline(xintercept = 0, linewidth = lw(0.5)) +
    geom_col(aes(fill = agree), width = 0.68, colour = EDGE, linewidth = lw(0.5)) +
    geom_linerange(aes(xmin = ci_low, xmax = ci_high), linewidth = lw(0.5), colour = EDGE) +
    geom_text(aes(x = fifelse(effect < 0, pmax(ci_high, 0), pmin(ci_low, 0)), label = fmt_p(q, prefix = "italic(q)"),
                  hjust = fifelse(effect < 0, -0.08, 1.08)), parse = TRUE, size = pt_size(6)) +
    facet_wrap(~regulator) +
    scale_fill_manual(values = c(Agrees = MOUSE, Disagrees = "white"), name = NULL, guide = if (last) "none" else "legend",
                      labels = c("Agrees with the prediction", "Disagrees")) +
    scale_x_continuous(labels = MINUS_LABELS, limits = LIMITS, expand = c(0, 0)) +
    scale_y_discrete(labels = setNames(scales::label_wrap(36)(d$label), d$row)) +
    labs(x = if (last) "Change of the m6A-marked target score" else NULL, y = NULL) +
    theme_box()
  if (last) p else p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}
rows <- s[, .N, keyby = regulator]$N
a <- block("FMR1 loss", FALSE) / block("PRRC2A loss", TRUE) + plot_layout(heights = rows, guides = "collect") &
  theme(legend.position = "bottom")

# b: median edge over the units of the primary source of each regulator, in units of the 95th percentile of the control
# proteins; *, m6A-specific edge in at least one unit; **, also shared across units; the YTHDF load edge counts as
# specific above the 95th percentile of the controls
m <- read_panel("sfig8b")
ROWS <- c("YTHDF1", "YTHDF2", "YTHDF3", "YTHDC1", "YTHDC2", "FMR1", "FXR1", "FXR2", "PRRC2A", "HNRNPC", "RBMX",
          "IGF2BP2", "YTHDF module", "FMRP family", "Writer complex", "Erasers",
          "Module, mouse", "Module, human", "YTHDF3, mouse", "YTHDF3, human")
m[, row := sub(" load,", ",", sub("^YTHDF module load", "Module", row))]
m[, `:=`(row = factor(row, rev(ROWS)), layer = factor(layer, LAYERS, c("Protein", "mRNA-adjusted\nprotein", "mRNA")),
         group = factor(group, c("Readers", "Modules", "Writers, erasers", "YTHDF load edge")),
         mark = fcase(shared, "**", m6A_specific, "*", default = ""))]
b <- ggplot(m, aes(layer, row, fill = value)) +
  geom_tile(colour = "white", linewidth = lw(0.6)) +
  geom_text(aes(label = number(value, 2)), size = pt_size(6)) +
  geom_text(aes(x = as.numeric(layer) + 0.47, label = mark), size = pt_size(6), hjust = 1, vjust = 0.3) +
  ggh4x::facet_grid2(group ~ ., scales = "free_y", space = "free_y", switch = "y", labeller = label_wrap_gen(9),
                     strip = ggh4x::strip_themed(background_y = ggh4x::elem_list_rect(fill = LIGHT, colour = LIGHT))) +
  scale_fill_signed(1.2, "Edge / control 95th percentile", light = TRUE) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_tile() +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 0, size = 6),
        panel.spacing.y = unit(1.2, "mm"), strip.switch.pad.grid = unit(1, "mm"), legend.position = "bottom")

# c: posterior mean of the edge of every regulator in every proteome unit after sharing across units (mashr); size,
# -log10 of the local false sign rate; open grey circles, lfsr >= 0.05
e <- read_panel("sfig8c")
STRIP <- tint(unname(SPECIES), 0.45)
REGS <- c("FMR1", "FMRP family", "FXR1", "FXR2", "PRRC2A", "HNRNPC", "RBMX", "YTHDC1", "YTHDF1", "YTHDF2", "YTHDF3",
          "YTHDF module", "Writer complex", "Erasers")
e[, `:=`(regulator = factor(regulator, REGS), species = factor(SPECIES_NAME[species], SPECIES_NAME))]
e[, unit := factor(unit, rev(unique(unit[order(species, unit)])))]
c_ <- fig_dotmatrix(e, "regulator", "unit", "posterior_mean", "lfsr", limits = c(-0.1, 0.1), cap = 6,
                    fill_lab = "Posterior mean edge", size_lab = expression(-log[10]~lfsr)) +
  ggh4x::facet_grid2(species ~ ., scales = "free_y", space = "free_y", switch = "y",
                     strip = ggh4x::strip_themed(background_y = ggh4x::elem_list_rect(fill = STRIP, colour = STRIP))) +
  scale_fill_signed(0.1, "Posterior mean edge", light = TRUE) +
  scale_size_area(max_size = 3.2, breaks = c(2, 4, 6), labels = c("2", "4", "≥6"), name = expression(-log[10]~lfsr),
                  guide = guide_legend(nrow = 1, title.position = "top", override.aes = list(fill = "#F7F7F7"))) +
  scale_x_discrete(position = "top", expand = expansion(add = 0.5)) +
  scale_y_discrete(expand = expansion(add = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5), strip.placement = "outside",
        strip.text.y.left = element_text(angle = 90, size = 6), strip.switch.pad.grid = unit(1, "mm"),
        panel.spacing.y = unit(1.4, "mm"), legend.position = "bottom", legend.box = "horizontal")

fig <- fig_compose(a, patchwork::wrap_elements(full = b), patchwork::wrap_elements(full = c_), design = "AC\nBC",
                   widths = c(80, 103), heights = c(88, 82))
save_fig(fig, "sfig8", 170)
