# Figure 4: the association is not explained by sequence, structure or cell composition
# reads source_data/fig4a-f (scripts 21 and 50); writes output/figures/fig4.png and .pdf
source("figures/style.R")

wrap_target <- function(x) gsub(" ", "\n", x, fixed = TRUE)

# a-c: age effect per log unit of the m6A load (or of the unmarked DRACH motifs given the load) under each check, by
# context (random effects over cohorts, 95% CI; filled, P < 0.05); rows from top to bottom
checks <- function(name, rows, xlab) {
  d <- read_panel(name)[check %in% names(rows)]
  d[, `:=`(panel = factor(panel, CONTEXTS), row = factor(rows[check], rev(rows)))]
  d[, key := KEYS[as.character(panel)]]
  d[, c("edge", "inside") := p_colours(p, key)]
  ggplot(d, aes(effect, row)) +
    ZERO_X +
    geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = CI, linewidth = 0.6) +
    geom_point(aes(colour = I(edge), fill = I(inside)), shape = 21, size = 3.8, stroke = 0.5) +
    facet_wrap(~panel, nrow = 1, scales = "free_x", labeller = label_wrap_gen(14), drop = TRUE) +
    scale_x_continuous(labels = MINUS_LABELS, breaks = scales::breaks_pretty(3)) +
    labs(x = xlab, y = NULL) +
    theme_rows() +
    BOXED +
    theme(panel.spacing = unit(4, "mm"), strip.text = element_text(lineheight = 0.95))
}
a <- checks("fig4a", c("m6A load" = "m6A load", "DRACH as covariate" = "+ all DRACH motifs",
                       "unmarked DRACH as covariate" = "+ unmarked DRACH motifs",
                       "unmarked DRACH given load" = "Unmarked motifs, given load"),
            "Age effect per log unit (%)")
b <- checks("fig4b", c("m6A load" = "m6A load, all maps", "GLORI sites" = "GLORI sites",
                       "GLORI stoichiometry" = "GLORI stoichiometry"),
            "Age effect of the m6A load per log unit (%)")
cc <- checks("fig4c", c("m6A load" = "All genes", "glia-enriched genes" = "Glia-enriched genes",
                        "ubiquitous genes" = "Ubiquitous genes"),
             "Age effect of the m6A load per log unit (%)")

# d: pooled load coefficient of 1,000 age permutations (histogram) and the observed coefficient (line); P, two-sided
# empirical (none of the 1,000 draws as extreme: P < 0.001)
null <- read_panel("fig4d")[, panel := factor(panel, CONTEXTS)]
obs <- read_panel("fig4d_observed")[, panel := factor(panel, CONTEXTS)]
dd <- ggplot(null, aes(effect_pct)) +
  geom_histogram(bins = 30, fill = "#E3E3E3", colour = NA) +
  geom_vline(data = obs, aes(xintercept = effect_pct, colour = I(LAYER[KEYS[as.character(panel)]])), linewidth = 1.2) +
  geom_text(data = obs, aes(x = Inf, y = Inf, label = fifelse(p_permutation < 1e-3, "P < 0.001", fmt_p(p_permutation))),
            hjust = 1.05, vjust = 1.3, size = PT, colour = INK) +
  facet_wrap(~panel, ncol = 2, scales = "free", labeller = label_wrap_gen(16)) +
  scale_x_continuous(labels = MINUS_LABELS, breaks = scales::breaks_pretty(3)) +
  scale_y_continuous(breaks = NULL, expand = expansion(mult = c(0, 0.25))) +
  labs(x = "Age effect of the m6A load per log unit (%)", y = NULL) +
  theme_fig() +
  BOXED +
  theme(panel.spacing = unit(4, "mm"), strip.text = element_text(lineheight = 0.95))

# e: gain in held-out R2 from the m6A site and binding features against 100 copies with the features permuted
e <- read_panel("fig4e")
e[, target := factor(target, TARGETS)]
obs <- unique(e[, .(target, observed_delta_r2, p)])
e <- ggplot(e, aes(target, permuted_delta_r2)) +
  ZERO_Y +
  geom_violin(fill = "#EDEDED", colour = "#BDBDBD", linewidth = 0.5, scale = "width", width = 0.8) +
  geom_point(data = obs, aes(y = observed_delta_r2), shape = 23, size = 4.6, fill = FOCAL, colour = EDGE,
             stroke = 0.5) +
  geom_text(data = obs, aes(y = Inf, label = fifelse(p < 0.0100001, "P ≤ 0.01", fmt_p(p))), vjust = 1.3, size = PT,
            colour = INK) +
  scale_x_discrete(labels = wrap_target) +
  scale_y_continuous(labels = MINUS_LABELS, expand = expansion(mult = c(0.05, 0.2))) +
  labs(x = NULL, y = "Gain in held-out R² from m6A sites and binding") +
  theme_fig() +
  theme(axis.text.x = element_text(lineheight = 0.9), axis.ticks.x = element_blank())

# f: in silico deletion of the m6A sites (all, or of one region): observed mean change of the prediction (point; filled,
# P < 0.05) against the 95% range of 100 permuted copies (bar)
f <- read_panel("fig4f")
f[, `:=`(target = factor(target, rev(TARGETS)), sig = factor(p < 0.05, c(TRUE, FALSE)),
         region = factor(region, c("all", "5UTR", "CDS", "3UTR"), c("All sites", "5′ UTR", "CDS", "3′ UTR")))]
f <- ggplot(f, aes(effect, target)) +
  ZERO_X +
  geom_linerange(aes(xmin = null_low, xmax = null_high), colour = LIGHT, linewidth = 4.4) +
  geom_point(aes(fill = sig, colour = sig), shape = 21, size = 3.4, stroke = 0.5) +
  facet_grid(. ~ region, scales = "free_x") +
  scale_fill_manual(values = c("TRUE" = FOCAL, "FALSE" = "white"), labels = c("TRUE" = "P < 0.05", "FALSE" = "n.s."),
                    name = NULL, drop = FALSE) +
  scale_colour_manual(values = c("TRUE" = EDGE, "FALSE" = GREY), labels = c("TRUE" = "P < 0.05", "FALSE" = "n.s."),
                      name = NULL, drop = FALSE) +
  scale_x_continuous(labels = MINUS_LABELS, breaks = function(l) pretty(l, n = 2)) +
  labs(x = "Change in predicted age z (with − without m6A sites)", y = NULL) +
  theme_rows() +
  BOXED +
  theme(legend.position = "bottom", legend.key.size = unit(4, "mm"), panel.spacing.x = unit(6, "mm"))

fig <- panel_grid(list(a = a, b = b, c = cc, d = dd, e = e, f = f), design = "AAB\nCCD\nEFF", widths = c(70, 52, 61),
                  heights = c(50, 54, 58))
save_fig(fig, "fig4", 324)
