# Fig. S6: m6A against the other marks (Ψ, m5C, m7G)
# reads source_data/sfig6a and sfig6b; writes output/figures/sfig6
source("figures/style.R")

MARKS <- c("m6A", "Ψ", "m5C", "m7G")
LIMIT <- 6
ROWS <- c("Mouse protein", "Mouse mRNA-adjusted protein", "Mouse mRNA", "Mouse neuronal mRNA", "Human protein",
          "Human striatum mRNA", "Human cortex mRNA", "Human hippocampus mRNA", "Human cerebellum mRNA",
          "Human neuronal mRNA")

# the ten contexts of the two panels, species in the strip and the unit in the row (first row on top)
NAMES <- sub("-adjusted ", "-adjusted\n", sub("^[A-Za-z]+ ", "", ROWS))
NAMES <- setNames(fifelse(startsWith(NAMES, "mRNA"), NAMES, paste0(toupper(substr(NAMES, 1, 1)), substring(NAMES, 2))),
                  ROWS)
read_marks <- function(name) {
  read_panel(name)[, `:=`(species = factor(sub(" .*", "", row), names(SPECIES)), mark = factor(modification, MARKS),
                          context = factor(row, rev(ROWS)))]
}
facets <- function(strip) {
  ggh4x::facet_grid2(species ~ ., scales = "free_y", space = "free_y", switch = "y",
                     strip = ggh4x::strip_themed(background_y = ggh4x::elem_list_rect(fill = strip, colour = strip)))
}
STRIP <- tint(unname(SPECIES), 0.45)
SPACING <- theme(panel.spacing.y = unit(1.4, "mm"), strip.placement = "outside", strip.switch.pad.grid = unit(1, "mm"))

# a: signature of the genes marked by each modification against unmarked genes (cameraPR z, combined over units); bold,
# q < 0.05
a <- read_marks("sfig6a")
a <- ggplot(a, aes(mark, context, fill = z)) +
  geom_tile(colour = "white", linewidth = lw(0.6)) +
  geom_text(aes(label = number(z, 1), fontface = fifelse(q < 0.05, "bold", "plain")), size = pt_size(6)) +
  facets(STRIP) +
  scale_fill_signed(LIMIT, "z", light = TRUE) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(labels = NAMES, expand = c(0, 0)) +
  theme_tile() +
  SPACING +
  theme(strip.text.y.left = element_text(angle = 90, size = 6))

# b: z of the age coefficient of each modification load per context (random effects over cohorts); size, -log10 q;
# open grey circles, q >= 0.05
b <- read_marks("sfig6b")
b <- fig_dotmatrix(b, "mark", "context", "z", "q", limits = c(-LIMIT, LIMIT), fill_lab = "z",
                   size_lab = expression(-log[10]~italic(q))) +
  facets(NA) +
  scale_fill_signed(LIMIT, "z", light = TRUE) +
  scale_size_area(max_size = 3.6, breaks = c(2, 5, 10), name = expression(-log[10]~italic(q)),
                  guide = guide_legend(nrow = 1, title.position = "top", override.aes = list(fill = "#F7F7F7"))) +
  scale_x_discrete(position = "top", expand = expansion(add = 0.5)) +
  scale_y_discrete(expand = expansion(add = 0.5)) +
  SPACING +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0), axis.text.y = element_blank(),
        strip.text.y.left = element_blank())

fig <- fig_compose(a, b, design = "AB", widths = c(1.25, 1), heights = 1, guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal", legend.box.spacing = unit(2, "mm"))
save_fig(fig, "sfig6", 84)
