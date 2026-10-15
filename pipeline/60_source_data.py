# Source data of Figures 1-5 and Additional file 1 Figs S1-S6: one small table per panel with the plotted values and
# their statistics, and Tables S1-S7 (Additional file 2); per-gene values only for the units shown
# reads the outputs of scripts 09, 10, 13, 20, 21, 25, 26, 28, 30, 33, 34, 37, 38 and 50, the prepared units, the site
# atlas, the regulator gene set and the CLIP provenance
# writes source_data/fig<N><panel>.tsv, sfig<N><panel>.tsv and supplementary_table<N>.tsv, replacing the tables there
import re

import numpy as np
import pandas as pd
from scipy import stats as sps

import utils

OUT = utils.ROOT / "source_data"
MODS = {"m6A": "m6A", "Psi": "Ψ", "m5C": "m5C", "m7G": "m7G"}
LAYERS = {"protein_only": "Protein", "mrna_adjusted": "mRNA-adjusted protein", "mrna_only": "mRNA"}
STUDIES = {"gtex": "GTEx", "brainseq2": "BrainSeq II", "brainseq3": "BrainSeq III", "chdi2016": "Langfelder 2016"}
STUDIES |= {"matlik2024": "Mätlik 2024", "colladotorres2019": "Collado-Torres 2019"}
PERCENT = {"effect_pct": "effect", "effect_pct_ci_low": "ci_low", "effect_pct_ci_high": "ci_high"}
ACRONYMS = dict(
    dlpfc="DLPFC", dpfc="DPFC", smtg="SMTG", ipl="IPL", ca1="CA1", ba9="BA9", ba24="BA24", c1="C1", nac="NAc", rna="RNA"
)


def result(script, name):
    return utils.read(utils.OUT / script / name)


def label(unit):
    # koopmans2026_dpfc_l1_3 -> Koopmans 2026 DPFC L1-3
    study, *words = unit.split("_")
    name = STUDIES.get(study) or re.sub(r"^(.)(\D*)(\d{4})$", lambda m: m[1].upper() + m[2] + " " + m[3], study)
    text = " ".join(ACRONYMS.get(word, word) for word in words if word not in ("bankadj", "expanded"))
    return f"{name} {re.sub(r'\bl(\d) (\d)', r'L\1-\2', text)}".strip()


def cohort_name(dataset):
    # benjamin2022_brainseq3_caudate -> Benjamin 2022
    return label(dataset.split("_")[0].replace("gtexv10", "gtex"))


def save(name, table):
    table.to_csv(OUT / f"{name}.tsv", sep="\t", index=False, float_format="%.4g")
    print(name, table.shape)


OUT.mkdir(exist_ok=True)
for stale in OUT.glob("*.tsv"):
    stale.unlink()
EFFECT_COLUMNS = ["unit", "gene_id", "symbol", "estimate", "t", "q", "effect_pct", "status"]


def age_effects(layer):
    # script 10: age effects of the genes each unit estimated
    table = utils.read(utils.OUT / "10_protein_age_effects" / f"{layer}.tsv.gz", usecols=EFFECT_COLUMNS)
    return table[table.status == "ok"]


effects = {layer: age_effects(layer) for layer in LAYERS}
cohort_of = utils.read(utils.PREPARED / "units.tsv", index_col="unit").dataset  # dataset of each unit

# Figure 1d: human protein age effects against mouse age effects on each layer
concordance = result("13_species_concordance", "species_concordance.tsv").query("human_layer == 'protein_only'")
mouse = {"chdi2016_expanded": "Mouse striatum (Langfelder 2016)", "takasugi2024": "Mouse cerebrum (Takasugi 2024)"}
fig1d = concordance.assign(
    human_unit=concordance.human_unit.map(label),
    mouse_cohort=concordance.mouse_unit.map(mouse),
    mouse_layer=concordance.mouse_layer.map(LAYERS),
)
assert not fig1d.duplicated(["human_unit", "mouse_cohort", "mouse_layer"]).any()
save("fig1d", fig1d[["human_unit", "mouse_cohort", "mouse_layer", "rho", "ci_low", "ci_high", "p"]])

# 1c: age t of the orthologs of one pair: Koopmans 2026 DPFC L1-3 protein, Takasugi 2024 mRNA-adjusted protein
PAIR = ("koopmans2026_dpfc_l1_3", "takasugi2024")
protein = effects["protein_only"]
human_t = utils.by_symbol(protein[protein.unit == PAIR[0]]).t
orthologs = utils.orthologs().set_index("human").mouse
human_t = human_t[human_t.index.isin(orthologs.index)].rename(index=orthologs)
paired = {layer: utils.by_symbol(table[table.unit == PAIR[1]]) for layer, table in effects.items()}
paired = {layer: table[["estimate", "t", "q"]] for layer, table in paired.items()}
paired = pd.concat(paired, axis=1, join="inner").dropna()  # genes of all three mouse layers, as in script 13
pairs = pd.concat({"human_t": human_t, "mouse_t": paired["mrna_adjusted", "t"]}, axis=1, join="inner")
row = concordance[(concordance.human_unit == PAIR[0]) & (concordance.mouse_unit == PAIR[1])]
row = row[row.mouse_layer == "mrna_adjusted"].iloc[0]
assert len(pairs) == row.n and np.isclose(sps.spearmanr(pairs.human_t, pairs.mouse_t).statistic, row.rho)
save("fig1c", pairs.rename_axis("gene").reset_index())

CONTEXT_UNITS = {(context, layer): listed.split() for context, layer, listed in utils.CONTEXTS}
screen = result("30_regulator_screen", "regulator_effects.tsv")
regulators = utils.read(utils.ROOT / "ref/gene_sets/rna_mod_regulators.tsv")
net = result("33_regulator_network", "edges.tsv")


def bridge(name):
    return result("38_reader_bridges", name)


CONTEXT_ROWS = {  # context and layer of the Additional file 1 panels: row label
    ("mouse_protein", "protein_only"): "Mouse protein",
    ("mouse_protein", "mrna_adjusted"): "Mouse protein given mRNA",
    ("mouse_mrna", "mrna_only"): "Mouse mRNA",
    ("mouse_neurons", "mrna_only"): "Mouse neuronal mRNA",
    ("human_protein", "protein_only"): "Human protein",
    ("human_striatum", "mrna_only"): "Human striatum mRNA",
    ("human_cortex", "mrna_only"): "Human cortex mRNA",
    ("human_hippocampus", "mrna_only"): "Human hippocampus mRNA",
    ("human_cerebellum", "mrna_only"): "Human cerebellum mRNA",
    ("human_neurons", "mrna_only"): "Human neuronal mRNA",
}


def rows_of(table):
    # the contexts and layers shown in the Additional file 1 panels, with their row label
    row = pd.Series(list(zip(table.context, table.layer, strict=True)), index=table.index).map(CONTEXT_ROWS)
    return table.assign(row=row)[row.notna()]


def site_counts(layer):
    # mean age effect of the genes by number of m6A sites (script 21)
    return result("21_modification_load", "site_count_effects.tsv").query("modification == 'm6A' and layer == @layer")


# Figure 2a: m6A signature (cameraPR z) of every mouse unit
sig = result("20_modification_signatures", "unit_signatures.tsv").query("species == 'mouse' and modification == 'm6A'")
save("fig2a", sig.assign(unit=sig.unit.map(label), layer=sig.layer.map(LAYERS))[["unit", "context", "layer", "z", "p"]])

# 2b: mean age effect of the protein by number of m6A sites, over the units of a species with a 95% interval across
# units (human in Figure 3a)
sites = site_counts("protein_only")
rows = []
for (species, count), group in sites.groupby(["species", "sites"], sort=False):
    mean, sem = group.effect_pct.mean(), group.effect_pct.sem()
    width = sps.t.ppf(0.975, len(group) - 1) * sem
    rows.append(
        dict(species=species, sites=count, n_units=len(group), effect=mean, ci_low=mean - width, ci_high=mean + width)
    )
by_sites = pd.DataFrame(rows)
save("fig2b", by_sites.query("species == 'mouse'"))

# 2b, right: the same for the protein given its mRNA, one line per unit (human units in Figure 3a)
adjusted = site_counts("mrna_adjusted").rename(columns=PERCENT)
adjusted = adjusted.assign(unit=adjusted.unit.map(label))
adjusted = adjusted[["species", "unit", "sites", "n_genes", "effect", "ci_low", "ci_high"]]
save("fig2b_given_mrna", adjusted.query("species == 'mouse'"))

# 2c: age effect of every protein of the largest wild-type mouse proteome cohort (Wang 2025, 40 mice) with the number
# of m6A sites of its gene
features_mouse = utils.OUT / "09_transcript_features" / "mouse_features.tsv.gz"
sites_mouse = utils.read(features_mouse, usecols=["gene_id", "n_sites_m6A"])
wang = effects["protein_only"].query("unit == 'wang2025_brain'").merge(sites_mouse, on="gene_id")
save("fig2c", wang[["gene_id", "symbol", "effect_pct", "n_sites_m6A"]])
unit_load = result("21_modification_load", "modification_load_effects.tsv")
unit_load = unit_load.query("unit == 'wang2025_brain' and modification == 'm6A' and layer == 'protein_only'")
unit_load = unit_load[unit_load.model == "without protein lifetime"]
save("fig2c_load", unit_load.assign(unit=unit_load.unit.map(label)).rename(columns=PERCENT)[["unit", "effect", "p"]])

# 2d: m6A load with the covariates added in turn (the fifth step is the primary model, the sixth adds the other loads)
PANELS = {
    ("mouse_protein", "protein_only"): "Mouse protein",
    ("mouse_protein", "mrna_adjusted"): "Mouse protein given mRNA",
    ("human_protein", "protein_only"): "Human protein",
    ("human_striatum", "mrna_only"): "Human striatum mRNA",
}


def panel_of(table):
    # the four contexts of Figures 2d, 3a and 4, as the label of the panel; the others are None
    return [PANELS.get(key) for key in zip(table.context, table.layer, strict=True)]


ladder = result("21_modification_load", "covariate_ladder.tsv")
ladder["panel"] = panel_of(ladder)
ladder = ladder.dropna(subset="panel").rename(columns=PERCENT)
save("fig2d", ladder[["panel", "step", "n_cohorts", "effect", "ci_low", "ci_high", "p", "retained_pct"]])

# 2e: average partial effect of the m6A load inside and outside the long-lived proteins and their difference (causal
# forests, mouse protein; the other classes are in Supplementary Data 1)
forest = result("25_load_heterogeneity", "combined_effects.tsv")
forest = forest.query("context == 'mouse_protein' and layer == 'protein_only' and gene_class == 'long-lived'")
columns = ["gene_class", "quantity", "n_cohorts", "effect", "ci_low", "ci_high", "p", "q", "i2"]
save("fig2e", forest.rename(columns=PERCENT)[columns])

# Figure 3a: mean age effect of human proteins by number of m6A peaks of the brain map, and the protein given its mRNA
save("fig3a", by_sites.query("species == 'human'"))
save("fig3a_given_mrna", adjusted.query("species == 'human'"))

# 3b: m6A load of the human units and its random-effects estimate over cohorts, brain m6A map
HUMAN = ["human_protein", "human_striatum", "human_cortex", "human_hippocampus", "human_cerebellum", "human_neurons"]
combined = result("21_modification_load", "combined_load_effects.tsv")
combined = combined[(combined.model == "without protein lifetime") & (combined.modification == "m6A")]
combined = combined[combined.context.isin(HUMAN) & (combined.layer != "mrna_adjusted")].rename(columns=PERCENT)
columns = ["context", "layer", "n_cohorts", "effect", "ci_low", "ci_high"]
save("fig3b", combined.assign(layer=combined.layer.map(LAYERS))[columns])
per_unit = result("21_modification_load", "modification_load_effects.tsv")
per_unit = per_unit[(per_unit.modification == "m6A") & (per_unit.species == "human")]
per_unit = per_unit[(per_unit.layer != "mrna_adjusted") & (per_unit.model == "without protein lifetime")]
human_sig = result("20_modification_signatures", "unit_signatures.tsv")
human_sig = human_sig.query("modification == 'm6A' and species == 'human'")
per_unit = per_unit.merge(human_sig[["unit", "layer", "context"]], on=["unit", "layer"])
per_unit = per_unit[per_unit.context.isin(HUMAN)]
per_unit = per_unit.assign(unit_id=per_unit.unit, unit=per_unit.unit.map(label), layer=per_unit.layer.map(LAYERS))
columns = ["context", "unit", "layer", "effect", "unit_id"]
save("fig3b_units", per_unit.rename(columns={"effect_pct": "effect"})[columns])

# 3c: validation tiers of the m6A load and signature
tiers = result("26_validation_tiers", "validation_tiers.tsv").query("test != 'YTHDF3 level'")
save("fig3c", tiers[["tier", "test", "validation", "validation_z", "p_one_sided", "passed"]])

# 3d: human m6A load counted from the brain map or the cell-line maps
maps = result("28_map_sensitivity", "combined_load_effects.tsv")
maps = maps[maps.context.isin(HUMAN) & (maps.layer != "mrna_adjusted")]
save("fig3d", maps[["context", "layer", "definition", "n_cohorts", "effect_pct", "ci_low_pct", "ci_high_pct"]])

# Figure 4: is it m6A? a-c, load checks of script 21 (DRACH motifs; GLORI stoichiometry, mouse; cell composition) beside
# the m6A load of the primary model; d, the age-permutation null of the load coefficient; e, f, gain and in silico
# deletion of the m6A features in the models of gene-level ageing (script 50)
main = result("21_modification_load", "combined_load_effects.tsv")
main = main.query("model == 'without protein lifetime' and modification == 'm6A'").assign(check="m6A load")
checks = pd.concat([main, result("21_modification_load", "load_checks.tsv")])
checks["panel"] = panel_of(checks)
checks = checks.dropna(subset="panel").rename(columns=PERCENT)
SHOWN = {
    "fig4a": ["m6A load", "DRACH as covariate", "unmarked DRACH as covariate", "unmarked DRACH given load"],
    "fig4b": ["m6A load", "GLORI sites", "GLORI stoichiometry"],
    "fig4c": ["m6A load", "glia-enriched genes", "ubiquitous genes"],
}
for name, shown in SHOWN.items():
    rows = checks[checks.check.isin(shown) & (checks.panel.str.startswith("Mouse") | (name != "fig4b"))]
    save(name, rows[["panel", "check", "n_cohorts", "effect", "ci_low", "ci_high", "p", "i2"]])
null = result("21_modification_load", "permutation_null_draws.tsv")
null["panel"] = panel_of(null)
save("fig4d", null.dropna(subset="panel")[["panel", "draw", "effect_pct"]])
observed = result("21_modification_load", "permutation_null.tsv")
observed["panel"] = panel_of(observed)
save("fig4d_observed", observed.dropna(subset="panel")[["panel", "n_draws", "effect_pct", "p", "p_permutation"]])


def model(name):
    return result("50_transcript_ageing_model", name)


TARGETS = {"mouse_protein": "Mouse protein", "human_mrna": "Human mRNA", "human_neurons": "Human neuronal mRNA"}
TARGETS |= {"mouse_mrna_adjusted": "Mouse protein given mRNA", "human_protein": "Human protein"}
ablation = model("block_ablation.tsv").dropna(subset=["p"])
observed_gain = ablation[["target", "delta_r2", "p"]].rename(columns={"delta_r2": "observed_delta_r2"})
draws = model("block_ablation_permuted.tsv").merge(observed_gain)
save("fig4e", draws.assign(target=draws.target.map(TARGETS)).rename(columns={"delta_r2": "permuted_delta_r2"}))
deletion = model("deletion.tsv").query("sites == 'main'")
columns = ["target", "region", "n_genes", "effect", "null_low", "null_high", "p"]
save("fig4f", deletion.assign(target=deletion.target.map(TARGETS))[columns])

# Figure 5b: age effect of every m6A regulator in mouse protein with its q over all regulators of the context
m6a_roles = regulators[regulators.modification.str.split(";").apply(lambda m: "m6A" in m)]


def role_of(roles):
    # a regulator with several roles is a reader first, then an eraser, else a writer (the writer complex included)
    return "reader" if "reader" in set(roles) else "eraser" if "eraser" in set(roles) else "writer"


m6a_roles = m6a_roles.groupby("human_symbol").role.agg(role_of).rename("role")
volcano = screen.query("context == 'mouse_protein' and layer == 'protein_only'").drop(columns=["ci_low", "ci_high"])
volcano = volcano.merge(m6a_roles, left_on="regulator", right_index=True)
save("fig5b", volcano.rename(columns=PERCENT)[["regulator", "role", "effect", "ci_low", "ci_high", "p", "q"]])

# 5c: median edge of the regulators and modules over the units of both species with their primary source, in
# units of the 95th percentile of the control RBPs (Additional file 1 Fig. S5b adds the YTHDF load edges of
# script 38)
NETWORK = ["YTHDF1", "YTHDF2", "YTHDF3", "YTHDC1", "YTHDC2", "FMR1", "FXR1", "FXR2", "PRRC2A", "HNRNPC", "RBMX"]
NETWORK += ["IGF2BP2"]
MODULES = {"YTHDF": "YTHDF module", "FMRP_family": "FMRP family", "writer_complex": "Writer complex"}
MODULES |= {"erasers": "Erasers"}
ARMS = ["Output", "Translation", "Stability"]


def primary_edges(edges, regulators):
    # the edges of the primary source of each regulator and species; FMR1: the brain CLIP of Darnell 2011 in both
    # species (human units through orthologs)
    darnell = edges.source.isin(["Darnell2011_brain", "Darnell2011_brain_mapped"])
    keep = np.where(edges.regulator == "FMR1", darnell, edges.primary.eq(True))
    return edges[keep & edges.regulator.isin(regulators)]


edges = primary_edges(net, [*NETWORK, *MODULES])
edges = edges.assign(
    value=edges.b_s / edges.control_p95,
    m6A_specific=edges.edge_class == "m6A-specific",
    binding_only=edges.edge_class == "Binding-only",
)
edges["shared"] = edges.m6A_specific & edges.shared_edge.fillna(False).astype(bool)
cells = edges.groupby(["regulator", "family", "arm"]).agg(
    value=("value", "median"),
    m6A_specific=("m6A_specific", "any"),
    shared=("shared", "any"),
    binding_only=("binding_only", "any"),
).reset_index()
GROUPS = {"reader": "Readers", "module": "Modules", "writer": "Writers, erasers", "eraser": "Writers, erasers"}
cells["group"] = cells.family.map(GROUPS)
cells["row"] = cells.regulator.replace(MODULES)
pooled = bridge("yth_coupling_pooled.tsv").query(
    "target_type == 'm6A_load' and regulator in ['YTHDF', 'YTHDF3'] and arm in @ARMS")
pooled = pooled.assign(
    row=pooled.regulator.map({"YTHDF": "YTHDF module", "YTHDF3": "YTHDF3"}) + " load, " + pooled.species,
    value=pooled.estimate / pooled.control_abs_p95,
    m6A_specific=pooled.pctile_signed >= 95,
    shared=False,
    binding_only=False,
    group="YTHDF load edge",
)
columns = ["group", "row", "arm", "value", "m6A_specific", "shared", "binding_only"]
edge_values = pd.concat([cells, pooled])[columns]
save("fig5c", edge_values[edge_values.group != "YTHDF load edge"])

# 5d: FMR1 edge to its m6A-marked (s) and unmarked (u) targets of the two brain CLIP lists (mouse; human through
# orthologs) in every unit and arm, and in the striatal mRNA unit of the allelic series refitted on its wild-type mice
FMRP_LISTS = {"Darnell2011_brain": "Darnell 2011", "Maurin2018_brain": "Maurin 2018"}
FMRP_LISTS |= {f"{source}_mapped": name for source, name in FMRP_LISTS.items()}
fmr1 = net[(net.regulator == "FMR1") & net.source.isin(FMRP_LISTS)].assign(mice="all")
wildtype = result("33_regulator_network", "wildtype_edges.tsv").query("regulator == 'FMR1' and source in @FMRP_LISTS")
half = utils.Z95 * wildtype.se_s
wildtype = wildtype.assign(ci_low=wildtype.b_s - half, ci_high=wildtype.b_s + half)
fmr1 = pd.concat([fmr1, wildtype]).assign(list=lambda d: d.source.map(FMRP_LISTS), unit=lambda d: d.unit.map(label))
save("fig5d", fmr1[["list", "species", "arm", "unit", "mice", "b_s", "ci_low", "ci_high", "b_u", "control_p95"]])

# 5e: percentile of the YTHDF module and YTHDF3 load edge among the load edges of the control RBPs (pooled over units),
# and in the striatal mRNA unit of the allelic series refitted on its wild-type mice
hits = bridge("yth_coupling_pooled.tsv").query("target_type == 'm6A_load' and regulator in ['YTHDF', 'YTHDF3']")
wild = bridge("yth_wildtype_coupling.tsv").assign(arm="Stability, wild-type mice")
save("fig5e", pd.concat([hits, wild])[["species", "arm", "regulator", "pctile_signed", "p_control_two_sided"]])

# Additional file 1 Fig. S1: m6A against the other marks and the other marks beside m6A
signatures = rows_of(result("20_modification_signatures", "combined_signatures.tsv"))
save("sfig1a", signatures.assign(modification=signatures.modification.map(MODS))[["row", "modification", "z", "q"]])
loads = result("21_modification_load", "combined_load_effects.tsv").query("model == 'without protein lifetime'")
loads = rows_of(loads)
save("sfig1b", loads.assign(modification=loads.modification.map(MODS))[["row", "modification", "z", "q"]])

# Fig. S2: genes with an m6A site in each map (a, the two human maps; b, the intersections of the four mouse maps), and
# the m6A load and signature by site map, human and mouse
MAP_NAMES = {
    "m6a_brain_shafik2026": "brain peaks",
    "m6a_hek293t_hela": "cell line",
    "m6a_hippocampus_cyto": "hippocampus cyto",
    "m6a_hippocampus_syn": "hippocampus syn",
    "m6a_cortex_scdart_tegowski2024": "cortex scDART",
    "m6a_neuron_dart_tegowski2024": "neuron DART",
}
sites = utils.read(utils.RESOURCES / "site_atlas" / "sites.parquet", columns=["map", "gene_id", "ref_ok"])
sites = sites[sites.ref_ok & sites["map"].isin(MAP_NAMES) & sites.gene_id.notna()]
genes = sites.groupby("map").gene_id.apply(set).rename(index=MAP_NAMES)
brain, cell = genes["brain peaks"], genes["cell line"]
regions = ["brain peaks only", "both", "cell line only"]
save("sfig2a", pd.DataFrame({"region": regions, "genes": [len(brain - cell), len(brain & cell), len(cell - brain)]}))
MOUSE_MAPS = ["hippocampus cyto", "hippocampus syn", "cortex scDART", "neuron DART"]
membership = pd.DataFrame({name: pd.Series(True, index=sorted(genes[name])) for name in MOUSE_MAPS}).fillna(False)
save("sfig2b", membership.value_counts().rename("genes").reset_index())
load_maps = result("28_map_sensitivity", "combined_load_effects.tsv")
signature_maps = result("28_map_sensitivity", "combined_signatures.tsv")
# c, the human signature by map; d, e, the mouse load and signature by map (the human load is Figure 3d)
for name, table, species, columns in (
    ("c", signature_maps, "human", ["z", "q"]),
    ("d", load_maps, "mouse", ["effect_pct", "ci_low_pct", "ci_high_pct"]),
    ("e", signature_maps, "mouse", ["z", "q"]),
):
    shown = rows_of(table[table.context.str.startswith(species)])
    save(f"sfig2{name}", shown[["context", "layer", "definition", *columns]])

# f: the four mouse m6A maps: sites, share in a DRACH motif, share placed on a canonical transcript, sites per marked
# gene and share at a GLORI position
columns = ["map", "chrom", "start", "strand", "motif_class", "gene_id", "tx_offset", "ref_ok"]
atlas = utils.read(utils.RESOURCES / "site_atlas" / "sites.parquet", columns=columns)
atlas = atlas[atlas["map"].isin([m for m, n in MAP_NAMES.items() if n in MOUSE_MAPS])]
glori = set(map(tuple, atlas[atlas["map"].str.contains("hippocampus")][["chrom", "start", "strand"]].to_numpy()))
quality = []
for name, group in atlas.groupby("map"):
    placed = group[group.gene_id.notna() & group.tx_offset.notna() & group.ref_ok]
    at_glori = np.mean([key in glori for key in map(tuple, group[["chrom", "start", "strand"]].to_numpy())])
    quality.append(
        dict(
            map=MAP_NAMES[name],
            sites=len(group),
            drach=(group.motif_class == "DRACH").mean(),
            placed=len(placed) / len(group),
            sites_per_gene=len(placed) / placed.gene_id.nunique(),
            at_glori=at_glori,
        )
    )
save("sfig2f", pd.DataFrame(quality))

# Fig. S3: the models of gene-level ageing; b, held-out rho of every model, and of the mouse protein model applied to
# human proteins (transfer)
performance = model("performance.tsv")
transfer = model("transfer.tsv").query("target == 'human_protein'")
transfer = transfer.assign(
    target="Human protein (mouse model)",
    model=transfer.model.map({"shared": "tabpfn", "shared_shuffled": "shuffled_target"}),
)
columns = ["target", "model", "n", "rho", "rho_ci_low", "rho_ci_high"]
save("sfig3b", pd.concat([performance.assign(target=performance.target.map(TARGETS)), transfer])[columns])
# c, per-gene effect of deleting all sites by the number of sites of the gene (mean, 95% CI)
genes, rows = model("deletion_genes.tsv.gz"), []
for target in ("mouse_protein", "human_protein"):
    inputs = utils.OUT / "50_transcript_ageing_model" / "inputs" / f"{target}.parquet"
    sites = utils.read(inputs, columns=["m6A_all"])
    effect = genes[genes.target == target].set_index("gene_id").effect_all
    count = np.expm1(sites.m6A_all.reindex(effect.index)).round()
    classes = pd.cut(count, [0, 1, 3, 7, np.inf], labels=["1", "2-3", "4-7", ">=8"])
    for bin_, values in effect.groupby(classes, observed=True):
        mean, half = values.mean(), utils.Z95 * values.std() / np.sqrt(len(values))
        rows.append(
            dict(target=TARGETS[target], sites=bin_, n=len(values), mean=mean, ci_low=mean - half, ci_high=mean + half)
        )
save("sfig3c", pd.DataFrame(rows))

# Fig. S4: YTHDF3 and the YTHDF readers; a, YTHDF3 age effect, mouse protein per cohort and combined, human mRNA per
# context
cohorts = result("30_regulator_screen", "regulator_cohorts.tsv")
cohorts = cohorts.query("regulator == 'YTHDF3' and context == 'mouse_protein' and layer == 'protein_only'")
cohorts = pd.DataFrame(
    dict(
        panel="Mouse protein",
        label=cohorts.cohort.map(cohort_name),
        effect=utils.pct(cohorts.estimate),
        ci_low=utils.pct(cohorts.estimate - utils.Z95 * cohorts.se),
        ci_high=utils.pct(cohorts.estimate + utils.Z95 * cohorts.se),
    )
)
ythdf3 = screen.query("regulator == 'YTHDF3'").drop(columns=["ci_low", "ci_high"])
mouse_protein = (ythdf3.context == "mouse_protein") & (ythdf3.layer == "protein_only")
ythdf3 = ythdf3[mouse_protein | (ythdf3.context.isin(HUMAN[1:]) & (ythdf3.layer == "mrna_only"))]
mouse_row = ythdf3.context == "mouse_protein"
ythdf3 = ythdf3.assign(
    panel=np.where(mouse_row, "Mouse protein", "Human mRNA"),
    label=np.where(mouse_row, "Combined", ythdf3.context.str.removeprefix("human_")),
).rename(columns=PERCENT)
save("sfig4a", pd.concat([cohorts, ythdf3[["panel", "label", "effect", "ci_low", "ci_high", "q"]]]))

# b: pooled load edge of the YTHDF module and YTHDF3 against the control RBPs (random effects over units)
loads = bridge("yth_coupling_edges.tsv").query("target_type == 'm6A_load'")
controls = loads[loads.is_control].groupby(["species", "arm", "regulator"]).apply(
    lambda d: utils.random_effects(d.b_s, d.se_s)["est"] if len(d) > 1 else d.b_s.iloc[0], include_groups=False
)
save("sfig4b_controls", controls.rename("edge").reset_index())
save("sfig4b", hits[["species", "arm", "regulator", "estimate", "pctile_signed"]])

# c: fraction of the age slope of the m6A load taken over by the YTHDF terms, mouse units
mediation = bridge("yth_mediation_units.tsv")
mediation = mediation[(mediation.target_type == "m6A_load") & (mediation["mode"] == "main")]
mediation = mediation[mediation.regulator.isin(["YTHDF", "YTHDF3"]) & mediation.arm.isin(ARMS)]
mediation = mediation[mediation.species == "mouse"]
percent = {f"{c}_pct": 100 * mediation[c] for c in ["fraction", "ci_low", "ci_high"]}
mediation = mediation.assign(unit=mediation.unit.map(label), **percent)
save("sfig4c", mediation[["unit", "arm", "regulator", "fraction_pct", "ci_low_pct", "ci_high_pct"]])

# d: change of the target score after loss of a YTHDF reader (log2 intervals converted to %)
loss = bridge("yth_perturbation.tsv").query("primary == True and predicted != 0 and effect_pct.notna()")
loss = loss.assign(ci_low_pct=utils.pct(loss.ci_low), ci_high_pct=utils.pct(loss.ci_high))
columns = ["dataset", "reader", "target_reader", "readout", "comparison", "effect_pct", "ci_low_pct", "ci_high_pct"]
save("sfig4d", loss[[*columns, "sign_matches"]])

# Fig. S5: FMRP family and PRRC2A; a, sign check against knockouts; b, the edge values of Figure 5c with the YTHDF load
# edges; c, edges shared across units
check = result("34_perturbation_signs", "sign_check.tsv")
check = check.query("regulator in ['FMR1', 'PRRC2A'] and layer == 'protein_only'")
columns = ["perturbation", "regulator", "tissue", "readout", "effect", "ci_low", "ci_high", "q", "agrees"]
save("sfig5a", check[columns])
save("sfig5b", edge_values)
SHARED = [*NETWORK[5:11], "YTHDC1", "YTHDF1", "YTHDF2", "YTHDF3", "FMRP_family", "YTHDF", "writer_complex", "erasers"]
shared = primary_edges(net[net.arm == "Output"], SHARED)
shared = shared.assign(unit=shared.unit.map(label), regulator=shared.regulator.replace(MODULES))
save("sfig5c", shared[["regulator", "unit", "species", "posterior_mean", "lfsr"]])

# Fig. S6: the targets of IGF2BP2 do not follow its decline; a, IGF2BP2 mRNA and the target scores of its primary
# tiers in GTEx (region x age) and the single-nucleus cohorts
regions = bridge("igf2bp2_striatum_x_age.tsv")
regions = regions[~regions.is_control & regions.source.isin(["IGF2BP2:primary", "IGF2BP2 mRNA"])]
nuclei = bridge("igf2bp2_nuclei_combined.tsv").query("not is_control")
columns = ["panel", "tier", "comparison", "estimate_pct", "ci_low_pct", "ci_high_pct"]
save("sfig6a", pd.concat([regions.assign(panel="GTEx"), nuclei.assign(panel="Single nucleus")])[columns])

# b: genes bound with m6A in each combination of the four neural binding sources (T1, at least two)
NEURAL = ["conway2016_H9_eCLIP", "degrauwe2016_MGG8_PARCLIP", "han2022_A172_eCLIP", "preitner2016_E14brain_CLIP"]
bound = result("37_reader_compendium", "bound_genes.tsv").query("reader == 'IGF2BP2' and T1 and source in @NEURAL")
membership = pd.DataFrame(
    {source: pd.Series(True, index=bound[bound.source == source].human_gene_id.dropna().unique()) for source in NEURAL}
).fillna(False)
save("sfig6b", membership.value_counts().rename("genes").reset_index())
coupling = bridge("igf2bp2_coupling.tsv")
columns = ["dataset", "population", "b_s", "ci_low", "ci_high", "pctile_abs"]
save("sfig6c", coupling.assign(dataset=coupling.dataset.map(cohort_name))[columns])
mediated = bridge("igf2bp2_mediation.tsv")
mediated = mediated.assign(dataset=mediated.dataset.map(cohort_name), fraction_pct=100 * mediated.fraction)
save("sfig6d", mediated[["dataset", "population", "tier", "comparison", "fraction_pct"]])
protein_score = bridge("igf2bp2_proteome_pooled.tsv").query("not is_control and source == 'IGF2BP2:primary'")
save("sfig6e", protein_score[["species", "tier", "comparison", "estimate_pct", "ci_low_pct", "ci_high_pct"]])
loss = bridge("igf2bp2_perturbation.tsv")
loss = loss.assign(ci_low_pct=utils.pct(loss.ci_low), ci_high_pct=utils.pct(loss.ci_high))
save("sfig6f", loss[["dataset", "tier", "comparison", "effect_pct", "ci_low_pct", "ci_high_pct", "q"]])
neuron_glia = bridge("igf2bp2_neuron_glia_combined.tsv").query("not is_control")
save("sfig6g", neuron_glia[["set", "comparison", "pctile_z"]])

# Supplementary Table 1 and Figure 1a: analysis units with dataset, layers, samples, ages, model terms and contexts
units = utils.read(utils.PREPARED / "units.tsv", keep_default_na=False)
sources = utils.read(utils.ROOT / "data_sources.tsv", keep_default_na=False)
doi = sources.assign(dataset=sources.id.str.split("/").str[0]).drop_duplicates("dataset").set_index("dataset").doi
contexts = {}
for context, layer, listed in utils.CONTEXTS:
    for unit in listed.split():
        contexts.setdefault(unit, []).append(f"{context} ({LAYERS[layer]})")
rows = []
for unit in units.itertuples():
    folder = utils.PREPARED / unit.unit
    samples = utils.read(folder / "samples.tsv")
    species = "mouse" if "age_months" in samples else "human"
    ages = samples.age_months if species == "mouse" else samples.filter(["age_years", "age_midpoint"]).iloc[:, 0]
    layers = [name for name, file in (("protein", "protein"), ("RNA", "rna")) if (folder / f"{file}.parquet").exists()]
    rows.append(
        dict(unit=label(unit.unit), unit_id=unit.unit, dataset=cohort_name(unit.dataset), species=species)
        | dict(layers=" and ".join(layers), samples=len(samples), donors=samples.donor_id.nunique())
        | dict(age_min=ages.min(), age_max=ages.max(), age_unit="months" if species == "mouse" else "years")
        | dict(age_term=utils.AGE_SCALE[species], covariates=",".join(filter(None, [unit.numeric, unit.categorical])))
        | dict(hidden_factors=unit.hidden_factors, permutation_blocks=unit.blocks)
        | dict(contexts="; ".join(contexts.get(unit.unit, [])), doi=doi.get(unit.dataset, ""))
    )
table1 = pd.DataFrame(rows)
save("supplementary_table1", table1)

# 1a: one row per dataset and molecular layer (protein, bulk RNA, single-nucleus RNA): units, donors of its largest unit
# and age range; mouse bulk and single-nucleus RNA make one group
NUCLEI = CONTEXT_UNITS["human_neurons", "mrna_only"] + CONTEXT_UNITS["mouse_neurons", "mrna_only"]
GROUPS = {
    ("human", "Protein"): "Human protein",
    ("human", "Bulk RNA"): "Human bulk RNA",
    ("human", "Single-nucleus RNA"): "Human single-nucleus RNA",
    ("mouse", "Protein"): "Mouse protein",
    ("mouse", "Bulk RNA"): "Mouse RNA",
    ("mouse", "Single-nucleus RNA"): "Mouse RNA",
}
long = table1.assign(layer=table1.layers.str.split(" and ")).explode("layer")
long["layer"] = np.where(
    long.layer == "protein", "Protein", np.where(long.unit_id.isin(NUCLEI), "Single-nucleus RNA", "Bulk RNA")
)
long = long.assign(group=[GROUPS[key] for key in zip(long.species, long.layer, strict=True)])
long["dataset_id"] = long.unit_id.map(cohort_of)
fig1a = long.groupby(["group", "layer", "dataset_id"], sort=False).agg(
    species=("species", "first"),
    n_units=("unit_id", "size"),
    donors=("donors", "max"),
    age_min=("age_min", "min"),
    age_max=("age_max", "max"),
    age_unit=("age_unit", "first"),
)
fig1a = fig1a.reset_index().assign(dataset=lambda d: d.dataset_id.map(cohort_name))
twice = fig1a.duplicated(["group", "dataset"], keep=False)  # Hahn 2023: bulk atlas and nuclei
fig1a["dataset"] += np.where(twice, fig1a.layer.map({"Bulk RNA": " bulk", "Single-nucleus RNA": " nuclei"}), "")
assert fig1a.dataset_id.nunique() == units.dataset.nunique() and not fig1a.duplicated(["group", "dataset"]).any()
save("fig1a", fig1a[["group", "layer", "dataset", "n_units", "donors", "age_min", "age_max", "age_unit"]])

# Supplementary Table 2: modification site maps with their source, resolution and numbers of sites and genes
columns = ["map", "modification", "species", "resolution", "source", "gene_id", "tx_offset", "ref_ok"]
atlas = utils.read(utils.RESOURCES / "site_atlas" / "sites.parquet", columns=columns)
first = {column: (column, "first") for column in ["modification", "species", "resolution", "source"]}
maps = atlas.groupby("map").agg(**first, sites=("map", "size"))
placed = atlas[atlas.gene_id.notna() & atlas.tx_offset.notna() & atlas.ref_ok]
maps["genes"] = placed.groupby("map").gene_id.nunique()
maps = maps.assign(modification=maps.modification.map(MODS))
save("supplementary_table2", maps.reset_index())

# Supplementary Table 3: CLIP libraries of the readers with the read counts of each processing step
clip = utils.read(utils.RESOURCES / "reader_clip/provenance.tsv")
clip = clip.assign(deduplicated_reads=clip.deduplicated_reads.astype("Int64"))  # libraries without UMIs have none
save("supplementary_table3", clip.drop(columns="tools"))

# Supplementary Table 4: modification load per context, layer and covariate set, cohorts combined (script 21)
load = result("21_modification_load", "combined_load_effects.tsv")
load = load.assign(layer=load.layer.map(LAYERS), modification=load.modification.map(MODS))
columns = ["context", "layer", "model", "modification", "n_units", "n_cohorts", "n_negative", "effect_pct"]
save("supplementary_table4", load[[*columns, "effect_pct_ci_low", "effect_pct_ci_high", "z", "p", "q"]])

# Supplementary Table 5: age effect of every regulator in every context and layer (script 30), with the aged cohorts
columns = ["context", "layer", "regulator", "n_cohorts", "effect_pct", "effect_pct_ci_low", "effect_pct_ci_high"]
save("supplementary_table5", screen.assign(layer=screen.layer.map(LAYERS))[[*columns, "p", "q", "cohorts"]])

# Supplementary Table 6: binding sources and perturbation datasets of the reader compendium (script 37) with the
# publication of each accession (PMID and DOI from PubMed through the GEO record; GSE171790, GSE266096 and GSE306850
# have no PubMed link in GEO and were matched by the paper's text, or its title and the GEO contributors; ENCODE: Van
# Nostrand 2020)
PAPERS = {  # accession: study, PMID, DOI
    "GSE104867": ("Li 2018", "29855337", "10.1186/s13059-018-1436-y"),
    "GSE106607": ("Shi 2018", "30401835", "10.1038/s41586-018-0666-1"),
    "GSE119801": ("Regué 2019", "30692269", "10.1128/MCB.00528-18"),
    "GSE133898": ("Heck 2020", "32169943", "10.1261/rna.073502.119"),
    "GSE146726": ("Okholm 2020", "33287884", "10.1186/s13073-020-00812-8"),
    "GSE148039": ("Lasman 2020", "32943573", "10.1101/gad.340695.120"),
    "GSE148536": ("Lasman 2020", "32943573", "10.1101/gad.340695.120"),
    "GSE153688": ("Yu 2021", "34643063", "10.1002/advs.202101329"),
    "GSE153689": ("Yu 2021", "34643063", "10.1002/advs.202101329"),
    "GSE158742": ("Dixit 2021", "33023892", "10.1158/2159-8290.CD-20-0331"),
    "GSE171790": ("Zhuang 2023", "36670199", "10.1038/s41380-023-01953-z"),
    "GSE188846": ("Ito-Kureha 2022", "35879451", "10.1038/s41590-022-01268-1"),
    "GSE192705": ("Alvarez-Dominguez 2022", "35036870", "10.1016/j.isci.2021.103680"),
    "GSE192792": ("Han 2022", "35992092", "10.1016/j.isci.2022.104823"),
    "GSE194207": ("Flamand 2022", "35438793", "10.1093/nar/gkac251"),
    "GSE199962": ("Niu 2022", "36039295", "10.1016/j.isci.2022.104908"),
    "GSE201070": ("Zou 2023", "36694229", "10.1186/s13059-023-02862-8"),
    "GSE214882": ("Zou 2023", "37949069", "10.1016/j.molcel.2023.10.028"),
    "GSE230846": ("Nabeel-Shah 2024", "39303720", "10.1016/j.molcel.2024.08.037"),
    "GSE264498": ("Zhang 2024", "39543097", "10.1038/s41467-024-54000-y"),
    "GSE266096": ("Zhang 2025", "39476834", "10.1016/j.stem.2024.10.002"),
    "GSE270417": ("Jang 2025", "39587835", "10.1002/advs.202410806"),
    "GSE292558": ("Schymik 2026", "41943844", "10.7150/ijbs.122142"),
    "GSE306850": ("Zhou 2026", "42697900", "10.1038/s41380-026-03840-9"),
    "GSE315406": ("Her 2026", "42025161", "10.1016/j.cels.2026.101588"),
    "GSE326611": ("Shi 2026", "41973927", "10.1073/pnas.2518250123"),
    "GSE49339": ("Wang 2014", "24284625", "10.1038/nature12730"),
    "GSE63591": ("Wang 2015", "26046440", "10.1016/j.cell.2015.05.014"),
    "GSE73845": ("Degrauwe 2016", "27184842", "10.1016/j.celrep.2016.04.086"),
    "GSE78507": ("Conway 2016", "27068461", "10.1016/j.celrep.2016.03.052"),
    "GSE83822": ("Preitner 2016", "27385015", "10.1242/dev.128348"),
    "GSE86214": ("Shi 2017", "28106072", "10.1038/cr.2017.15"),
    "GSE90639": ("Huang 2018", "29476152", "10.1038/s41556-018-0045-z"),
    "GSE90684": ("Huang 2018", "29476152", "10.1038/s41556-018-0045-z"),
    "ENCSR": ("Van Nostrand 2020", "32728246", "10.1038/s41586-020-2077-3"),
}
PAPERS["GSE147849"] = PAPERS["GSE148039"]
PAPERS["Yu 2021 Table S3"] = PAPERS["GSE153688"]
ACCESSION = r"GSE\d+|ENCSR|Yu 2021 Table S3"  # a GEO series, an ENCODE experiment or a table of the study


def publication(accession):
    first = re.search(ACCESSION, accession)[0]
    return pd.Series(PAPERS.get(first, ("", "", "")), index=["study", "PMID", "DOI"])


compendium = utils.OUT / "37_reader_compendium"
columns = ["reader", "source", "species", "cell", "tissue_class", "rank", "method", "accession", "T1", "gene_id"]
bound = utils.read(compendium / "bound_genes.tsv", usecols=[*columns, "note"])
binding = bound.groupby(["reader", "source"], sort=False).agg(
    species=("species", "first"),
    cell_type=("cell", "first"),
    tissue_class=("tissue_class", "first"),
    rank=("rank", "first"),
    method=("method", "first"),
    accession=("accession", "first"),
    n_genes=("gene_id", "nunique"),
    n_genes_with_m6A=("T1", "sum"),
    note=("note", "first"),
).reset_index().rename(columns={"source": "dataset"}).assign(kind="binding")
columns = ["dataset", "reader", "accession", "species", "cell_type", "direction", "readout", "gene_id", "padj", "note"]
trial = utils.read(compendium / "perturbation.tsv", usecols=columns)
perturbation = trial.groupby(["dataset", "reader", "readout"], sort=False).agg(
    species=("species", "first"),
    cell_type=("cell_type", "first"),
    method=("direction", "first"),
    accession=("accession", "first"),
    n_genes=("gene_id", "nunique"),
    n_genes_padj_below_0_05=("padj", lambda p: (p < 0.05).sum()),
    note=("note", "first"),
).reset_index().assign(kind="perturbation")
table6 = pd.concat([binding, perturbation], ignore_index=True)
table6 = pd.concat([table6, table6.accession.apply(publication)], axis=1)
columns = ["kind", "reader", "dataset", "species", "cell_type", "method", "readout", "tissue_class", "rank", "n_genes"]
columns += ["n_genes_with_m6A", "n_genes_padj_below_0_05", "accession", "study", "PMID", "DOI", "note"]
save("supplementary_table6", table6[columns])

# Supplementary Table 7: the 90 control RBPs of the network (ENCODE eCLIP proteins that do not read m6A, quantified in
# the proteome units) and the RBPs left out, with the reason
controls = result("33_regulator_network", "control_rbps.tsv").sort_values("selected", ascending=False, kind="stable")
save("supplementary_table7", controls[["rbp", "selected", "reason", "cell_lines", "n_units_quantified", "n_units"]])
