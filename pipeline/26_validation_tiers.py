# Validation tiers of the m6A discovery results: a test passes when the validation estimate has the discovery sign and
# a one-sided P below 0.05. mouse_to_human: brain-map m6A load, m6A signature and YTHDF3 level in each human context
# against the matching mouse context; within_human: BrainGVEX (luo2024_bankadj) against the other human proteomes
# reads outputs of scripts 20, 21 and 30; writes output/26_validation_tiers/
import numpy as np
import pandas as pd
import scipy.stats as sps

import utils

OUT = utils.OUT / "26_validation_tiers"
P_ONE_SIDED = 0.05
DISCOVERY_UNIT = "luo2024_bankadj"
# human context and layer, mouse context and layer of the discovery
PAIRS = [
    ("human_protein", "protein_only", "mouse_protein", "protein_only"),
    ("human_striatum", "mrna_only", "mouse_mrna", "mrna_only"),
    ("human_cortex", "mrna_only", "mouse_mrna", "mrna_only"),
    ("human_hippocampus", "mrna_only", "mouse_mrna", "mrna_only"),
    ("human_cerebellum", "mrna_only", "mouse_mrna", "mrna_only"),
    ("human_neurons", "mrna_only", "mouse_neurons", "mrna_only"),
]

load = utils.read(utils.OUT / "21_modification_load" / "combined_load_effects.tsv")
load = load[(load.modification == "m6A") & (load.model == "without protein lifetime")]
signature = utils.read(utils.OUT / "20_modification_signatures" / "combined_signatures.tsv")
signature = signature[signature.modification == "m6A"]
regulator = utils.read(utils.OUT / "30_regulator_screen" / "regulator_effects.tsv").query("regulator == 'YTHDF3'")
regulator = regulator.assign(z=regulator.estimate / regulator.se)
# test: table, estimate column
TESTS = {
    "m6A load (brain map)": (load, "effect_pct"),
    "m6A signature": (signature, "diff_pct"),
    "YTHDF3 level": (regulator, "effect_pct"),
}


def verdict(direction, z):
    # one-sided P in the discovery direction
    p = sps.norm.sf(direction * z)
    same = bool(np.sign(z) == direction)
    return dict(p_one_sided=p, same_sign=same, passed=same and p < P_ONE_SIDED)


def get(table, context, layer):
    return table[(table.context == context) & (table.layer == layer)].iloc[0]


rows = []
for test, (table, estimate) in TESTS.items():
    for context, layer, mouse_context, mouse_layer in PAIRS:
        mouse, human = get(table, mouse_context, mouse_layer), get(table, context, layer)
        row = dict(tier="mouse_to_human", test=test, discovery=mouse_context, validation=context, layer=layer)
        row |= dict(n_discovery=mouse.n_cohorts, n_validation=human.n_cohorts)
        row |= dict(discovery_pct=mouse[estimate], discovery_p=mouse.p, validation_pct=human[estimate])
        rows.append(row | verdict(np.sign(mouse.z), human.z) | dict(validation_z=human.z))

# within human: Luo 2024 against the other protein cohorts, each cohort the mean of its units
units = utils.read(utils.OUT / "21_modification_load" / "modification_load_effects.tsv")
units = units[(units.species == "human") & (units.layer == "protein_only") & (units.modification == "m6A")]
units = units[units.model == "without protein lifetime"]
discovery = units[units.unit == DISCOVERY_UNIT].iloc[0]
cohorts = units[units.unit != DISCOVERY_UNIT].groupby("cohort")[["estimate", "se"]].mean()
validation = utils.random_effects(cohorts.estimate, cohorts.se)
row = dict(tier="within_human", test="m6A load (brain map)", discovery=DISCOVERY_UNIT, validation="other cohorts")
row |= dict(layer="protein_only", n_discovery=1, n_validation=len(cohorts), discovery_pct=discovery.effect_pct)
row |= dict(discovery_p=discovery.p, validation_pct=utils.pct(validation["est"]), validation_z=validation["z"])
row |= dict(validation_i2=validation["i2"])
rows.append(row | verdict(np.sign(discovery.estimate), validation["z"]))

tiers = pd.DataFrame(rows)
utils.write(tiers, OUT / "validation_tiers.tsv")
shown = ["tier", "test", "validation", "n_validation", "discovery_pct", "validation_pct", "p_one_sided", "passed"]
print(tiers[shown].round(4).to_string(index=False))
