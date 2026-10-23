# Query website of the per-gene age effects: the static site of website/ with the data cut into one small JSON per gene
# (mouse genes join their one-to-one human ortholog), for GitHub Pages; every unit, m6A map and reader binding source
# carries its reference
# reads the units of Table S1, Table S2, data_sources.tsv, the outputs of scripts 09 and 10, the reader binding of
# script 33 (primary sources, rebuilt with network.binding_links) and the weights of the m6Brain clock of script 55
# the landscape page combines every gene over the cohorts of each context with the rule of the regulator screen
# (script 30)
# writes output/62_website/: the site, data/meta.json, data/genes.json (search index: key, human symbol, mouse symbol,
# mouse id of a gene joined to its human ortholog), data/genes/<key>.json and data/landscape/ (one JSON per context and
# layer, and landscape.tsv.gz with every column) and data/clock_<species>.json (columns gene_id, symbol and weight:
# the non-zero weights of the m6Brain clock fitted on all cohorts, and its intercept, to 7 significant digits, which the
# Clock page turns into predicted ages in the browser); the key of a gene is its human gene id, or its mouse id without
# a one-to-one human ortholog
import json
import re
import shutil

import numpy as np
import pandas as pd

import utils
from network import LOAD, READERS, atlas, binding_links, features, site_map, symbols, to_human

OUT = utils.OUT / "62_website"
LAYERS = {"protein_only": "P", "mrna_adjusted": "A", "mrna_only": "R"}  # single-nucleus mRNA becomes S below
TISSUE = {"takasugi2024": "cerebrum", "luo2024_bankadj": "frontal cortex", "matlik2024_neurons": "striatum",
          "hahn2023_msn": "caudate putamen", "hahn2023_neurons": "caudate putamen"}  # units whose name has no tissue
REGIONS = [("Cerebellum", r"cerebell"), ("Striatum", r"caudate|putamen|nac|accumbens|striatum"),
           ("Hippocampus", r"hippocampus|ca1"), ("Amygdala", r"amygdala"), ("Thalamus and hypothalamus", r"thalamus"),
           ("Midbrain, brainstem and spinal cord", r"nigra|tegmental|brainstem|pons|medulla|spinal"),
           ("Olfactory bulb", r"olfactory"), ("Cerebrum", r"cerebrum"),
           ("White matter and ventricles", r"white|callosum|choroid|subventricular"),
           ("Cortex", r"cortex|dlpfc|dpfc|ba9|ba24|cingulate|temporal|smtg|ipl")]
MAP_SOURCE = {"m6a_hippocampus_cyto": "sun2025_gse270643", "m6a_hippocampus_syn": "sun2025_gse270643",
              "m6a_brain_shafik2026": "shafik2026_natneurosci"}
LISTS = {  # binding lists of the reader network (network.LISTS; ref/target_sets of the data release)
    "Darnell2011_brain": ("Darnell JC et al. FMRP stalls ribosomal translocation on mRNAs linked to synaptic "
                          "function and autism. Cell 146, 247-261 (2011); HITS-CLIP of mouse brain polyribosomes, "
                          "842 high-confidence targets", "10.1016/j.cell.2011.06.013"),
    "Wu2019_PRRC2A_brain": ("Wu R et al. A novel m6A reader Prrc2a controls oligodendroglial specification and "
                            "myelination. Cell Res 29, 23-41 (2019); mouse brain RIP-seq targets that change in the "
                            "Nestin-Cre knockout brain", "10.1038/s41422-018-0113-8"),
    "Shi2026_striatum_CLIP": ("Shi Z et al. mRNA m6A modifications and the RNA-binding protein YTHDF1 affect "
                              "translational control in both normal and pathological learning. Proc Natl Acad Sci "
                              "USA 123, e2518250123 (2026); CLIP of mouse striatum", "10.1073/pnas.2518250123"),
    "Ascano2012": ("Ascano M et al. FMRP targets distinct mRNA sequence elements to regulate protein expression. "
                   "Nature 492, 382-386 (2012); PAR-CLIP of HEK293 cells", "10.1038/nature11737"),
    "Liu2015_HNRNPC": ("Liu N et al. N6-methyladenosine-dependent RNA structural switches regulate RNA-protein "
                       "interactions. Nature 518, 560-564 (2015); genes with high-confidence m6A switches, HEK293T "
                       "cells", "10.1038/nature14234"),
    "Liu2017_RBMX": ("Liu N et al. N6-methyladenosine alters RNA structure to regulate binding of a low-complexity "
                     "protein. Nucleic Acids Res 45, 6051-6063 (2017); genes with high-confidence HNRNPG-bound m6A "
                     "sites, HEK293T cells", "10.1093/nar/gkx141"),
}


def region(unit):
    name = TISSUE.get(unit.unit_id, unit.unit.lower())
    return next((group for group, pattern in REGIONS if re.search(pattern, name)), "Whole brain")


def sig(x):
    r = float(f"{x:.3g}")  # three significant digits, but a value never crosses the q cutoff, which classifies the gene
    return float(x) if (r < utils.Q) != (x < utils.Q) else r


effects = pd.concat([utils.read(utils.OUT / "10_protein_age_effects" / f"{layer}.tsv.gz").assign(layer=code)
                     for layer, code in LAYERS.items()]).query("status == 'ok'")
effects["symbol"] = effects.symbol.fillna(effects.gene_id)  # a gene without a symbol is named by its id
single_nucleus = effects.unit.str.endswith(("neurons", "msn"))
effects.loc[single_nucleus & (effects.layer == "R"), "layer"] = "S"

# references: the first entry of data_sources.tsv of each dataset, map and CLIP dataset; the binding lists
sources_tsv = utils.read(utils.ROOT / "data_sources.tsv", keep_default_na=False)
sources_tsv["dataset"] = sources_tsv.id.str.split("/").str[0]
first = sources_tsv.drop_duplicates("dataset").set_index("dataset")
refs = {k: dict(text=r.citation, doi=r.doi, accession=r.accession) for k, r in first.iterrows()}
refs |= {k: dict(text=text, doi=doi, accession="") for k, (text, doi) in LISTS.items()}

units = utils.read(utils.ROOT / "source_data" / "supplementary_table1_units.tsv", keep_default_na=False)
units["region"] = [region(u) for u in units.itertuples()]
units["ref"] = units.unit_id.map(utils.read(utils.PREPARED / "units.tsv", index_col="unit").dataset)
genes_per_layer = effects.groupby(["unit", "layer"]).size()
units["layers"] = [genes_per_layer[u].to_dict() for u in units.unit_id]
index = {u: i for i, u in enumerate(units.unit_id)}

maps = utils.read(utils.ROOT / "source_data" / "supplementary_table2.tsv").set_index("map")
load_maps = {"mouse": sorted(atlas["map"][atlas.species == "Mus musculus"].unique()), "human": [utils.BRAIN_M6A]}
load_maps = {s: [dict(source=maps.source[m], sites=int(maps.sites[m]), ref=MAP_SOURCE[m]) for m in ms]
             for s, ms in load_maps.items()}

# one key per gene: the human gene id, which a mouse gene takes over through its one-to-one ortholog
species = effects.gene_id.str.startswith("ENSMUSG").map({True: "mouse", False: "human"})
human_ids = set(effects.gene_id[species == "human"])
ortholog = effects.symbol.str.upper().map(to_human).map(symbols["human"])
effects["key"] = effects.gene_id.where((species == "human") | ~ortholog.isin(human_ids), ortholog)
key_of = effects.drop_duplicates("gene_id").set_index("gene_id").key

# reader binding of each gene in the primary source of every m6A reader: 1, m6A within 50 nt of the binding (lists: a
# listed gene with m6A); 0, bound without m6A
binding, reader_sources = {}, {}
for s in ("mouse", "human"):
    S, U, table = binding_links(s, site_map(s), [], {})
    primary = table[table.primary].set_index("reader")
    other = "human" if s == "mouse" else "mouse"
    reader_sources[s] = {
        r: dict(source=p.source, kind=p.kind, mapped=p.source.endswith("_mapped"),
                ref=re.sub(r"_FXR[12]$", "", p.source.removesuffix("_mapped")) if p.kind == "list"
                else f"clip_{other if p.kind == 'mapped_peaks' else s}_{r}")
        for r, p in primary.iterrows()}
    marked = pd.DataFrame({r: S[f"{r}|{primary.source[r]}"] > 0 for r in READERS})
    bound = marked | pd.DataFrame({r: U[f"{r}|{primary.source[r]}"] > 0 for r in READERS})
    binding[s] = {g: {r: int(marked.at[g, r]) for r in READERS if bound.at[g, r]}
                  for g in bound.index[bound.any(axis=1)]}

# landscape: unit estimates and standard errors averaged within a cohort, DerSimonian-Laird random effects over the
# cohorts, BH within the context and layer, as in script 30 for the regulators; the aged-cohort repeat of mouse
# protein (a sensitivity analysis) is left out
CODES = {"protein_only": ["P"], "mrna_adjusted": ["A"], "mrna_only": ["R", "S"]}
cohort_of, cohort_name = units.set_index("unit_id").ref, units.drop_duplicates("ref").set_index("ref").dataset


def combine(group):
    pooled = utils.random_effects(group.estimate, group.se)
    same = int((np.sign(group.estimate) == np.sign(pooled["est"])).sum())
    return pd.Series(pooled | dict(symbol=group.symbol.iat[0], same=same))


landscape, contexts = [], []
for context, layer, listed in utils.CONTEXTS:
    rows = effects[effects.unit.isin(listed.split()) & effects.layer.isin(CODES[layer])]
    rows = rows.assign(cohort=rows.unit.map(cohort_of))
    cohorts = rows.groupby(["gene_id", "cohort"]).agg(symbol=("symbol", "first"), estimate=("estimate", "mean"),
                                                      se=("se", "mean")).reset_index()
    table = cohorts.groupby("gene_id")[["symbol", "estimate", "se"]].apply(combine).reset_index()
    s, code = context.split("_")[0], LAYERS[layer]
    table = table.assign(context=context, layer=code, q=utils.bh(table.p), key=table.gene_id.map(key_of),
                         m6a=table.gene_id.map(features[s][LOAD[s]]))
    for column, source in [("effect_pct", "est"), ("effect_pct_ci_low", "ci_low"), ("effect_pct_ci_high", "ci_high")]:
        table[column] = utils.pct(table[source])
    landscape.append(table.sort_values("p"))
    hits = table[table.q < utils.Q]
    contexts.append(dict(id=f"{context}_{code}", species=s, genes=len(table),
                         cohorts=sorted({cohort_name[c] for c in rows.cohort}), down=int((hits.est < 0).sum()),
                         up=int((hits.est > 0).sum())))
    print(f"landscape {context} {code}: {len(table)} genes, {contexts[-1]['down']} down and {contexts[-1]['up']} up")

used = set(units.ref) | {m["ref"] for ms in load_maps.values() for m in ms}
used |= {x["ref"] for rs in reader_sources.values() for x in rs.values()}
meta = dict(units=units.drop(columns=["contexts", "samples", "permutation_blocks"]).to_dict("records"), maps=load_maps,
            readers=READERS, landscape=contexts, sources=reader_sources, refs={k: refs[k] for k in sorted(used)})

shutil.rmtree(OUT, ignore_errors=True)
shutil.copytree(utils.ROOT / "website", OUT)
(OUT / "data" / "genes").mkdir(parents=True)
(OUT / "data" / "meta.json").write_text(json.dumps(meta, separators=(",", ":"), ensure_ascii=False))
search = []
for key, rows in effects.groupby("key", sort=False):
    gene = {}
    for (gene_id, symbol), _ in rows.groupby(["gene_id", "symbol"]):
        s = "mouse" if gene_id.startswith("ENSMUSG") else "human"
        m6a = features[s][LOAD[s]].get(gene_id)
        gene[s] = dict(id=gene_id, symbol=symbol, m6a=None if pd.isna(m6a) else int(m6a),
                       readers=binding[s].get(gene_id, {}))
    gene["rows"] = [[index[r.unit], r.layer, round(r.effect_pct, 2), round(r.effect_pct_ci_low, 2),
                     round(r.effect_pct_ci_high, 2), sig(r.p), sig(r.q), int(r.n)] for r in rows.itertuples()]
    (OUT / "data" / "genes" / f"{key}.json").write_text(json.dumps(gene, separators=(",", ":"), allow_nan=False))
    human, mouse = gene.get("human", {}), gene.get("mouse", {})
    search.append([key, human.get("symbol", ""), mouse.get("symbol", ""), mouse.get("id", "") if human else ""])
search.sort(key=lambda g: (g[1] or g[2]).upper())
(OUT / "data" / "landscape").mkdir()
for table, entry in zip(landscape, contexts, strict=True):
    rows = [[r.key, r.symbol, round(r.effect_pct, 2), round(r.effect_pct_ci_low, 2), round(r.effect_pct_ci_high, 2),
             sig(r.p), sig(r.q), int(r.k), int(r.same), None if np.isnan(r.i2) else round(r.i2),
             None if np.isnan(r.m6a) else int(r.m6a)] for r in table.itertuples()]
    path = OUT / "data" / "landscape" / f"{entry['id']}.json"
    path.write_text(json.dumps(rows, separators=(",", ":"), allow_nan=False))
names = dict(est="estimate", k="n_cohorts", same="n_same_direction", m6a="m6a_sites")
columns = ["context", "layer", "gene_id", "symbol", "n_cohorts", "n_same_direction", "estimate", "se", "ci_low",
           "ci_high", "effect_pct", "effect_pct_ci_low", "effect_pct_ci_high", "p", "q", "tau2", "i2", "m6a_sites"]
utils.write(pd.concat(landscape).rename(columns=names)[columns], OUT / "data" / "landscape" / "landscape.tsv.gz")
(OUT / "data" / "genes.json").write_text(json.dumps(search, separators=(",", ":")))

# m6Brain clock fitted on all cohorts (script 55: non-zero genes and the intercept); symbols from script 52
for s in ("mouse", "human"):
    file = utils.OUT / "55_m6brain_clock" / "weights" / f"{s}_protein__all__m6brain_clock.tsv"
    clock = utils.read(file, index_col="gene_id").weight
    intercept = clock.pop("(intercept)")
    symbol = pd.read_parquet(utils.OUT / "52_proteome_clock" / "runs" / f"{s}_protein__all.parquet").symbol[clock.index]
    clock = dict(gene_id=clock.index.tolist(), symbol=symbol.tolist(), weight=[float(f"{w:.7g}") for w in clock],
                 intercept=float(f"{intercept:.7g}"))
    (OUT / "data" / f"clock_{s}.json").write_text(json.dumps(clock, separators=(",", ":"), allow_nan=False))
print(f"{len(search)} genes, {len(effects)} age effects in {len(units)} units, {len(used)} references")
