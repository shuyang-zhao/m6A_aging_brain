# Proteome clock in disease: the m6Brain clock (script 55) in the cases and controls of AD (Ping 2020, Merrihew 2023),
# PD (Blumenreich 2024) and ALS (Umoh 2018) cohorts and along the CAG length of the HD knock-in allelic series
# (Langfelder 2016 striatum, the 6-month proteome of the article). Predicted age = intercept + the proteins standardised
# in the region (cases and controls together, missing = 0) @ weights, years (HD: months), the mean over the regions of
# a donor; age acceleration = residual of the predicted on the chronological age in the controls (HD: deviation from
# the wild type, all of one age). OLS of the acceleration on disease (HD: CAG length) with age, sex and batch (HD: sex
# and batch, one age) = years (HD: months per CAG) of acceleration; the same model on the individual m6A index of
# script 54 (SD; utils.load_design). The clock is linear and the control fit adds no gene, so the
# disease effect splits exactly over gene groups: the part of a group is its weighted sum of the standardised proteins,
# the parts of a partition add up to the predicted age less the intercept and their effects to the effect on the
# acceleration; partitions by cell source (ubiquitous, else the cell type of largest share) and by m6A sites (the
# classes of Figure 3a), their cross and cell source by the sign of the weight (proteins that decline or rise with
# age); per donor, per region in the multi-region cohorts, and for AD combined with the random-effects weights of the
# total; n_genes and weight_share (share of sum |w|) over the genes in the score. Comparators (comparators.tsv): the
# 27 published clocks of utils.clock_scores, the unguided elastic net (script 55) and the label-free score scored in
# every study, the score standardised in the region (donor mean over regions) and the same OLS, effects in SD of the
# score (HD: per CAG), ranked by effect; the m6Brain clock in SD and in years. Severity (severity.tsv): the same models
# on the Braak stage and CERAD score of the AD cohorts (ordinal as numeric, per stage; cases and controls, with
# age, sex and batch); the PD, ALS and HD sample tables hold no severity measure. Composition: the clock effect of AD
# and PD with the neuron index of the region (mean z of SNAP25, SYT1, STMN2, NEFL, NEFM and SYN1) added.
# reads the prepared dataset folders (01), data/raw/umoh2018_emmm201708202, script 09, the weights of 50 and 52
# (label-free) and 55 (clock: the leave-cohort-out model of the study; Umoh: the model without Ping, Langfelder: all
# mouse cohorts);
# writes output/53_disease_clock/ (tests, comparators, severity, scores, sample_scores, decomposition and
# decomposition_regions)
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats as sps

import utils

OUT = utils.OUT / "53_disease_clock"
CLOCK = utils.OUT / "52_proteome_clock" / "runs"
WEIGHTS = utils.OUT / "55_m6brain_clock" / "weights"
UMOH = utils.DATA / "raw" / "umoh2018_emmm201708202"
# study: disease, dataset, tissue column values, case and control queries; weights of the leave-cohort-out model of
# the study (Umoh and Langfelder are in no training set: the models of all human and of all mouse cohorts; for the
# clock, Umoh takes the model fitted without Ping, with which it shares Emory donors)
STUDIES = {
    "ping2020": ("AD", "ping2020_pxd020296", ["DLPFC"], "diagnosis == 'AD'", "diagnosis == 'control'"),
    "merrihew2023": (
        "AD",
        "merrihew2023_pxd034525",
        ["Superior and middle temporal gyri (SMTG)", "Hippocampus", "Caudate", "Inferior parietal lobule (IPL)"],
        "diagnosis == 'Sporadic ADD'",
        "diagnosis == 'HCF/Low ADNC'",
    ),
    "blumenreich2024": (
        "PD",
        "blumenreich2024_pxd047134",
        ["Striatum", "Substantia nigra", "Cingulate gyrus", "Middle temporal gyrus", "Occipital cortex"],
        "diagnosis != 'Control'",
        "diagnosis == 'Control'",
    ),
    "langfelder2016": ("HD", "langfelder2016_pxd006302", None, None, None),
    "umoh2018": ("ALS", None, None, None, None),
}
SEVERITY = {  # measure: column of the sample table
    "ping2020": dict(braak="braak", cerad="cerad"),
    "merrihew2023": dict(braak="Braak Stage", cerad="CERAD Score"),
}
NEURON = "SNAP25 SYT1 STMN2 NEFL NEFM SYN1".split()
OURS = {
    "m6Brain clock": "m6Brain clock", "unguided elastic net": "elastic net", "label-free m6Brain score": "label-free"
}
OUTCOMES = dict(clock="m6Brain clock", m6a_index="m6A index")  # outcome: mode
SOURCES = dict(
    ct_frac_neuron="neuron", ct_frac_astro="astrocyte", ct_frac_oligo="oligodendrocyte", ct_frac_micro="microglia"
)
UBIQUITOUS = 0.1  # at least this share in each of the four cell types (script 21)
SITE_BINS, SITE_CLASSES = [-1, 0, 1, 3, np.inf], ["0", "1", "2-3", ">=4"]  # m6A sites as in Figure 3a (script 21)


def gene_matrix(matrix, labels, species):
    # median-normalised on proteins seen in 70% of samples, one protein group per gene (most often seen, then most
    # abundant), genes seen in half of the samples
    detected = matrix.notna().mean(axis=1) >= 0.7
    matrix = matrix.sub(matrix[detected].sub(matrix[detected].median(axis=1), axis=0).median(), axis=1)
    rank = pd.DataFrame(dict(symbol=labels.gene_symbol.reindex(matrix.index), seen=matrix.notna().mean(axis=1)))
    rank = rank.assign(mean=matrix.mean(axis=1)).dropna(subset=["symbol"])
    rank = rank[~rank.symbol.str.contains(";")].sort_values(["seen", "mean"], ascending=False).drop_duplicates("symbol")
    rank["gene_id"] = rank.symbol.map(utils.symbol_steps(species)["to_ensembl"])
    rank = rank.dropna(subset=["gene_id"]).drop_duplicates("gene_id")
    matrix = matrix.loc[rank.index].set_axis(rank.gene_id)
    return matrix[matrix.notna().mean(axis=1) >= 0.5]


def cohort(study):
    # (matrix, samples) per region; samples: donor_id, disease (HD: CAG length), age, sex, batch
    disease, dataset, tissues, case, control = STUDIES[study]
    if study == "umoh2018":
        # log2 LFQ from which the authors regressed age, sex and PMI
        lfq = pd.read_excel(UMOH / "emmm201708202-sup-0003-DatasetEV2.xlsb", sheet_name="Dataset EV2", header=2)
        lfq = lfq.set_index("Unique ID").filter(regex=r"^\d+_mfon_\d+$").astype(float)
        info = pd.read_excel(UMOH / "emmm201708202-sup-0002-DatasetEV1.xlsb", sheet_name="Dataset EV1", header=1)
        info = info.set_index("RAW Filename").loc[lfq.columns].query("Grouping in ['ALS', 'Control']")
        labels = pd.DataFrame(dict(gene_symbol=lfq.index.str.split("|").str[0]), index=lfq.index)
        samples = pd.DataFrame(
            dict(donor_id=info["Case ID #"].astype(str), disease=(info.Grouping == "ALS").astype(int)),
            index=info.index,
        ).assign(age=info["(years).1"], sex=info.Sex.str.lower(), batch="one")
        return [(gene_matrix(lfq[samples.index], labels, "human"), samples)]
    folder = utils.PREPARED / dataset
    table = utils.read(folder / "protein_samples.tsv", index_col="sample_id", dtype={"sample_id": str})
    protein, labels = utils.read(folder / "protein.parquet"), utils.read(folder / "protein_features.tsv", index_col=0)
    if study == "langfelder2016":
        # the 6-month proteome of the article; CAG: median measured repeat length of the genotype, wild type 7 (Htt)
        table = table[table.age_months == 6]
        cag = table.groupby("genotype")["Seq. CAG Repeat #"].median().dropna().to_dict() | {"WT": 7}
        samples = table[["donor_id", "sex", "batch"]].assign(disease=table.genotype.map(cag), age=table.age_months)
        return [(gene_matrix(protein[samples.index], labels, "mouse"), samples)]
    excluded = table.author_excluded.eq(True) if "author_excluded" in table else False
    out = []
    for tissue in tissues:
        rows = table[(table.tissue == tissue) & ~excluded]
        rows = rows[rows.eval(case) | rows.eval(control)]
        samples = rows[["donor_id", "sex"]].assign(disease=rows.eval(case).astype(int), age=rows.age_years)
        samples = samples.assign(**{m: rows[c] for m, c in SEVERITY.get(study, {}).items()})
        samples["region"] = tissue
        samples["batch"] = rows.batch
        out.append((gene_matrix(protein[samples.index], labels, "human"), samples))
    return out


def weights(study):
    if study == "umoh2018":
        runs = utils.read(utils.OUT / "50_transcript_ageing_model" / "runs" / "human_protein.parquet")
        return runs[["tabpfn", "without m6A sites and binding"]].set_axis(["model", "without_m6a"], axis=1)
    name = "mouse_protein__all" if STUDIES[study][0] == "HD" else f"human_protein__{STUDIES[study][1]}"
    return utils.read(CLOCK / f"{name}.parquet")[["model", "without_m6a"]]


def clock(study, variant="m6brain_clock"):
    disease, dataset = STUDIES[study][:2]
    name = "mouse_protein__all" if disease == "HD" else f"human_protein__{dataset or STUDIES['ping2020'][1]}"
    w = utils.read(WEIGHTS / f"{name}__{variant}.tsv", index_col="gene_id").weight
    return w.drop("(intercept)"), w["(intercept)"]


def standardised(x):
    return (x - x.mean()) / x.std()


def partitions(w, table):
    genes = w.index
    shares = table.loc[genes, list(SOURCES)].rename(columns=SOURCES).dropna()
    source = shares.idxmax(axis=1).where(shares.min(axis=1) < UBIQUITOUS, "ubiquitous")
    source = source.reindex(genes).fillna("no cell-type data")
    sites = pd.cut(table.loc[genes, "n_sites_m6A"], SITE_BINS, labels=SITE_CLASSES).astype(str)
    direction = pd.Series(np.where(w < 0, "declining", "rising"), index=genes)
    return pd.DataFrame(
        {
            "cell source": source,
            "m6A sites": sites,
            "cell source x m6A sites": source + " | " + sites,
            "cell source x direction": source + " | " + direction,
        }
    )


def neuron_index(matrix, species):
    symbol = utils.ensg_to_symbol(species).str.upper()
    rows = matrix.loc[matrix.index.intersection(symbol.index[symbol.isin(NEURON)])]
    rows = rows[rows.notna().mean(axis=1) >= 0.8]
    return rows.sub(rows.mean(axis=1), axis=0).div(rows.std(axis=1), axis=0).mean().rename("neuron")


def disease_effect(table, y, covariates, term="disease"):
    fit = smf.ols(f"y ~ {term}{covariates}", table.assign(y=y)).fit()
    low, high = fit.conf_int().loc[term]
    est, se, p = fit.params[term], fit.bse[term], fit.pvalues[term]
    return dict(n=int(fit.nobs), est=est, se=se, ci_low=low, ci_high=high, p=p)


def decompose(table, total, parts, labels, mass, covariates):
    # disease effect of the total score and of the part of every group; labels and mass |w| of the genes in the score
    rows = [dict(partition="total", group="all genes", n_genes=len(labels), weight_share=1.0)]
    rows[0] |= disease_effect(table, total, covariates)
    for partition in labels:
        for group, (n, m) in mass.groupby(labels[partition]).agg(["size", "sum"]).iterrows():
            row = dict(partition=partition, group=group, n_genes=int(n), weight_share=m / mass.sum())
            rows.append(row | disease_effect(table, parts[(partition, group)], covariates))
    rows = pd.DataFrame(rows)
    return rows.assign(share_of_total=rows.est / rows.est.iloc[0])


features = {s: utils.read(utils.OUT / "09_transcript_features" / f"{s}_features.tsv.gz").set_index("gene_id") for s in
            ("mouse", "human")}
scores, samples_out, tests, decomposition, by_region, comparators, severity = [], [], [], [], [], [], []
for study, (disease, *_) in STUDIES.items():
    mouse = disease == "HD"
    species = "mouse" if mouse else "human"
    w, (wc, intercept), (we, intercept_en) = weights(study), clock(study), clock(study, "elastic_net")
    groups = partitions(wc, features[species])
    regions, partials, genes, compared = [], [], [], []
    for matrix, samples in cohort(study):
        # cases and controls standardised together, as one unit of scripts 52 and 54
        z = matrix.reindex(w.index).T.pipe(standardised).fillna(0)
        zc = matrix.reindex(wc.index).T.pipe(standardised).fillna(0)
        _, x, y = utils.load_design(matrix, features[species])
        index = pd.Series(np.linalg.pinv(x)[1] @ y, index=matrix.columns, name="m6a_index").pipe(standardised)
        predicted = (intercept + zc @ wc).rename("predicted_age")
        neuron = neuron_index(matrix, species)
        regions.append(samples.join((z @ w).pipe(standardised)).join(index).join(predicted).join(neuron))
        elastic = intercept_en + matrix.reindex(we.index).T.pipe(standardised).fillna(0) @ we
        compared.append(
            utils.clock_scores(matrix, species)
            .assign(**dict(zip(OURS, [predicted, elastic, (z @ w).model], strict=True)))
            .pipe(standardised)
        )
        quantified = matrix.index.intersection(wc.index)  # the other genes have z = 0
        contribution = zc[quantified] * wc[quantified]
        by_group = {p: contribution.T.groupby(groups.loc[quantified, p]).sum().T for p in groups}
        partials.append(pd.concat(by_group, axis=1))
        genes.append(quantified)
    first = ["disease", "age", "sex", "batch", *SEVERITY.get(study, {})]
    donors = pd.concat(regions).groupby("donor_id").agg(
        m6a_index=("m6a_index", "mean"), predicted_age=("predicted_age", "mean"), neuron=("neuron", "mean"),
        **{c: (c, "first") for c in first},
    )
    controls = donors[donors.disease == donors.disease.min()]  # HD: wild type
    control_fit = smf.ols("predicted_age ~ 1" if mouse else "predicted_age ~ age", controls).fit()
    for table in (donors, *regions):  # the fit is linear: the acceleration of a donor is the mean of its regions
        table["clock"] = table.predicted_age - control_fit.predict(table)
    per_sample = pd.concat(regions).assign(study=study, disease_name=disease)
    samples_out.append(per_sample.rename_axis("sample_id").reset_index())
    scores.append(donors.assign(study=study, disease_name=disease).reset_index())
    age = "" if mouse else " + age"
    covariates = f"{age} + sex" + (" + C(batch)" if donors.batch.nunique() > 1 else "")
    for outcome, mode in OUTCOMES.items():
        fit = disease_effect(donors, donors[outcome], covariates)
        unit = ("months" if mouse else "years") if outcome == "clock" else "SD"
        tests.append(dict(study=study, disease=disease, outcome=outcome, mode=mode, unit=unit) | fit)
    if disease in ("AD", "PD"):
        fit = disease_effect(donors, donors.clock, covariates + " + neuron")
        row = dict(study=study, disease=disease, outcome="clock_neuron", mode=OUTCOMES["clock"], unit="years")
        tests.append(row | fit)
    # comparators: scores standardised in the region, donor mean over the regions, in SD of the score
    scored = pd.concat(compared).groupby(pd.concat(regions).donor_id.to_numpy()).mean()
    table = donors.join(scored)
    for name in scored:
        fit = disease_effect(table, scored[name], covariates)
        group = OURS.get(name, "published")
        comparators.append(dict(study=study, disease=disease, clock=name, group=group, unit="SD") | fit)
    for measure in SEVERITY.get(study, {}):
        d = table.dropna(subset=[measure])
        for name, y, unit in [("m6Brain clock", d.clock, "years"), *[(c, d[c], "SD") for c in scored]]:
            fit = disease_effect(d, y, covariates, measure)
            severity.append(dict(study=study, disease=disease, measure=measure, clock=name, unit=unit) | fit)
    # a group absent from a region adds nothing to its predicted age there (missing = 0)
    parts = pd.concat(partials).fillna(0).groupby(pd.concat(regions).donor_id.to_numpy()).mean()
    union, mass = genes[0].append(genes[1:]).unique(), wc.abs()
    fit = decompose(donors, donors.clock, parts, groups.loc[union], mass[union], covariates)
    decomposition.append(fit.assign(study=study, disease=disease))
    if len(regions) > 1:
        for table, parts, quantified in zip(regions, partials, genes, strict=True):
            covariates = f"{age} + sex" + (" + C(batch)" if table.batch.nunique() > 1 else "")
            fit = decompose(table, table.clock, parts, groups.loc[quantified], mass[quantified], covariates)
            by_region.append(fit.assign(study=study, disease=disease, region=table.region.iloc[0]))

tests = pd.DataFrame(tests)
# AD: the two studies combined by random effects
for outcome, rows in tests[tests.disease == "AD"].groupby("outcome", sort=False):
    fit = utils.random_effects(rows.est, rows.se)
    tests.loc[len(tests)] = dict(study="AD combined", disease="AD", outcome=outcome, mode=rows["mode"].iloc[0]) | dict(
        unit=rows.unit.iloc[0], n=rows.n.sum(), est=fit["est"]
    ) | {k: fit[k] for k in ("se", "ci_low", "ci_high", "p")}
# comparators: AD combined by random effects, the rank by effect among the clocks of a study (1 = largest acceleration),
# the m6Brain clock in years added
comparators = pd.DataFrame(comparators)
for name, rows in comparators[comparators.disease == "AD"].groupby("clock", sort=False):
    fit = utils.random_effects(rows.est, rows.se)
    row = dict(study="AD combined", disease="AD", clock=name, group=rows.group.iloc[0], unit="SD", n=rows.n.sum())
    comparators.loc[len(comparators)] = row | {k: fit[k] for k in ("est", "se", "ci_low", "ci_high", "p")}
comparators["rank"] = comparators.groupby("study").est.rank(ascending=False)
years = tests.query("outcome == 'clock'").assign(clock="m6Brain clock", group="m6Brain clock")
comparators = pd.concat([comparators, years[comparators.columns.drop("rank")]])
severity = pd.DataFrame(severity)
severity["rank"] = severity.query("unit == 'SD'").groupby(["study", "measure"]).est.rank(ascending=False)
# AD combined: every group with the random-effects weights of the total, so that the groups add up to its total
total = tests[(tests.disease == "AD") & (tests.outcome == "clock") & (tests.study != "AD combined")].set_index("study")
omega = 1 / (total.se**2 + utils.random_effects(total.est, total.se)["tau2"])
ad = pd.concat(decomposition).query("disease == 'AD'").assign(omega=lambda d: d.study.map(omega))
ad = ad.assign(est=ad.omega * ad.est, weight_share=ad.omega * ad.weight_share, se=ad.omega**2 * ad.se**2)
sums = ad.groupby(["partition", "group"], sort=False)[["omega", "est", "weight_share", "se"]].sum()
combined = pd.DataFrame(dict(weight_share=sums.weight_share, est=sums.est, se=np.sqrt(sums.se))).div(sums.omega, axis=0)
combined["ci_low"], combined["ci_high"] = combined.est - utils.Z95 * combined.se, combined.est + utils.Z95 * combined.se
combined["p"] = 2 * sps.norm.sf(abs(combined.est / combined.se))
combined["share_of_total"] = combined.est / combined.est["total", "all genes"]
decomposition.append(combined.reset_index().assign(study="AD combined", disease="AD"))
columns = ["study", "disease", "partition", "group", "n_genes", "weight_share", "est", "se", "ci_low", "ci_high", "p",
           "share_of_total"]
utils.write(pd.concat(decomposition)[columns].astype({"n_genes": "Int64"}), OUT / "decomposition.tsv")
utils.write(pd.concat(by_region)[[*columns[:2], "region", "n", *columns[2:]]], OUT / "decomposition_regions.tsv")
utils.write(pd.concat(scores), OUT / "scores.tsv.gz")
utils.write(pd.concat(samples_out), OUT / "sample_scores.tsv.gz")  # per sample and region
utils.write(tests, OUT / "tests.tsv")
utils.write(comparators, OUT / "comparators.tsv")
utils.write(severity, OUT / "severity.tsv")
print(tests.round(4).to_string(index=False))
print(comparators.query("clock == 'm6Brain clock'").round(4).to_string(index=False))
print(severity.query("clock == 'm6Brain clock'").round(4).to_string(index=False))
