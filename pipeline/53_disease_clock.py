# Proteome clock in disease: the score of script 52 in the cases and controls of AD (Ping 2020, Merrihew 2023), PD
# (Blumenreich 2024) and ALS (Umoh 2018) cohorts and along the CAG length of the HD knock-in allelic series (Langfelder
# 2016 striatum, the 6-month proteome of the article); OLS of the donor score on disease (HD: CAG length) with age, sex
# and batch (HD: sex and batch, one age), and the same model on the difference between the scores with and without the
# m6A features
# reads the prepared dataset folders (01), data/raw/umoh2018_emmm201708202 and the weights of scripts 50 and 52;
# writes output/53_disease_clock/
import pandas as pd
import statsmodels.formula.api as smf

import utils

OUT = utils.OUT / "53_disease_clock"
CLOCK = utils.OUT / "52_proteome_clock" / "runs"
UMOH = utils.DATA / "raw" / "umoh2018_emmm201708202"
# study: disease, dataset, tissue column values, case and control queries; weights of the leave-cohort-out model of
# the study (Umoh and Langfelder are in no training set: the models of all human and of all mouse cohorts)
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


def standardised(x):
    return (x - x.mean()) / x.std()


scores, samples_out, tests = [], [], []
for study, (disease, *_) in STUDIES.items():
    w = weights(study)
    regions = []
    for matrix, samples in cohort(study):
        # cases and controls standardised together, as one unit of script 52
        z = matrix.reindex(w.index).T.pipe(standardised).fillna(0)
        regions.append(samples.join((z @ w).pipe(standardised)))
    per_sample = pd.concat(regions).assign(study=study, disease_name=disease)
    samples_out.append(per_sample.rename_axis("sample_id").reset_index())
    first = ["disease", "age", "sex", "batch"]
    donors = pd.concat(regions).groupby("donor_id").agg(
        model=("model", "mean"), without_m6a=("without_m6a", "mean"), **{c: (c, "first") for c in first}
    )
    donors["m6a_part"] = donors.model - donors.without_m6a
    scores.append(donors.assign(study=study, disease_name=disease).reset_index())
    age = "" if disease == "HD" else " + age"
    batch = " + C(batch)" if donors.batch.nunique() > 1 else ""
    for outcome in ("model", "without_m6a", "m6a_part"):
        fit = smf.ols(f"{outcome} ~ disease{age} + sex{batch}", donors).fit()
        low, high = fit.conf_int().loc["disease"]
        row = dict(study=study, disease=disease, outcome=outcome, n=int(fit.nobs), est=fit.params["disease"])
        tests.append(row | dict(se=fit.bse["disease"], ci_low=low, ci_high=high, p=fit.pvalues["disease"]))

tests = pd.DataFrame(tests)
# AD: the two studies combined by random effects
for outcome, rows in tests[tests.disease == "AD"].groupby("outcome", sort=False):
    fit = utils.random_effects(rows.est, rows.se)
    tests.loc[len(tests)] = dict(study="AD combined", disease="AD", outcome=outcome, n=rows.n.sum(), est=fit["est"]) | {
        k: fit[k] for k in ("se", "ci_low", "ci_high", "p")
    }
utils.write(pd.concat(scores), OUT / "scores.tsv.gz")
utils.write(pd.concat(samples_out), OUT / "sample_scores.tsv.gz")  # per sample and region
utils.write(tests, OUT / "tests.tsv")
print(tests.round(4).to_string(index=False))
