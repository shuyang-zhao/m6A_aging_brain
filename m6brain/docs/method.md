# Method

## Label-free weights

For each species, TabPFN-3.5 (a tabular foundation model) was fitted to the combined age z of the genes in public brain proteome aging cohorts (Stouffer combination, weighted by √n, of the per-cohort age effects of each protein) using the transcript features of the gene: Orthrus embedding, RiboNN predictions, m6A sites per region (mouse: the two GLORI hippocampus maps; human: brain m6A-seq peaks) and reader binding, transcript structure (lengths, exons, GC, junction density), expression, and the cell-type origin of the transcript (`pipeline/50_transcript_ageing_model.py`). Each gene is predicted by a model fitted on genes of other chromosome groups (five groups), and this out-of-fold prediction is the label-free weight of the gene (`pipeline/52_proteome_clock.py`, run `all`). It is the column `prediction` of `src/m6brain/data/{mouse,human}_weights.tsv`; `observed_z` is the target and `m6A_sites` the site count of the gene.

## m6Brain clock

The clock is an elastic net of chronological age (months in mouse, years in human) on the proteins standardized within their unit (missing = 0), fitted on the stacked samples of all training units (`pipeline/55_m6brain_clock.py`). The penalty factor of protein j is p_j = exp(−γ s_j), where s_j is the percentile rank (0–1) of the absolute label-free weight of the gene and p is normalized to mean 1; the protein is divided by p_j before the fit and the coefficients are divided back, which is exact for the L1 term. The penalty is therefore lighter for the proteins whose age effect m6Brain predicts to be large, in either direction. The l1 ratio (0.1, 0.5, 0.9), the penalty strength α (30 values on a path) and γ (0.5, 1, 2, 4) minimize the mean squared error of an inner cross-validation that leaves out one training cohort at a time (5-fold by donor within a cohort); the model is then refitted on all training samples. The fit shipped here uses all training cohorts and keeps 279 proteins in mouse and 421 in human (l1 ratio 0.5, γ = 1 in both). The label-free weights that set the penalty are the ones of the model fitted on all cohorts.

`build_weights.py` is the single source of the tables: it reads the clock weights `output/55_m6brain_clock/weights/{species}_protein__all__m6brain_clock.tsv` (non-zero genes and the `(intercept)` row), `output/52_proteome_clock/runs/{species}_protein__all.parquet` (`model` = label-free weight, `target` = observed z) and `output/50_transcript_ageing_model/inputs/{species}_protein.parquet` (`m6A_all` is log1p of the site count; the count is stored). It writes the tables with `clock_weight` (0 outside the clock) and `data/version.json` with the build date, the git commit, the date of the clock fit, the number of genes and clock genes, the intercept, the age unit and a hash of each table (`m6brain.weights_version()`).

## Predicted age

For a proteome with genes g and samples j, log abundance x(g,j) (missing allowed), clock weight w(g) for the clock genes (a gene in several rows counts once) and intercept b:

1. Standardize each protein across the samples in which it was quantified:
   z(g,j) = (x(g,j) − mean_j x(g,·)) / sd_j x(g,·), with the sample SD (n − 1). Missing values, and proteins quantified in fewer than 2 samples or with SD 0, get z = 0.
2. Predicted age: a(j) = b + Σ_g w(g) · z(g,j).
3. With chronological ages t(j): Pearson r of a and t over the samples with an age; age acceleration d(j) = a(j) − (c + m t(j)), where c and m are the intercept and slope of the least-squares line of a on t over the control samples (all samples with an age if none are named).

Because z is centered across the samples, the mean of a over the samples is b (exactly without missing values), the training mean age. The units are those of the clock papers: months (mouse) and years (human).

`n_genes_matched(j)` is the number of clock proteins with a defined z in sample j.

This is the calculation of `donor_scores` in `pipeline/55_m6brain_clock.py` (before donor averaging) and of the disease clock in `pipeline/53_disease_clock.py` (acceleration fitted in the controls), and the same as `website/clock.js`; `tests/test_score.py` checks it against a re-implementation of the pipeline function and the Pearson r of the example proteome (Scifo et al. 2025, PXD058684) against `output/55_m6brain_clock/cohorts.tsv` (external set, weights `m6brain_clock`).

## Label-free score

The same standardization (step 1), then the raw score r(j) = Σ_g p(g) · z(g,j) over the genes with a label-free weight p(g) (`prediction`), and the score s(j) = (r(j) − mean_j r) / sd_j r (sample SD). This is `donor_scores` of `pipeline/52_proteome_clock.py`; `tests/test_score.py` checks it against `output/52_proteome_clock/cohorts.tsv` (external set, weights `model`).

## Validation (manuscript)

The Pearson correlation of predicted (or label-free) and chronological age per cohort was pooled by random effects on the Fisher scale.

| | m6Brain clock | label-free |
|---|---|---|
| Mouse, 6 held-out cohorts (weights refitted without the scored cohort) | r 0.97 (95% CI 0.96–0.98) | 0.85 (0.61–0.95) |
| Mouse, 11 external brain proteomes (weights shipped here) | 0.94 (0.88–0.97) | 0.63 (0.52–0.72) |
| Human, 8 held-out cohorts | 0.78 (0.66–0.86) | 0.34 (0.13–0.52) |

Within a cohort (10-fold cross-validation by donor, 10 seeds, the penalty and α chosen in an inner 5-fold cross-validation) the mean absolute error was 5.2 years in BrainGVEX (human, 193 donors; 10.0 years for predicting the training mean) and 1.55 months in Lu 2025 (mouse, 27 donors; 5.95 months for the training mean). In Alzheimer's disease cohorts the clock, with the acceleration fitted in the controls, reads the cases +10.2 years older.

## Limits

The predicted age is relative to the uploaded samples (step 1 standardizes within the upload): the level of a cohort is not an estimate of its age, and the label-free score is not an age. The clock was validated in brain proteomes only; proteins outside the clock do not contribute, and a clock protein that is not quantified contributes 0. Batch, sex, region or cell-composition differences between samples that act on the same proteins move the prediction like age does.
