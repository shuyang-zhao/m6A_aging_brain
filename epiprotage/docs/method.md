# Method

## Weights

For each species, TabPFN-3.5 (a tabular foundation model) was fitted to the combined age z of the genes in public brain proteome aging cohorts (Stouffer combination, weighted by √n, of the per-cohort age effects of each protein) using the transcript features of the gene: Orthrus embedding, RiboNN predictions, m6A sites per region (mouse: the two GLORI hippocampus maps; human: brain m6A-seq peaks) and reader binding, transcript structure (lengths, exons, GC, junction density), expression, and the cell-type origin of the transcript (`pipeline/50_transcript_ageing_model.py`). Each gene is predicted by a model fitted on genes of other chromosome groups (five groups), and this out-of-fold prediction is the weight of the gene (`pipeline/52_proteome_clock.py`, run `all`). The weights are the column `weight` of `src/epiprotage/data/{mouse,human}_weights.tsv`; `observed_z` is the target and `m6A_sites` the site count of the gene.

`build_weights.py` is the single source of the tables: it reads `output/52_proteome_clock/runs/{species}_protein__all.parquet` (`model` = weight, `target` = observed z) and `output/50_transcript_ageing_model/inputs/{species}_protein.parquet` (`m6A_all` is log1p of the site count; the count is stored) and writes `data/version.json` with the build date, the git commit, the date of the run and a hash of each table (`epiprotage.weights_version()`).

## Score

For a proteome with genes g and samples j, log abundance x(g,j) (missing allowed), and weight w(g) for the genes in the table (one row per gene; a gene in several rows counts once):

1. Standardize each protein across the samples in which it was quantified:
   z(g,j) = (x(g,j) − mean_j x(g,·)) / sd_j x(g,·), with the sample SD (n − 1). Missing values, and proteins quantified in fewer than 2 samples or with SD 0, get z = 0.
2. Raw score of a sample: r(j) = Σ_g w(g) · z(g,j), over the genes with a weight.
3. Score: s(j) = (r(j) − mean_j r) / sd_j r (sample SD).

`n_genes_matched(j)` is the number of weighted proteins with a defined z in sample j.

This is the calculation of `donor_scores` in `pipeline/52_proteome_clock.py` (before donor averaging) and of `website/clock.js`; `tests/test_score.py` checks it against a re-implementation of the pipeline function and the Spearman correlation of the example proteome (Scifo et al. 2025, PXD058684) against `output/52_proteome_clock/cohorts.tsv` (external set, model weights).

## Validation (manuscript)

The Spearman correlation of score and age per cohort was pooled by random effects on the Fisher scale. Mouse, 6 held-out cohorts (each scored with weights refitted without it): rho 0.83 (95% CI 0.59–0.93). Mouse, 11 external cohorts scored with the weights shipped here: rho 0.66 (0.53–0.76). Human, 8 held-out cohorts: rho 0.29 (0.11–0.46), at the level of a length-only model. The shipped human weights are therefore of limited use.

## Limits

The score is relative to the uploaded samples (step 1 and 3 standardize within the upload), is not an age in years, and was validated in brain proteomes only. Weights are a prediction for ~10,000 genes per species; proteins outside the table do not contribute. Batch, sex, region or cell-composition differences between samples that act on the same proteins move the score like age does.
