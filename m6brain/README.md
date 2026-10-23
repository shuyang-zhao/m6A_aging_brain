# m6Brain

The m6Brain clock predicts the age of brain proteome samples from their protein abundances, in months for mouse and in years for human. It is an elastic net of age on the standardized abundances of 279 mouse and 421 human proteins, fitted on public brain aging cohorts. The penalty of each protein is guided by m6Brain, an m6A-aware model of brain proteome aging that predicts the age effect of every gene from the features of its transcript (TabPFN: Orthrus embedding, RiboNN, m6A sites and reader binding, transcript structure, expression and cell-type origin). In leave-one-cohort-out tests the pooled Pearson correlation of predicted and chronological age was r = 0.97 in 6 held-out mouse cohorts and 0.78 in 8 held-out human cohorts; the clock fitted on all training cohorts, as shipped here, gave r = 0.94 in 11 external mouse brain proteomes. Within a cohort the mean absolute error was 1.55 months (Lu 2025, mouse) and 5.2 years (BrainGVEX, human). The clock was validated on brain tissue only.

A second, label-free mode (`score`) weights the proteins by the age effect that m6Brain predicts from the transcript alone, with nothing fitted on ages (pooled r 0.85, 0.63 and 0.34 in the same three sets).

Manuscript: *Transcript foundation models reveal m6A-linked protein decline and guide a proteomic clock of the aging brain* (in preparation).

## Install

```bash
pip install "git+https://github.com/shuyang-zhao/m6A_aging_brain#subdirectory=m6brain"
```

Python >= 3.10; dependencies pandas, numpy, scipy.

## Usage

```python
import m6brain as mb

ages = mb.predict_age("proteome.csv", species="mouse")  # sample, predicted_age (months), n_genes_matched
ages = mb.predict_age(df, "human", age=years, controls=names)  # adds age, age_acceleration; ages.attrs["pearson_r"], ["pearson_p"]
mb.score("proteome.csv", "mouse")  # label-free mode: sample, score, n_genes_matched
mb.gene("Trmt6", "mouse")  # clock_weight, prediction (label-free weight), observed_z, m6A_sites
mb.weights("human")  # the full table: gene_id, symbol, clock_weight, prediction, observed_z, m6A_sites
mb.weights_version()  # date, git commit, intercept and size of the weights
```

```bash
m6brain age proteome.csv --species mouse --out ages.csv  # [--controls S1,S2,...]
m6brain score proteome.csv --species mouse --out scores.csv
```

A worked example (a public mouse whole-brain proteome, 45 samples) is in `examples/quickstart.ipynb`.

## Input

A table of log-scale protein abundances, one row per protein, one column per sample (CSV, or TSV with the extension `.tsv`/`.txt`; a DataFrame works too).

- First column: gene symbol or Ensembl gene id, detected automatically (`ids="symbol"` or `"ensembl"` to force). Version suffixes (`ENSMUSG00000000001.5`) are removed and `A;B` protein groups use the first entry.
- Symbols are case-sensitive in the species convention: mouse `Gnai3`, human `GNAI3`. A warning is raised when few ids match, which usually means the wrong species or the wrong case.
- Missing values: empty or `NA`; they contribute 0 after standardizing. At least 3 samples; the standardization is across the uploaded samples, so the table should hold the samples you want to compare (one brain region, one quantification).
- Optional: a row with the id `age` holds the sample ages (months in mouse and years in human for the predicted age to be comparable), or pass `age=`; the Pearson correlation of prediction and age is reported.

## What the prediction means

The predicted age of a sample is the intercept plus the sum over the clock proteins of weight x standardized abundance, in the unit of the training ages. The proteins are standardized across your samples (mean 0, SD 1), so the mean predicted age of a table is the training mean age (mouse 13.2 months, human 69.4 years; exactly when no value is missing). Differences between samples carry the information; the level of a table does not, and a mean predicted age far from the true mean age is not an error of the clock. Three outputs are meaningful: the predicted age of each sample relative to the others, its Pearson r with chronological age, and the age acceleration, the residual of the predicted on the chronological age in a linear fit. For a disease comparison, fit the acceleration in the controls (`controls=` the names of the control samples) and apply it to all samples; this is how the manuscript measures disease acceleration (Alzheimer's disease +10.2 years).

The clock was validated on brain tissue only, mostly for ages that span months in mouse; with few samples, a narrow age range, or age confounded with batch, region or sex, the correlation says little. Clock proteins that are not in the table contribute 0 (`attrs["genes_matched"]` of `attrs["clock_genes"]` are found); the label-free score is standardized within the upload and is not an age.

## Layout

`src/m6brain/` package and weight tables (`build_weights.py` regenerates them from the analysis pipeline), `docs/` usage and method, `examples/` notebook and data, `tests/` pytest. Run the tests with `pip install pytest && pytest tests`.

MIT license.
