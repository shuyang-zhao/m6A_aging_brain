# m6Brain

m6Brain, an m6A-aware model of brain proteome aging, provides a zero-shot clock of brain proteome age. Every gene carries a weight, the age effect (z) that a tabular foundation model (TabPFN) predicts from the gene's transcript: Orthrus embedding, RiboNN, m6A sites and reader binding, transcript structure, expression and cell-type origin. The weights were fitted on public brain aging cohorts and are applied to a new proteome as they are: nothing is fitted on your samples. The score of a sample is the weighted sum of its standardized protein abundances.

Manuscript: *Transcript foundation models reveal m6A-linked protein decline and a zero-shot clock of brain proteome aging* (in preparation).

## Install

```bash
pip install "git+https://github.com/shuyang-zhao/m6A_aging_brain#subdirectory=m6brain"
```

Python >= 3.10; dependencies pandas, numpy, scipy.

## Usage

```python
import m6brain as mb

scores = mb.score("proteome.csv", species="mouse")  # sample, score, n_genes_matched
scores = mb.score(df, "human", age=ages)  # adds an age column; scores.attrs["spearman_rho"], ["spearman_p"]
mb.gene("Trmt6", "mouse")  # weight (predicted age z), observed_z, m6A_sites
mb.weights("human")  # the full table: gene_id, symbol, weight, observed_z, m6A_sites
mb.weights_version()  # date and git commit of the weights
```

```bash
m6brain score proteome.csv --species mouse --out scores.csv
```

A worked example (a public mouse whole-brain proteome, 45 samples) is in `examples/quickstart.ipynb`.

## Input

A table of log-scale protein abundances, one row per protein, one column per sample (CSV, or TSV with the extension `.tsv`/`.txt`; a DataFrame works too).

- First column: gene symbol or Ensembl gene id, detected automatically (`ids="symbol"` or `"ensembl"` to force). Version suffixes (`ENSMUSG00000000001.5`) are removed and `A;B` protein groups use the first entry.
- Symbols are case-sensitive in the species convention: mouse `Gnai3`, human `GNAI3`. A warning is raised when few ids match, which usually means the wrong species or the wrong case.
- Missing values: empty or `NA`; they contribute 0 after standardizing. At least 3 samples; the standardization is across the uploaded samples, so the table should hold the samples you want to compare (one brain region, one quantification).
- Optional: a row with the id `age` holds the sample ages (any unit), or pass `age=`; the Spearman correlation of score and age is reported.

## What the score means

The score places each sample on the axis along which proteins predicted to rise with age are high and proteins predicted to fall are low (mean 0, SD 1 across your samples by construction). Higher = older-like proteome within the uploaded set.

It is a zero-shot score: the weights come from transcript features, not from your data. In the manuscript the pooled Spearman correlation of the score with age was rho 0.83 in six held-out mouse brain cohorts (weights fitted without the scored cohort) and 0.66 in eleven external mouse brain cohorts; in eight human brain cohorts it was weaker (rho 0.29).

It is not an age in years, and scores of different uploads are not on a common scale (each is standardized within its upload). It has been validated on brain tissue only, mostly for ages that span months in mouse; with few samples, a narrow age range, or age confounded with batch, region or sex, the correlation says little.

## Layout

`src/m6brain/` package and weight tables (`build_weights.py` regenerates them from the analysis pipeline), `docs/` usage and method, `examples/` notebook and data, `tests/` pytest. Run the tests with `pip install pytest && pytest tests`.

MIT license.
