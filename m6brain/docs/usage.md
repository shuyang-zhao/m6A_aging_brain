# Usage

## Python

```python
import m6brain as mb

result = mb.predict_age("proteome.csv", species="mouse")
```

`predict_age(table, species, age=None, controls=None, ids="auto")`

| argument | |
|---|---|
| `table` | DataFrame or path to a CSV/TSV. Gene ids in the first column (a DataFrame with a non-default index uses the index), one column per sample, log abundances. Paths ending in `.tsv`, `.txt` or `.tab` (also `.gz`) are read as tab-separated, others as comma-separated. |
| `species` | `"mouse"` (predicted age in months) or `"human"` (years). |
| `age` | Ages of the samples: a Series indexed by sample name or a list in column order. Overrides a row named `age` (case-insensitive) in the table. |
| `controls` | Names of the samples on which the age acceleration is fitted; default all samples. |
| `ids` | `"auto"` (default), `"symbol"` or `"ensembl"`. Auto takes Ensembl ids when more than half of the rows look like `ENSG...` / `ENSMUSG...`. |

Returns a DataFrame with `sample`, `predicted_age`, `n_genes_matched` (clock proteins quantified in that sample) and, with ages, `age` and `age_acceleration` (residual of `predicted_age` on `age` in a linear fit over the controls, or over all samples with an age). `result.attrs` holds `species`, `ids`, `genes_input` (rows of the table), `genes_matched` (clock proteins found), `clock_genes` (279 mouse, 421 human), `unit` and, with ages, `pearson_r` and `pearson_p` (all samples with an age).

Warnings: fewer than half of the ids match the species table (wrong species or case; the message says how many match the other species and how many match ignoring case); Errors: fewer than 3 samples, no id matches a clock protein. A gene in several rows counts once (the first row).

```python
mb.score(table, species, ids="auto", age=None)  # label-free mode
```

Returns `sample`, `score` (mean 0, SD 1 across the samples), `n_genes_matched` (proteins with a label-free weight quantified in the sample, about 10,000 genes per species) and, with ages, `age`; `result.attrs["pearson_r"]` and `["pearson_p"]` hold the correlation of score and age.

## Gene weights

```python
mb.gene(
    "Trmt6", "mouse"
)  # case-insensitive fallback for single lookups; Ensembl ids (with or without version) work too
mb.weights("human")  # DataFrame: gene_id, symbol, clock_weight, prediction, observed_z, m6A_sites
mb.weights_version()  # {"date": ..., "git_commit": ..., "mouse": {"intercept": ..., "unit": ...}, "human": {...}}
```

- `clock_weight`: coefficient of the m6Brain clock in months (mouse) or years (human) per SD of the protein; 0 for the genes outside the clock. The intercept of each species is in `weights_version()[species]["intercept"]`.
- `prediction`: label-free weight, the age z of the gene that the TabPFN model predicts from its transcript (positive = predicted to increase with age). This is what `score` uses.
- `observed_z`: combined age z of the gene over the training cohorts (the target of the model).
- `m6A_sites`: number of m6A sites of the gene (mouse: the two GLORI hippocampus maps; human: brain m6A-seq peaks).

The table holds the genes of the model (about 10,000 per species); the clock proteins are among them, and genes without a row are ignored.

## Command line

```bash
m6brain age FILE --species mouse|human [--ids auto|symbol|ensembl] [--controls S1,S2,...] [--out ages.csv]
m6brain score FILE --species mouse|human [--ids auto|symbol|ensembl] [--out scores.csv]
m6brain --version
```

The result goes to `--out` or stdout (CSV); the number of matched proteins and, if the file has an `age` row, the Pearson correlation go to stderr. `--controls` takes the sample names, separated by commas.

## Checking a result

- `result.attrs["genes_matched"]` should be a large part of the 279 (mouse) or 421 (human) clock proteins; the example proteome of 4,714 proteins holds 190 of the 279 mouse clock proteins.
- Plot `predicted_age` against `age` (see the notebook). The mean predicted age equals the training mean age by construction, so judge the slope and the correlation, not the offset.
