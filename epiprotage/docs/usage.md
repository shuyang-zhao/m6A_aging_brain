# Usage

## Python

```python
import epiprotage as ep

result = ep.score("proteome.csv", species="mouse")
```

`score(table, species, ids="auto", age=None)`

| argument | |
|---|---|
| `table` | DataFrame or path to a CSV/TSV. Gene ids in the first column (a DataFrame with a non-default index uses the index), one column per sample, log abundances. Paths ending in `.tsv`, `.txt` or `.tab` (also `.gz`) are read as tab-separated, others as comma-separated. |
| `species` | `"mouse"` or `"human"`. |
| `ids` | `"auto"` (default), `"symbol"` or `"ensembl"`. Auto takes Ensembl ids when more than half of the rows look like `ENSG...` / `ENSMUSG...`. |
| `age` | Ages of the samples: a Series indexed by sample name or a list in column order. Overrides a row named `age` (case-insensitive) in the table. |

Returns a DataFrame with `sample`, `score` (mean 0, SD 1), `n_genes_matched` (weighted proteins quantified in that sample) and, with ages, `age`. `result.attrs` holds `species`, `ids`, `genes_input`, `genes_matched` and, with ages, `spearman_rho` and `spearman_p`.

Warnings: fewer than half of the ids match (wrong species or case; the message says how many match the other species and how many match ignoring case); Errors: fewer than 3 samples, no id matches. A gene in several rows counts once (the first row).

## Gene weights

```python
ep.gene(
    "Trmt6", "mouse"
)  # case-insensitive fallback for single lookups; Ensembl ids (with or without version) work too
ep.weights("human")  # DataFrame: gene_id, symbol, weight, observed_z, m6A_sites
ep.weights_version()  # {"date": ..., "git_commit": ..., "mouse": {...}, "human": {...}}
```

- `weight`: predicted age z of the gene (positive = predicted to increase with age). This is what `score` uses.
- `observed_z`: combined age z of the gene over the training cohorts (the target of the model).
- `m6A_sites`: number of m6A sites of the gene (mouse: the two GLORI hippocampus maps; human: brain m6A-seq peaks).

Genes without a row in the table (about 10,000 per species, the genes of the model) are ignored by `score`.

## Command line

```bash
epiprotage score FILE --species mouse|human [--ids auto|symbol|ensembl] [--out scores.csv]
epiprotage --version
```

The scores go to `--out` or stdout (CSV); the number of matched rows and, if the file has an `age` row, the Spearman correlation go to stderr.

## Checking a result

- `result.attrs["genes_matched"]` should be a large part of the proteins you quantified (thousands in a typical brain proteome).
- Plot score against age (see the notebook); the score is only meaningful relative to the other samples of the same upload.
