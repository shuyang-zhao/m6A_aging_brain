# Transcript foundation models reveal m6A-linked protein decline and guide a proteomic clock of the aging brain

Analysis code of the article of the same title (Genome Biology, in preparation). The study combines 39 public brain
datasets (101 analysis units of mouse and human proteomes and transcriptomes: mouse discovery cohorts, human validation
cohorts and eleven external mouse proteomes) with transcriptome-wide m6A maps. In m6Brain, an m6A-aware model of brain
proteome aging, a tabular foundation model (TabPFN) predicts the age effect of every gene from the features of its
transcript; the m6A sites of a transcript mark the decline of its protein with age, and the predictions guide an
elastic net of age on brain proteomes, the m6Brain clock, which reads neurodegeneration as years of aging. The m6A
readers and their marked targets are tested alongside. All inputs are public; `data_sources.tsv` lists every dataset with its citation, DOI,
accession and location, and `ref/datasets.tsv` (Table S1) every dataset considered, with its set or the reason it was
not used.

## Contents

- `pipeline/`: numbered scripts, run in order from the repository root; each begins with what it does, what it reads
  and what it writes
  - `00_download.sh`: the public inputs, into `data/raw/`
  - `01`-`08`: proteomes, transcriptomes, single-nucleus pseudobulks, m6A and other modification site maps, gene
    annotation, m6A reader CLIP peaks and the analysis units (`output/prepared/`)
  - `09`: per-gene transcript features
  - `10`, `12`, `13`: age effect of every gene in every unit; age permutations; concordance between species
  - `20`, `21`: signatures and load of the modified transcripts, with the successive adjustment, the DRACH,
    stoichiometry and cell-composition checks and the age-permutation null of the m6A load
  - `28`: dependence of the m6A results on the site map
  - `30`, `33`, `34`, `37`, `38`, `40`: age effects of the RNA-modification regulators, the m6A regulator-target
    coupling network, its sign check against knockouts, the compendium of reader binding and perturbation data, the
    reader tests and the reader CLIP targets
  - `50`: models of gene-level ageing from transcript features (TabPFN, with LightGBM and an elastic net)
  - `52`: the label-free proteome score of m6Brain: leave-one-cohort-out weights, scores of the held-out and external
    cohorts and the control weights
  - `53`: the m6Brain clock in Alzheimer's, Parkinson's, Huntington's and amyotrophic lateral sclerosis cohorts
  - `54`: the individual m6A index and its titration against the reader levels
  - `55`: the m6Brain clock (an elastic net of age guided by the m6Brain weights of script 52): leave-one-cohort-out
    and all-cohort weights, the benchmark against 27 published clocks and the within-cohort errors; run after 52
  - `60`: source data of every figure panel and Tables S1-S7 (`source_data/`)
  - `61`: Additional files 2-5 as Excel workbooks (`output/61_supplementary_data/`): the Tables S1-S7, the per-gene age
    effects, the clock weights and gene-level model predictions, and the clock scores of every sample with the
    individual m6A index
  - `62`: the query website of the per-gene age effects and the browser m6Brain clock (`output/62_website/`)
  - `extra/`: steps run once on other hardware: CLIP peaks with nf-core/clipseq on a SLURM cluster, transcript
    embeddings and the TabPFN steps of scripts 50 and 52 on a GPU; `select_clip_datasets.py` chose the CLIP libraries
    from the library catalogue of a separate project (clipflow), which is not part of this repository
  - `utils.py`, `model.py`, `network.py`, `single_nucleus.py` and the `.R` files: functions shared by the scripts
- `figures/`: one R script per figure (`fig1.R`-`fig6.R`; Additional file 1: `sfig1.R`-`sfig10.R`) and their common
  style (`style.R`); the figures read only `source_data/`
- `m6brain/`: the m6Brain clock Python package (weights, command line, tests and a worked example)
- `website/`: the pages of the query website (static HTML, CSS and JavaScript; script 62 adds the data)
- `source_data/`: the plotted values of every panel and Tables S1-S7, as written by script 60
- `ref/datasets.tsv`: every dataset considered (Table S1); `ref/gene_sets/rna_mod_regulators.tsv`: the curated
  RNA-modification regulators
- `docs/supplementary_data_columns.tsv`: meaning and unit of every column of the Additional file workbooks

Datasets that were prepared but fail the inclusion criteria as aging cohorts (a preprint and a knock-in allelic series)
are listed in `utils.EXCLUDED`: script 08 writes no units for them, so no later script sees them, and removing an entry
and rerunning from script 08 includes the dataset again.

A new proteome cohort, such as an in-house aging cohort, enters as an external cohort: a block in
`01_prepare_proteomes.py`, its units in `08_define_units.py` and their names in `utils.EXTERNAL`; then scripts 10, 12,
21, 52 and 55 score it with the clock of the training cohorts and test it for the m6A load, and nothing else uses it.

## m6Brain

The m6Brain clock predicts the age of brain proteome samples (months in mouse, years in human) with an elastic net of
age on the standardized abundances of 279 mouse and 421 human proteins, guided by the m6Brain weights of script 52.
Pooled over the held-out cohorts the Pearson correlation of predicted and chronological age is 0.97 in mouse and 0.78 in
human, and 0.94 in eleven external mouse brain proteomes. A label-free mode (`score`) uses the m6Brain weights alone.

```bash
pip install "git+https://github.com/shuyang-zhao/m6A_aging_brain#subdirectory=m6brain"
```

```python
import m6brain as mb

ages = mb.predict_age("proteome.csv", species="mouse")  # sample, predicted_age, n_genes_matched
```

```bash
m6brain age proteome.csv --species mouse --out ages.csv
```

The weights, input format, the age acceleration fitted in control samples and a worked example are in
`m6brain/README.md`; `m6brain/build_weights.py` rebuilds the weight tables from the output of scripts 52 and 55.

## System requirements

- Linux on x86-64 (tested on Ubuntu 22.04, kernel 5.15). Python 3.12 with the packages and versions of
  `environment.yml` (pandas 3.0, numpy 2.5, scipy 1.18, statsmodels 0.15, scikit-learn 1.9, LightGBM 4.7).
- R 4.6.1 on the PATH (the scripts call `Rscript`) with Bioconductor limma 3.68.5, edgeR 4.10.5 and
  SummarizedExperiment 1.42.0 and CRAN lme4 2.0-6, mashr 0.2.79, data.table 1.18.6.1, ggplot2 4.0.3,
  patchwork 1.3.2, ggrepel 0.9.8, ggbeeswarm 0.7.3, ggforce 0.5.0, png 0.1-9 and systemfonts 1.3.2 (the versions
  the results were produced with). mashr needs the GSL library (script 33 also looks for it in `env/lib`). The
  figures use the Arial font, which must be installed.
- The results were produced on one machine with 56 cores and 251 GB of memory. The pipeline needs about 45 GB for the
  raw inputs (the files listed in `data_sources.tsv`) and 30 GB for the outputs.
- No special hardware for the figures and the CPU steps. The TabPFN steps of scripts 50 and 52 and the transcript
  embeddings ran on NVIDIA RTX 6000 Ada GPUs (48 GB); the CLIP peaks ran with nf-core/clipseq 1.0.0 (Nextflow 22.10.6,
  Apptainer 1.5.4 with the Singularity profile) on a SLURM cluster. Their products (the CLIP binding regions in
  `data/raw/clip_standard`, the Orthrus and RiboNN predictions in `data/raw/gpu_precomputed` and the TabPFN runs in
  `output/50_transcript_ageing_model/runs` and `output/52_proteome_clock/runs`) belong to the data release, so the rest
  of the pipeline runs without them.

## Installation

```bash
mamba env create -p env -f environment.yml          # Python
Rscript -e 'install.packages(c("lme4", "mashr", "data.table", "ggplot2", "patchwork", "ggrepel", "ggbeeswarm",
  "ggforce", "png", "systemfonts", "BiocManager"))'
Rscript -e 'BiocManager::install(c("limma", "edgeR", "SummarizedExperiment"))'
```

## Demo: the figures from the source data

`source_data/` holds the plotted values of every panel, so the figures can be drawn without running the pipeline:

```bash
for figure in figures/fig*.R figures/sfig*.R; do Rscript "$figure"; done
```

This writes `output/figures/fig1-6` and `sfig1-10` (PDF and PNG), the figures of the article and of Additional file 1,
in about 2 minutes on one core.

## Reproducing the analysis

```bash
bash pipeline/00_download.sh                         # public inputs; files kept with the data release are marked
for script in pipeline/0[1-9]_*.py; do env/bin/python "$script"; done
for script in 10 12 13 20 21 28 30 33 34 37 40 38; do env/bin/python pipeline/${script}_*.py; done
env/bin/python pipeline/50_transcript_ageing_model.py inputs
WEIGHTS=/path/to/tabpfn_3_5 sbatch pipeline/extra/transcript_model.sbatch  # TabPFN, 15 GPU tasks of 3-37 minutes
env/bin/python pipeline/50_transcript_ageing_model.py
env/bin/python pipeline/52_proteome_clock.py inputs
WEIGHTS=/path/to/tabpfn_3_5 sbatch pipeline/extra/proteome_clock.sbatch    # TabPFN, 16 GPU tasks of a few minutes
env/bin/python pipeline/52_proteome_clock.py
OMP_NUM_THREADS=1 PYTHONWARNINGS=ignore env/bin/python pipeline/55_m6brain_clock.py    # hours on the first run
env/bin/python pipeline/53_disease_clock.py
env/bin/python pipeline/54_reader_titration.py
env/bin/python pipeline/60_source_data.py
env/bin/python pipeline/61_supplementary_data.py
env/bin/python pipeline/62_website.py
for figure in figures/fig*.R figures/sfig*.R; do Rscript "$figure"; done
```

The TabPFN-3.5 regressor weights (Hugging Face `Prior-Labs/tabpfn_3_5`, non-commercial licence) are not part of this
repository. The inputs that no script can download come with the data release: the files marked in `00_download.sh`
and the reader binding lists curated from published supplements in `ref/target_sets/` (read by scripts 33, 34, 37, 38
and 62). Scripts 09-54 and 61 write `output/<script name>/`, scripts 01-08 write `output/prepared/` and script 60
writes `source_data/`. The sbatch files name the partitions and Python of the cluster they ran on (UCR HPCC); adapt
them to another cluster. Run times on the machine above: the preparation steps 01-08 about 35 minutes, the age
permutations (12) about 3 hours, the coupling network (33) about 7 minutes, the reader tests (38) about 40 minutes,
the CPU steps of script 50 about 25 minutes, the Additional file workbooks (61) about 3 minutes and the website (62)
about 6 minutes; every other script takes a few minutes.

## Website

https://shuyang-zhao.github.io/m6A_aging_brain/

Script 62 writes the website to `output/62_website/`: the pages of `website/` and one JSON file per gene (29,524
genes, about 60 MB). Searching a gene shows its age effect in every unit, grouped by species and layer (protein,
mRNA-adjusted protein and mRNA), with its m6A sites and m6A reader binding; the landscape view combines every gene
over the cohorts of each context with the rule of the regulator screen (script 30) and ranks the genes
(`data/landscape/landscape.tsv.gz` holds every column); the cohort view is the interactive form of Fig. 1b and
Table S1, and every dataset, m6A map and binding source carries its reference. The Clock page predicts the age of the
samples of a proteome table with the m6Brain clock in the browser (nothing is uploaded; `website/clock.js`). Preview it
with `env/bin/python -m http.server -d output/62_website`. GitHub Pages serves the built site from the `gh-pages`
branch of this repository.

## Licence

MIT (`LICENSE`; `m6brain/LICENSE` for the package).
