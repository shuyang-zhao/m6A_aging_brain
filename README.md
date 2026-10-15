# m6A load marks proteins that decline in the ageing mammalian brain

Analysis code for the article of the same title (Genome Biology, in preparation). The study combines 30 public brain
ageing datasets (90 analysis units of mouse and human proteomes and transcriptomes) with transcriptome-wide m6A maps,
tests whether the number of m6A sites of a transcript marks the decline of its protein with age, and relates the m6A
readers to the output of their marked transcripts. All inputs are public; `data_sources.tsv` lists every dataset with
its citation, DOI, accession and location.

## Contents

- `pipeline/`: numbered scripts, run in order from the repository root; each begins with what it does, what it reads
  and what it writes
  - `00_download.sh`: the public inputs, into `data/raw/`
  - `01`-`08`: proteomes, transcriptomes, single-nucleus pseudobulks, m6A and other modification site maps, gene
    annotation, m6A reader CLIP peaks and the analysis units (`output/prepared/`)
  - `09`: per-gene transcript features
  - `10`, `12`, `13`: age effect of every gene in every unit; age permutations; concordance between species
  - `20`, `21`, `25`: signatures and load of the modified transcripts, with the successive adjustment, the DRACH,
    stoichiometry and cell-composition checks and the age-permutation null of the m6A load; causal forests
  - `26`, `28`: validation of the m6A results in human and their dependence on the site map
  - `30`, `33`, `34`, `37`, `38`, `40`: age effects of the RNA-modification regulators, the m6A regulator-target
    coupling network, its sign check against knockouts, the compendium of reader binding and perturbation data, the
    reader tests and the reader CLIP targets
  - `50`: models of gene-level ageing from transcript features (TabPFN, with LightGBM and an elastic net)
  - `60`: source data of every figure panel and Tables S1-S7 (`source_data/`)
  - `61`: the results workbooks, Additional files 3 and 4 (`output/61_supplementary_data/supplementary_data_1.xlsx`
    and `supplementary_data_2.xlsx`)
  - `62`: the query website of the per-gene age effects (`output/62_website/`)
  - `extra/`: steps run once on other hardware: CLIP peaks with nf-core/clipseq on a SLURM cluster, transcript
    embeddings and the TabPFN step of script 50 on a GPU; `select_clip_datasets.py` chose the CLIP libraries from the
    library catalogue of a separate project (clipflow), which is not part of this repository
  - `utils.py`, `model.py`, `network.py`, `single_nucleus.py` and the `.R` files: functions shared by the scripts
- `figures/`: one R script per figure (`fig1.R`-`fig5.R`; Additional file 1: `sfig1.R`-`sfig6.R`) and their common
  style (`style.R`); the figures read only `source_data/`
- `website/`: the pages of the query website (static HTML, CSS and JavaScript; script 62 adds the data)
- `source_data/`: the plotted values of every panel and Tables S1-S7, as written by script 60
- `ref/gene_sets/rna_mod_regulators.tsv`: the curated RNA-modification regulators
- `docs/supplementary_data_columns.tsv`: meaning and unit of every column of the results workbooks

## System requirements

- Linux on x86-64 (tested on Ubuntu 22.04, kernel 5.15). Python 3.12 with the packages and versions of
  `environment.yml` (pandas 3.0, numpy 2.5, scipy 1.18, statsmodels 0.15, scikit-learn 1.9, LightGBM 4.7).
- R 4.6.1 on the PATH (the scripts call `Rscript`) with Bioconductor limma 3.68.5, edgeR 4.10.5 and
  SummarizedExperiment 1.42.0 and CRAN lme4 2.0-6, mashr 0.2.79, grf 2.6.1, data.table 1.18.6.1, ggplot2 4.0.3,
  patchwork 1.3.2, ggrepel 0.9.8, ggbeeswarm 0.7.3, ggforce 0.5.0, png 0.1-9 and systemfonts 1.3.2 (the versions
  the results were produced with). mashr needs the GSL library (script 33 also looks for it in `env/lib`). The
  figures use the Arial font, which must be installed.
- The results were produced on one machine with 56 cores and 251 GB of memory. The pipeline needs about 45 GB for the
  raw inputs (the files of the raw-data manifest in Additional file 3) and 30 GB for the outputs.
- No special hardware for the figures and the CPU steps. The TabPFN step of script 50 and the transcript embeddings
  ran on NVIDIA RTX 6000 Ada GPUs (48 GB); the CLIP peaks ran with nf-core/clipseq 1.0.0 (Nextflow 22.10.6, Apptainer
  1.5.4 with the Singularity profile) on a SLURM cluster. Their products (the CLIP binding regions in
  `data/raw/clip_standard`, the Orthrus and RiboNN predictions in `data/raw/gpu_precomputed` and the TabPFN runs in
  `output/50_transcript_ageing_model/runs`) belong to the data release, so the rest of the pipeline runs without
  them.

## Installation

```bash
mamba env create -p env -f environment.yml          # Python
Rscript -e 'install.packages(c("lme4", "mashr", "grf", "data.table", "ggplot2", "patchwork", "ggrepel", "ggbeeswarm",
  "ggforce", "png", "systemfonts", "BiocManager"))'
Rscript -e 'BiocManager::install(c("limma", "edgeR", "SummarizedExperiment"))'
```

## Demo: the figures from the source data

`source_data/` holds the plotted values of every panel, so the figures can be drawn without running the pipeline:

```bash
for figure in figures/fig*.R figures/sfig*.R; do Rscript "$figure"; done
```

This writes `output/figures/fig1-5` and `sfig1-6` (PDF and PNG), the figures of the article and of Additional file 1,
in about 2 minutes on one core.

## Reproducing the analysis

```bash
bash pipeline/00_download.sh                         # public inputs; files kept with the data release are marked
for script in pipeline/0[1-9]_*.py; do env/bin/python "$script"; done
for script in 10 12 13 20 21 25 28 30 26 33 34 37 40 38; do env/bin/python pipeline/${script}_*.py; done
env/bin/python pipeline/33_regulator_network.py wildtype
env/bin/python pipeline/50_transcript_ageing_model.py inputs
WEIGHTS=/path/to/tabpfn_3_5 sbatch pipeline/extra/transcript_model.sbatch  # TabPFN, 18 GPU tasks of 3-37 minutes
env/bin/python pipeline/50_transcript_ageing_model.py
env/bin/python pipeline/60_source_data.py
env/bin/python pipeline/61_supplementary_data.py
env/bin/python pipeline/62_website.py
for figure in figures/fig*.R figures/sfig*.R; do Rscript "$figure"; done
```

The inputs that no script can download come with the data release: the files marked in `00_download.sh` and the
reader binding lists curated from published supplements in `ref/target_sets/` (read by scripts 33, 34, 37, 38 and 62).
Scripts 09-50 and 61 write `output/<script name>/`, scripts 01-08 write `output/prepared/` and script 60 writes
`source_data/`. The sbatch files name the partitions and Python of the cluster they ran on (UCR HPCC); adapt them to
another cluster. Run times on the machine above: the preparation steps 01-08 about 35 minutes, the age permutations (12)
about 3 hours, the coupling network (33) about 25 minutes and its wild-type refit 5, the reader tests (38) about 50
minutes, the CPU steps of script 50 about 25 minutes, the results workbooks (61) about 6 minutes and the website (62)
about 8 minutes; every other script takes a few minutes.

## Website

https://shuyang-zhao.github.io/m6A_aging_brain/

Script 62 writes the website to `output/62_website/`: the pages of `website/` and one JSON file per gene (33,229
genes, about 80 MB). Searching a gene shows its age effect in every unit, grouped by species and layer, with its m6A
sites and m6A reader binding; the landscape view combines every gene over the cohorts of each context with the rule
of the regulator screen (script 30) and ranks the genes (`data/landscape/landscape.tsv.gz` holds every column); the
cohort view is the interactive form of Fig. 1a and Table S1, and every dataset, m6A map and binding source carries its
reference. Preview it with `env/bin/python -m http.server -d output/62_website`. GitHub Pages serves the built site
from the `gh-pages` branch of this repository.

## Licence

MIT (`LICENSE`).
