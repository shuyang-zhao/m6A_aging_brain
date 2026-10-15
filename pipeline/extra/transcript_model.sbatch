#!/bin/bash
#SBATCH -p short_gpu
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 02:00:00
#SBATCH --array=1-18
#SBATCH -J transcript_model
#SBATCH -o transcript_model_%A_%a.log
# The GPU steps of script 50, one array task per GPU (UCR HPCC, RTX 6000 Ada; at most about 40 minutes per task): the
# five targets, their permuted copies in two halves and the human targets with the cell-line sites. Run from a folder
# holding pipeline/utils.py, pipeline/50_transcript_ageing_model.py, ref/gene_sets and the output of
# `50_transcript_ageing_model.py inputs`, with WEIGHTS the folder of the TabPFN-3.5 regressor weights (Hugging Face
# Prior-Labs/tabpfn_3_5, non-commercial licence); the environment is built in node-local scratch
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
TARGETS=(mouse_protein human_mrna human_neurons mouse_mrna_adjusted human_protein)
TASKS=("${TARGETS[@]}")
for target in "${TARGETS[@]}"; do TASKS+=("$target 1-50" "$target 51-100"); done
TASKS+=(human_mrna_cell_line_sites human_neurons_cell_line_sites human_protein_cell_line_sites)

W=$(mktemp -d "${TMPDIR:-/tmp}/transcript_model_XXXXXX")  # per task: tasks on one node share /tmp
trap 'rm -rf "$W"' EXIT
export TABPFN_MODEL_CACHE_DIR=$W/models TABPFN_NO_BROWSER=1 PIP_CACHE_DIR=$W/pip
mkdir "$W/models" && cp "$WEIGHTS"/tabpfn-v3.5-*.safetensors "$W/models/"
/opt/linux/rocky/8.x/x86_64/pkgs/miniconda3/py312_24.7.1/bin/python -m venv "$W/env"
"$W/env/bin/pip" install -q tabpfn==9.1.0 pyarrow
"$W/env/bin/pip" freeze | grep -iE "^(tabpfn|torch|numpy|pandas|scipy|scikit-learn)=="
nvidia-smi --query-gpu=name --format=csv,noheader
"$W/env/bin/python" pipeline/50_transcript_ageing_model.py ${TASKS[$SLURM_ARRAY_TASK_ID - 1]}
