#!/bin/bash -l
#SBATCH -p short_gpu
#SBATCH --gres=gpu:ada6000:1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -t 01:30:00
#SBATCH --array=1-16%2
#SBATCH -J proteome_clock
#SBATCH -o proteome_clock_%A_%a.log
# The GPU step of script 52, one leave-one-cohort-out table per array task, or the branch given as argument
# (`sbatch --array=1 proteome_clock.sbatch unseen`) (UCR HPCC, RTX 6000 Ada). Run from a folder holding
# pipeline/utils.py, pipeline/52_proteome_clock.py and the output of `52_proteome_clock.py inputs`, with WEIGHTS the
# folder of the TabPFN-3.5 regressor weights; the environment is built in node-local scratch
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
TABLES=(output/52_proteome_clock/inputs/*.parquet)
name=${1:-$(basename "${TABLES[$SLURM_ARRAY_TASK_ID - 1]}" .parquet)}
out=output/52_proteome_clock/runs/$name.parquet
[ "$name" = unseen ] && out=output/52_proteome_clock/unseen/predictions.parquet
[ -e "$out" ] && { echo "done already: $name"; exit 0; }

W=${TMPDIR:?}
export TABPFN_MODEL_CACHE_DIR=$W/models TABPFN_NO_BROWSER=1 PIP_CACHE_DIR=$W/pip XDG_CACHE_HOME=$W/xdg PYTHONUNBUFFERED=1
mkdir "$W/models" && cp "$WEIGHTS"/tabpfn-v3.5-*.safetensors "$W/models/"
/opt/linux/rocky/8.x/x86_64/pkgs/miniconda3/py312_24.7.1/bin/python -m venv "$W/env"
"$W/env/bin/pip" install -q tabpfn==9.1.0 pyarrow
"$W/env/bin/pip" freeze | grep -iE "^(tabpfn|torch|numpy|pandas|scipy|scikit-learn)==" || true
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader || true
"$W/env/bin/python" pipeline/52_proteome_clock.py "$name"; rc=$?
echo "=== $name rc=$rc $(date)"
exit $rc
