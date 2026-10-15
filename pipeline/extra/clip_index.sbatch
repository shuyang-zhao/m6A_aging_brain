#!/usr/bin/env bash
#SBATCH --job-name=clip_index
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=16
#SBATCH --mem=90G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.out
# STAR indices of the human and mouse genomes on one node, built by the index step of nf-core/clipseq and kept in
# node-local scratch for the array tasks of clip_peaks.sbatch (the step needs about 80 GB, the tasks 50 GB)
# usage: sbatch --nodelist=NODE clip_index.sbatch SELECTION     (from the directory of clip_peaks.sh)
set -euo pipefail

export THREADS=$SLURM_CPUS_PER_TASK MEMORY_GB=84
export CACHE=/scratch/clip_peaks_$USER/cache
for species in human mouse; do
  row=$(awk -F'\t' -v s=$species \
    'NR == 1 {for (i = 1; i <= NF; i++) col[$i] = i; next} $col["species"] == s {print NR - 1; exit}' "$1")
  WORK=/scratch/clip_peaks_$USER/index_$SLURM_JOB_ID INDEX_ONLY=1 \
    bash "$SLURM_SUBMIT_DIR/clip_peaks.sh" "$1" "$row" /tmp/clip_index_unused
done
