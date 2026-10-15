#!/usr/bin/env bash
#SBATCH --job-name=clip_peaks
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --output=%x_%A_%a.out
# One row of the selection table per array task, on the UCR HPCC cluster; run clip_index.sbatch on the nodes first
# usage: sbatch --array=1-23%23 --exclude=NODES clip_peaks.sbatch SELECTION RESULTS
#        (from the directory of clip_peaks.sh; NODES: the epyc nodes without STAR indices)
# PureCLIP, which learns over all contigs, needs more than the 50 GB of this header for the large human libraries:
# submit those rows with --mem=150G (the pipeline's memory limit follows the job's memory)
# Node-local scratch holds everything but the peaks: STAR indices and tools in /scratch/clip_peaks_$USER/cache,
# FASTQ files and the Nextflow work directory in one directory per task, deleted when the task ends. RESULTS (shared
# disk) receives the peaks and the provenance.
set -euo pipefail

export THREADS=$SLURM_CPUS_PER_TASK MEMORY_GB=$((SLURM_MEM_PER_NODE / 1024 - 6))  # the job's memory less 6 GB
NODE_DIR=/scratch/clip_peaks_$USER  # /scratch/$USER belongs to root on the compute nodes
export CACHE=$NODE_DIR/cache
export WORK=$NODE_DIR/task_$SLURM_JOB_ID
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

# the array starts at once: spread the first requests to ENA, ENCODE and GitHub
sleep $((RANDOM % 90))
bash "$SLURM_SUBMIT_DIR/clip_peaks.sh" "$1" "$SLURM_ARRAY_TASK_ID" "$2"
