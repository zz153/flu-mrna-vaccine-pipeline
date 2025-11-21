#!/bin/bash
#SBATCH --job-name=filter_seqs
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/filter_%x_%j.out
#SBATCH --error=logs/filter_%x_%j.err

# Usage:
#   sbatch --export=LINEAGE=H1N1 scripts/slurm/01_filter_sequences.slurm
#   sbatch --export=LINEAGE=H3N2 scripts/slurm/01_filter_sequences.slurm
#   sbatch --export=LINEAGE=VicB  scripts/slurm/01_filter_sequences.slurm

LINEAGE=${LINEAGE:-H1N1}

echo "Starting sequence filtering for ${LINEAGE}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Running on: ${SLURMD_NODENAME}"

# Ensure logs directory exists
mkdir -p logs

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate flu-vaccine-pipeline

# Run Script 01 (filter + clean)
bash scripts/01_filter_sequences.sh "${LINEAGE}"

echo "Filtering complete for ${LINEAGE}"
