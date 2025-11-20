#!/bin/bash
#SBATCH --job-name=calc_distances
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/distances_%x_%j.out
#SBATCH --error=logs/distances_%x_%j.err

# Usage: sbatch --export=LINEAGE=H1N1 scripts/slurm/05_calculate_distances.slurm
# Run from the repository root directory

# Get lineage from environment variable (default to H1N1)
LINEAGE=${LINEAGE:-H1N1}

echo "Starting distance calculations for ${LINEAGE}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Running on: ${SLURMD_NODENAME}"

# Create logs directory if it doesn't exist
mkdir -p logs

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate flu-vaccine-pipeline

# Run the distance calculation
bash scripts/05_calculate_distances.sh ${LINEAGE} 2009 2025

echo "Distance calculations complete for ${LINEAGE}"
