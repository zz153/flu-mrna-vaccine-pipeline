#!/bin/bash
#SBATCH --job-name=per_year_analysis
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=logs/per_year_%x_%j.out
#SBATCH --error=logs/per_year_%x_%j.err

# Usage: sbatch --export=LINEAGE=H3N2 scripts/slurm/04_per_year_analysis.slurm
# Run from the repository root directory

# Get lineage from environment variable (default to H1N1)
LINEAGE=${LINEAGE:-H1N1}

echo "Starting per-year analysis for ${LINEAGE}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Running on: ${SLURMD_NODENAME}"

# Create logs directory if it doesn't exist
mkdir -p logs

# Activate conda environment
# Note: User must have created the environment using env/environment.yml
eval "$(conda shell.bash hook)"
conda activate flu-vaccine-pipeline

# Run the analysis (assumes script is run from repo root)
bash scripts/04_per_year_analysis.sh ${LINEAGE} 2009 2025 16

echo "Analysis complete for ${LINEAGE}"
