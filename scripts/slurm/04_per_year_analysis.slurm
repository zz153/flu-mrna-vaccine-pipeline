#!/bin/bash
#SBATCH --job-name=per_year_analysis
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/per_year_%x_%j.out
#SBATCH --error=logs/per_year_%x_%j.err

# Usage: 
#   For clustered (H1N1, H3N2):
#     sbatch --export=LINEAGE=H1N1 scripts/slurm/04_per_year_analysis.slurm
#   For unclustered (VicB):
#     sbatch --export=LINEAGE=VicB,USE_UNCLUSTERED=true scripts/slurm/04_per_year_analysis.slurm

# Get parameters from environment variables
LINEAGE=${LINEAGE:-H1N1}
YEAR_START=${YEAR_START:-2009}
YEAR_END=${YEAR_END:-2025}
USE_UNCLUSTERED=${USE_UNCLUSTERED:-false}
THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "Starting per-year analysis for ${LINEAGE}"
echo "Years: ${YEAR_START} to ${YEAR_END}"
echo "Use unclustered: ${USE_UNCLUSTERED}"
echo "Threads: ${THREADS}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Running on: ${SLURMD_NODENAME}"

# Create logs directory if it doesn't exist
mkdir -p logs

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate flu-vaccine-pipeline

# Run the per-year analysis
bash scripts/04_per_year_analysis.sh ${LINEAGE} ${YEAR_START} ${YEAR_END} ${THREADS} ${USE_UNCLUSTERED}

echo "Per-year analysis complete for ${LINEAGE}"
