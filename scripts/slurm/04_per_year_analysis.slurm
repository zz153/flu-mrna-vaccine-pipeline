#!/bin/bash
#SBATCH --job-name=per_year_analysis
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=logs/per_year_%x_%j.out
#SBATCH --error=logs/per_year_%x_%j.err

# Usage: 
#   sbatch --export=LINEAGE=H3N2 scripts/slurm/04_per_year_analysis.slurm
#   sbatch --export=LINEAGE=VicB,USE_UNCLUSTERED=true scripts/slurm/04_per_year_analysis.slurm

# Get parameters from environment variables
LINEAGE=${LINEAGE:-H1N1}
USE_UNCLUSTERED=${USE_UNCLUSTERED:-false}

echo "Starting per-year analysis for ${LINEAGE}"
echo "Use unclustered: ${USE_UNCLUSTERED}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Running on: ${SLURMD_NODENAME}"

# Create logs directory if it doesn't exist
mkdir -p logs

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate flu-vaccine-pipeline

# Run the analysis
bash scripts/04_per_year_analysis.sh ${LINEAGE} 2009 2025 16 ${USE_UNCLUSTERED}

echo "Analysis complete for ${LINEAGE}"
