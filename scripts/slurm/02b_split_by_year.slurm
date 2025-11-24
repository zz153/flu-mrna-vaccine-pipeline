#!/bin/bash
#SBATCH --job-name=split_years
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/split_years_%x_%j.out
#SBATCH --error=logs/split_years_%x_%j.err

# Usage (default years 2009–2025):
#   sbatch --export=LINEAGE=H1N1 scripts/slurm/02b_split_by_year.slurm
#
# Or override years:
#   sbatch --export=LINEAGE=H1N1,START_YEAR=2010,END_YEAR=2024 scripts/slurm/02b_split_by_year.slurm

LINEAGE=${LINEAGE:-H1N1}
START_YEAR=${START_YEAR:-2009}
END_YEAR=${END_YEAR:-2025}

echo "Starting per-year splitting for ${LINEAGE} (${START_YEAR}–${END_YEAR})"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Running on: ${SLURMD_NODENAME}"

# Ensure logs directory exists
mkdir -p logs

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate flu-vaccine-pipeline

# Run Script 02b (split clustered sequences by year)
bash scripts/02b_split_by_year.sh "${LINEAGE}" "${START_YEAR}" "${END_YEAR}"

echo "Per-year splitting complete for ${LINEAGE} (${START_YEAR}–${END_YEAR})"
