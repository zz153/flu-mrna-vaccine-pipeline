#!/bin/bash
#SBATCH --job-name=calc_distances
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/distances_%x_%j.out
#SBATCH --error=logs/distances_%x_%j.err

# Usage: 
#   For clustered (H1N1, H3N2):
#     sbatch --export=LINEAGE=H1N1 scripts/slurm/05_calculate_distances.slurm
#   For unclustered (VicB):
#     sbatch --export=LINEAGE=VicB,SOURCE_TYPE=unclustered scripts/slurm/05_calculate_distances.slurm

# Get parameters from environment variables
LINEAGE=${LINEAGE:-H1N1}
SOURCE_TYPE=${SOURCE_TYPE:-clustered}

echo "Starting distance calculations for ${LINEAGE}"
echo "Source type: ${SOURCE_TYPE}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Running on: ${SLURMD_NODENAME}"

# Create logs directory if it doesn't exist
mkdir -p logs

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate flu-vaccine-pipeline

# Run the distance calculation
bash scripts/05_calculate_distances.sh ${LINEAGE} 2009 2025 ${SOURCE_TYPE}

echo "Distance calculations complete for ${LINEAGE}"
