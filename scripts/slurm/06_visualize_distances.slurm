#!/bin/bash
#SBATCH --job-name=visualize
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/visualize_%x_%j.out
#SBATCH --error=logs/visualize_%x_%j.err

# Usage: 
#   For clustered (H1N1, H3N2):
#     sbatch --export=LINEAGE=H1N1 scripts/slurm/06_visualize_distances.slurm
#   For unclustered (VicB):
#     sbatch --export=LINEAGE=VicB,SOURCE_TYPE=unclustered scripts/slurm/06_visualize_distances.slurm

# Get parameters from environment variables
LINEAGE=${LINEAGE:-H1N1}
SOURCE_TYPE=${SOURCE_TYPE:-clustered}

echo "Starting visualization for ${LINEAGE}"
echo "Source type: ${SOURCE_TYPE}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Running on: ${SLURMD_NODENAME}"

# Create logs directory if it doesn't exist
mkdir -p logs

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate flu-vaccine-pipeline

# Run the visualization
bash scripts/06_visualize_distances.sh ${LINEAGE} ${SOURCE_TYPE}

echo "Visualization complete for ${LINEAGE}"
