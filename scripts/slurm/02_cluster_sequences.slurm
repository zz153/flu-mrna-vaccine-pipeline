#!/bin/bash
#SBATCH --job-name=cluster_seqs
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/cluster_%x_%j.out
#SBATCH --error=logs/cluster_%x_%j.err

# Usage:
#   sbatch --export=LINEAGE=H1N1 scripts/slurm/02_cluster_sequences.slurm
#   sbatch --export=LINEAGE=H3N2 scripts/slurm/02_cluster_sequences.slurm
#   sbatch --export=LINEAGE=VicB  scripts/slurm/02_cluster_sequences.slurm

LINEAGE=${LINEAGE:-H1N1}

echo "Starting CD-HIT clustering for ${LINEAGE}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Running on: ${SLURMD_NODENAME}"

# Ensure logs directory exists
mkdir -p logs

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate flu-vaccine-pipeline

# Run Script 02 (CD-HIT 99% clustering)
bash scripts/02_cluster_sequences.sh "${LINEAGE}"

echo "Clustering complete for ${LINEAGE}"
