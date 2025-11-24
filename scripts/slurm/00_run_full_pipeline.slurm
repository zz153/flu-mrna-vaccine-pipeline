#!/bin/bash
#SBATCH --job-name=flu_pipeline
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=logs/pipeline_master_%j.out
#SBATCH --error=logs/pipeline_master_%j.err

# Master SLURM Script: Complete Influenza Vaccine Design Pipeline
#
# This script orchestrates the entire analysis pipeline by submitting
# all jobs with proper dependencies for all three lineages.
#
# Usage:
#   sbatch scripts/slurm/00_run_full_pipeline.slurm
#
# This will submit jobs for:
#   - H1N1 (clustered)
#   - H3N2 (clustered)
#   - VicB (unclustered)

set -e

# Navigate to project root (one level up from scripts directory)
cd ..
echo "Working directory: $(pwd)"
echo ""

echo "=========================================="
echo "FluHub: Master Pipeline Orchestrator"
echo "=========================================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo ""

# Create logs directory if it doesn't exist
mkdir -p logs

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate flu-vaccine-pipeline

# Function to submit pipeline for a lineage with dependencies
submit_lineage_pipeline() {
    local lineage=$1
    local use_unclustered=$2
    
    echo "----------------------------------------"
    echo "Submitting pipeline for ${lineage}"
    echo "Unclustered mode: ${use_unclustered}"
    echo "----------------------------------------"
    
    local export_base="LINEAGE=${lineage}"
    
    # Step 1: Filter sequences
    local job1=$(sbatch --parsable --export=${export_base} \
        scripts/slurm/01_filter_sequences.slurm)
    echo "  [1/6] Filter sequences: Job ${job1}"
    
    if [ "$use_unclustered" = "false" ]; then
        # CLUSTERED PIPELINE (H1N1, H3N2)
        
        # Step 2: Cluster sequences
        local job2=$(sbatch --parsable --dependency=afterok:${job1} \
            --export=${export_base} \
            scripts/slurm/02_cluster_sequences.slurm)
        echo "  [2/6] Cluster sequences: Job ${job2}"
        
        # Step 3: Split by year
        local job3=$(sbatch --parsable --dependency=afterok:${job2} \
            --export=${export_base} \
            scripts/slurm/02b_split_by_year.slurm)
        echo "  [3/6] Split by year: Job ${job3}"
        
        # Step 4: Per-year analysis & vaccine design (does its own alignment)
        local job4=$(sbatch --parsable --dependency=afterok:${job3} \
            --export=${export_base} \
            scripts/slurm/04_per_year_analysis.slurm)
        echo "  [4/6] Vaccine design: Job ${job4}"
        
        # Step 5: Calculate distances
        local job5=$(sbatch --parsable --dependency=afterok:${job4} \
            --export=${export_base} \
            scripts/slurm/05_calculate_distances.slurm)
        echo "  [5/6] Calculate distances: Job ${job5}"
        
        # Step 6: Visualizations
        local job6=$(sbatch --parsable --dependency=afterok:${job5} \
            --export=${export_base} \
            scripts/slurm/06_visualize_distances.slurm)
        echo "  [6/6] Visualizations: Job ${job6}"
        
    else
        # UNCLUSTERED PIPELINE (VicB)
        
        echo "  [2/6] Skipping clustering (unclustered mode)"
        
        # Step 3: Split by year
        local job3=$(sbatch --parsable --dependency=afterok:${job1} \
            --export=${export_base},USE_UNCLUSTERED=true \
            scripts/slurm/02b_split_by_year.slurm)
        echo "  [3/6] Split by year: Job ${job3}"
        
        # Step 4: Per-year analysis & vaccine design
        local job4=$(sbatch --parsable --dependency=afterok:${job3} \
            --export=${export_base},USE_UNCLUSTERED=true \
            scripts/slurm/04_per_year_analysis.slurm)
        echo "  [4/6] Vaccine design: Job ${job4}"
        
        # Step 5: Calculate distances
        local job5=$(sbatch --parsable --dependency=afterok:${job4} \
            --export=${export_base},SOURCE_TYPE=unclustered \
            scripts/slurm/05_calculate_distances.slurm)
        echo "  [5/6] Calculate distances: Job ${job5}"
        
        # Step 6: Visualizations
        local job6=$(sbatch --parsable --dependency=afterok:${job5} \
            --export=${export_base},SOURCE_TYPE=unclustered \
            scripts/slurm/06_visualize_distances.slurm)
        echo "  [6/6] Visualizations: Job ${job6}"
    fi
    
    echo ""
}

# Check for data availability and submit pipelines
echo "Checking data availability..."
echo ""

# H1N1 (clustered)
if [ -f "data/raw/H1N1_raw.fasta" ]; then
    echo "✓ H1N1 data found"
    submit_lineage_pipeline H1N1 false
else
    echo "⚠ H1N1 data not found - skipping"
    echo ""
fi

# H3N2 (clustered)
if [ -f "data/raw/H3N2_raw.fasta" ]; then
    echo "✓ H3N2 data found"
    submit_lineage_pipeline H3N2 false
else
    echo "⚠ H3N2 data not found - skipping"
    echo ""
fi

# VicB (unclustered)
if [ -f "data/raw/VicB_raw.fasta" ]; then
    echo "✓ VicB data found"
    submit_lineage_pipeline VicB true
else
    echo "⚠ VicB data not found - skipping"
    echo ""
fi

echo "=========================================="
echo "All pipeline jobs submitted!"
echo "=========================================="
echo ""
echo "Monitor progress with:"
echo "  squeue -u \$USER"
echo ""
echo "Check logs in:"
echo "  logs/"
echo ""
echo "Completed: $(date)"
