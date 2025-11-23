#!/bin/bash
# Master Script: Complete Influenza Vaccine Design Pipeline
#
# This script orchestrates the entire analysis pipeline from raw sequences
# to final visualizations for all three lineages (H1N1, H3N2, VicB)
#
# Usage:
#   bash scripts/run_full_pipeline.sh [EXECUTION_MODE]
#   
# Execution modes:
#   local  - Run all scripts locally (sequential)
#   slurm  - Submit all jobs to SLURM (parallel, recommended for HPC)
#   check  - Only check prerequisites and data availability

set -e

EXECUTION_MODE=${1:-local}

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo ""
echo "=========================================="
echo "FluHub: Influenza Vaccine Design Pipeline"
echo "=========================================="
echo ""

# Function to print colored messages
print_status() {
    echo -e "${BLUE}[$(date +%H:%M:%S)]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[✓]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[!]${NC} $1"
}

print_error() {
    echo -e "${RED}[✗]${NC} $1"
}

# Function to check prerequisites
check_prerequisites() {
    print_status "Checking prerequisites..."
    
    local all_good=true
    
    # Check conda environment
    if conda env list | grep -q "flu-vaccine-pipeline"; then
        print_success "Conda environment 'flu-vaccine-pipeline' found"
    else
        print_error "Conda environment 'flu-vaccine-pipeline' not found"
        echo "         Create it with: conda env create -f env/environment.yml"
        all_good=false
    fi
    
    # Check for raw data files
    for lineage in H1N1 H3N2 VicB; do
        if [ -f "data/raw/${lineage}_raw.fasta" ]; then
            print_success "Raw data found: ${lineage}_raw.fasta"
        else
            print_warning "Raw data missing: ${lineage}_raw.fasta"
            echo "         Place GISAID sequences in data/raw/"
        fi
    done
    
    # Check required tools
    if command -v mafft &> /dev/null; then
        print_success "MAFFT installed"
    else
        print_error "MAFFT not found"
        all_good=false
    fi
    
    if command -v iqtree &> /dev/null; then
        print_success "IQ-TREE installed"
    else
        print_error "IQ-TREE not found"
        all_good=false
    fi
    
    if command -v cd-hit &> /dev/null; then
        print_success "CD-HIT installed"
    else
        print_error "CD-HIT not found"
        all_good=false
    fi
    
    if ! $all_good; then
        print_error "Prerequisites check failed!"
        exit 1
    fi
    
    print_success "All prerequisites satisfied!"
    echo ""
}

# Function to run pipeline locally
run_local_pipeline() {
    local lineage=$1
    local use_unclustered=$2
    
    print_status "Running pipeline for ${lineage}..."
    
    # Script 01: Filter sequences
    print_status "Step 1/6: Filtering sequences..."
    bash scripts/01_filter_sequences.sh ${lineage}
    print_success "Filtering complete"
    
    # Script 02: Cluster sequences (skip for unclustered)
    if [ "$use_unclustered" = "false" ]; then
        print_status "Step 2/6: Clustering sequences at 99% identity..."
        bash scripts/02_cluster_sequences.sh ${lineage}
        print_success "Clustering complete"
    else
        print_status "Step 2/6: Skipping clustering (unclustered mode)"
    fi
    
    # Script 02b: Split by year
    print_status "Step 3/6: Splitting sequences by year..."
    bash scripts/02b_split_by_year.sh ${lineage} ${use_unclustered}
    print_success "Year splitting complete"
    
    # Script 03: Combined alignment (skip for unclustered)
    if [ "$use_unclustered" = "false" ]; then
        print_status "Step 4/6: Creating combined alignment..."
        bash scripts/03_combined_alignment.sh ${lineage}
        print_success "Combined alignment complete"
    else
        print_status "Step 4/6: Skipping combined alignment (unclustered mode)"
    fi
    
    # Script 04: Per-year analysis and vaccine design
    print_status "Step 5/6: Generating vaccine designs for all years..."
    bash scripts/04_per_year_analysis.sh ${lineage} 2009 2025 4 ${use_unclustered}
    print_success "Vaccine designs complete"
    
    # Script 05: Calculate distances
    print_status "Step 6/6: Calculating evolutionary distances..."
    if [ "$use_unclustered" = "true" ]; then
        bash scripts/05_calculate_distances.sh ${lineage} 2009 2025 unclustered
    else
        bash scripts/05_calculate_distances.sh ${lineage} 2009 2025 clustered
    fi
    print_success "Distance calculations complete"
    
    # Script 06: Visualizations
    print_status "Step 7/6: Creating visualizations..."
    if [ "$use_unclustered" = "true" ]; then
        bash scripts/06_visualize_distances.sh ${lineage} unclustered
    else
        bash scripts/06_visualize_distances.sh ${lineage} clustered
    fi
    print_success "Visualizations complete"
    
    print_success "Pipeline complete for ${lineage}!"
    echo ""
}

# Function to submit SLURM jobs
run_slurm_pipeline() {
    local lineage=$1
    local use_unclustered=$2
    local export_args="LINEAGE=${lineage}"
    
    if [ "$use_unclustered" = "true" ]; then
        export_args="${export_args},USE_UNCLUSTERED=true"
    fi
    
    print_status "Submitting SLURM jobs for ${lineage}..."
    
    # Submit jobs with dependencies
    local job1=$(sbatch --parsable --export=${export_args} scripts/slurm/01_filter_sequences.slurm)
    print_status "Job ${job1}: Filter sequences"
    
    if [ "$use_unclustered" = "false" ]; then
        local job2=$(sbatch --parsable --dependency=afterok:${job1} --export=${export_args} scripts/slurm/02_cluster_sequences.slurm)
        print_status "Job ${job2}: Cluster sequences"
        
        local job2b=$(sbatch --parsable --dependency=afterok:${job2} --export=${export_args} scripts/slurm/02b_split_by_year.slurm)
        print_status "Job ${job2b}: Split by year"
        
        local job3=$(sbatch --parsable --dependency=afterok:${job2b} --export=${export_args} scripts/slurm/03_combined_alignment.slurm)
        print_status "Job ${job3}: Combined alignment"
        
        local job4=$(sbatch --parsable --dependency=afterok:${job3} --export=${export_args} scripts/slurm/04_per_year_analysis.slurm)
        print_status "Job ${job4}: Vaccine design"
        
        local job5=$(sbatch --parsable --dependency=afterok:${job4} --export=${export_args} scripts/slurm/05_calculate_distances.slurm)
        print_status "Job ${job5}: Calculate distances"
        
        local job6=$(sbatch --parsable --dependency=afterok:${job5} --export=${export_args} scripts/slurm/06_visualize_distances.slurm)
        print_status "Job ${job6}: Visualize results"
    else
        local job2b=$(sbatch --parsable --dependency=afterok:${job1} --export=${export_args} scripts/slurm/02b_split_by_year.slurm)
        print_status "Job ${job2b}: Split by year"
        
        local job4=$(sbatch --parsable --dependency=afterok:${job2b} --export=${export_args} scripts/slurm/04_per_year_analysis.slurm)
        print_status "Job ${job4}: Vaccine design"
        
        local job5=$(sbatch --parsable --dependency=afterok:${job4} --export=${export_args},SOURCE_TYPE=unclustered scripts/slurm/05_calculate_distances.slurm)
        print_status "Job ${job5}: Calculate distances"
        
        local job6=$(sbatch --parsable --dependency=afterok:${job5} --export=${export_args},SOURCE_TYPE=unclustered scripts/slurm/06_visualize_distances.slurm)
        print_status "Job ${job6}: Visualize results"
    fi
    
    print_success "All jobs submitted for ${lineage}!"
    echo ""
}

# Main execution
case $EXECUTION_MODE in
    check)
        check_prerequisites
        exit 0
        ;;
    local)
        check_prerequisites
        
        print_status "Running pipeline in LOCAL mode (sequential execution)"
        echo ""
        
        # H1N1 (clustered)
        if [ -f "data/raw/H1N1_raw.fasta" ]; then
            run_local_pipeline H1N1 false
        else
            print_warning "Skipping H1N1 (no raw data)"
        fi
        
        # H3N2 (clustered)
        if [ -f "data/raw/H3N2_raw.fasta" ]; then
            run_local_pipeline H3N2 false
        else
            print_warning "Skipping H3N2 (no raw data)"
        fi
        
        # VicB (unclustered)
        if [ -f "data/raw/VicB_raw.fasta" ]; then
            run_local_pipeline VicB true
        else
            print_warning "Skipping VicB (no raw data)"
        fi
        
        print_success "Full pipeline complete!"
        ;;
    slurm)
        check_prerequisites
        
        print_status "Running pipeline in SLURM mode (parallel execution)"
        echo ""
        
        # H1N1 (clustered)
        if [ -f "data/raw/H1N1_raw.fasta" ]; then
            run_slurm_pipeline H1N1 false
        else
            print_warning "Skipping H1N1 (no raw data)"
        fi
        
        # H3N2 (clustered)
        if [ -f "data/raw/H3N2_raw.fasta" ]; then
            run_slurm_pipeline H3N2 false
        else
            print_warning "Skipping H3N2 (no raw data)"
        fi
        
        # VicB (unclustered)
        if [ -f "data/raw/VicB_raw.fasta" ]; then
            run_slurm_pipeline VicB true
        else
            print_warning "Skipping VicB (no raw data)"
        fi
        
        print_success "All jobs submitted!"
        print_status "Monitor with: squeue -u \$USER"
        ;;
    *)
        print_error "Invalid execution mode: ${EXECUTION_MODE}"
        echo ""
        echo "Usage: bash scripts/run_full_pipeline.sh [EXECUTION_MODE]"
        echo ""
        echo "Execution modes:"
        echo "  local  - Run all scripts locally (sequential)"
        echo "  slurm  - Submit all jobs to SLURM (parallel, recommended for HPC)"
        echo "  check  - Only check prerequisites and data availability"
        exit 1
        ;;
esac

echo ""
echo "=========================================="
echo "Pipeline Execution Complete!"
echo "=========================================="
