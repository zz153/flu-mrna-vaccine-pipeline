#!/bin/bash
# Script 05: Calculate Evolutionary Distances for All Vaccine Designs
#
# Calculates distances between vaccine designs and circulating strains
# For each design type:
#   - Consensus
#   - Medoid
#   - Ancestral
# 
# Uses aligned versions of designs to ensure proper comparison
# Distance methods:
# 1. Simple p-distance (observed differences, ignoring gap-gap pairs)
# 2. LG+G4 evolutionary distance (from IQ-TREE ML estimation, only for medoid)

set -e

LINEAGE=$1
YEAR_START=${2:-2009}
YEAR_END=${3:-2025}
SOURCE_TYPE=${4:-clustered}

if [ -z "$LINEAGE" ]; then
    echo "Usage: bash scripts/05_calculate_distances.sh [LINEAGE] [START_YEAR] [END_YEAR] [SOURCE_TYPE]"
    echo "Example: bash scripts/05_calculate_distances.sh H1N1 2009 2025 clustered"
    echo ""
    echo "SOURCE_TYPE: clustered or unclustered"
    exit 1
fi

echo "=========================================="
echo "Calculating Distances for ${LINEAGE}"
echo "Years: ${YEAR_START} to ${YEAR_END}"
echo "Source: ${SOURCE_TYPE}"
echo "=========================================="

INPUT_DIR="results/per_year_${SOURCE_TYPE}/${LINEAGE}"
OUTPUT_DIR="${INPUT_DIR}/distances"
mkdir -p ${OUTPUT_DIR}

# Design types to analyze
DESIGNS=("consensus" "medoid" "ancestral")

# Summary file for each design type
for DESIGN in "${DESIGNS[@]}"; do
    SUMMARY="${OUTPUT_DIR}/distance_summary_${DESIGN}.csv"
    echo "Year,N_Strains,Mean_P_Distance,Median_P_Distance,SD_P_Distance,Mean_ML_Distance,Median_ML_Distance,SD_ML_Distance" > ${SUMMARY}
done

for YEAR in $(seq ${YEAR_START} ${YEAR_END}); do
    
    ALIGNMENT="${INPUT_DIR}/alignments/${LINEAGE}_${YEAR}_aligned.fasta"
    MLDIST="${INPUT_DIR}/trees/${LINEAGE}_${YEAR}.mldist"
    
    # Skip if files don't exist
    if [ ! -f "$ALIGNMENT" ] || [ ! -f "$MLDIST" ]; then
        echo "⏭️  ${YEAR}: Missing required files"
        continue
    fi
    
    echo ""
    echo "📅 Processing ${YEAR}..."
    
    # Process each design type
    for DESIGN in "${DESIGNS[@]}"; do
        
        DESIGN_FILE="${INPUT_DIR}/designs/${LINEAGE}_${YEAR}_${DESIGN}_aligned.fasta"
        OUTPUT_CSV="${OUTPUT_DIR}/${LINEAGE}_${YEAR}_${DESIGN}_distances.csv"
        
        if [ ! -f "$DESIGN_FILE" ]; then
            echo "  ⏭️  ${DESIGN}: Design file not found"
            continue
        fi
        
        echo "  🔄 Calculating distances for ${DESIGN}..."
        
        # Calculate distances using Python
        python3 << EOF
import numpy as np
from Bio import AlignIO, SeqIO
import sys

design_type = "${DESIGN}"

# Read alignment
try:
    alignment = AlignIO.read("${ALIGNMENT}", "fasta")
except Exception as e:
    print(f"  ❌ Error reading alignment: {e}")
    sys.exit(1)

# Read design sequence (aligned version)
try:
    design_record = SeqIO.read("${DESIGN_FILE}", "fasta")
    design_seq_aligned = str(design_record.seq)
except Exception as e:
    print(f"  ❌ Error reading design: {e}")
    sys.exit(1)

# Find design in alignment (for medoid - it's an actual strain)
design_idx = None
alignment_list = list(alignment)

for i, record in enumerate(alignment_list):
    # For medoid, match by checking if IDs are similar
    if design_type == "medoid" and "|" in design_record.id:
        # Extract strain ID from medoid header
        design_strain_id = design_record.id.split("|")[1] if "|" in design_record.id else design_record.id
        if design_strain_id in record.id:
            design_idx = i
            design_seq_aligned = str(record.seq)  # Use the one from alignment
            break

# Calculate p-distances using aligned sequences
p_distances = []
strain_ids = []

for i, record in enumerate(alignment_list):
    # Skip if this is the design itself (medoid case)
    if design_idx is not None and i == design_idx:
        continue
    
    strain_seq_aligned = str(record.seq)
    
    # Compare aligned sequences, excluding gap-gap positions
    valid_positions = []
    for j in range(len(design_seq_aligned)):
        design_char = design_seq_aligned[j]
        strain_char = strain_seq_aligned[j]
        
        # Skip positions where both have gaps
        if design_char == '-' and strain_char == '-':
            continue
        # Skip positions where either has a gap (for p-distance)
        if design_char == '-' or strain_char == '-':
            continue
        
        valid_positions.append((design_char, strain_char))
    
    if len(valid_positions) == 0:
        continue
    
    # Calculate p-distance
    differences = sum(1 for d, s in valid_positions if d != s)
    p_dist = differences / len(valid_positions)
    
    p_distances.append(p_dist)
    strain_ids.append(record.id)

print(f"  Calculated p-distances for {len(p_distances)} strains")

# Read ML distances (only available for medoid since it's in the tree)
ml_dist_to_design = []

if design_idx is not None:
    try:
        with open("${MLDIST}", 'r') as f:
            lines = f.readlines()
        
        n_seqs = int(lines[0].strip())
        
        for i in range(n_seqs):
            if i == design_idx:
                continue
            
            if i > design_idx:
                dist_line = lines[i + 1].strip().split()
                if len(dist_line) > design_idx + 1:
                    ml_dist_to_design.append(float(dist_line[design_idx + 1]))
                else:
                    ml_dist_to_design.append(np.nan)
            elif i < design_idx:
                dist_line = lines[design_idx + 1].strip().split()
                if len(dist_line) > i + 1:
                    ml_dist_to_design.append(float(dist_line[i + 1]))
                else:
                    ml_dist_to_design.append(np.nan)
        
        print(f"  Extracted {len(ml_dist_to_design)} ML distances")
    
    except Exception as e:
        print(f"  ⚠️  Could not parse ML distances: {e}")
        ml_dist_to_design = [np.nan] * len(p_distances)
else:
    print(f"  ⚠️  Design not in tree, ML distances unavailable")
    ml_dist_to_design = [np.nan] * len(p_distances)

# Ensure arrays are same length
min_len = min(len(p_distances), len(ml_dist_to_design), len(strain_ids))
p_distances = p_distances[:min_len]
ml_dist_to_design = ml_dist_to_design[:min_len]
strain_ids = strain_ids[:min_len]

# Write distances to CSV
with open("${OUTPUT_CSV}", 'w') as f:
    f.write("Strain_ID,P_Distance,ML_Distance_LG_G4,Difference\n")
    for sid, pd, mld in zip(strain_ids, p_distances, ml_dist_to_design):
        diff = mld - pd if not np.isnan(mld) else np.nan
        f.write(f"{sid},{pd:.6f},{mld:.6f},{diff:.6f}\n")

# Calculate summary statistics
if len(p_distances) > 0:
    p_mean = np.mean(p_distances)
    p_median = np.median(p_distances)
    p_std = np.std(p_distances)
    
    ml_valid = [x for x in ml_dist_to_design if not np.isnan(x)]
    if ml_valid:
        ml_mean = np.mean(ml_valid)
        ml_median = np.median(ml_valid)
        ml_std = np.std(ml_valid)
    else:
        ml_mean = ml_median = ml_std = np.nan
    
    print(f"  P-distance: {p_mean:.6f} ± {p_std:.6f} (median: {p_median:.6f})")
    if not np.isnan(ml_mean):
        print(f"  ML distance: {ml_mean:.6f} ± {ml_std:.6f} (median: {ml_median:.6f})")
    
    # Append to summary
    summary_file = "${OUTPUT_DIR}/distance_summary_${DESIGN}.csv"
    with open(summary_file, 'a') as f:
        f.write(f"${YEAR},{len(p_distances)},{p_mean:.6f},{p_median:.6f},{p_std:.6f},")
        f.write(f"{ml_mean:.6f},{ml_median:.6f},{ml_std:.6f}\n")
    
    print(f"  ✅ Distances saved")
else:
    print(f"  ❌ No distances calculated")

EOF

        if [ $? -eq 0 ]; then
            echo "  ✅ ${DESIGN} complete"
        else
            echo "  ❌ ${DESIGN} failed"
        fi
        
    done
    
    echo "✅ ${YEAR} complete"
    
done

echo ""
echo "=========================================="
echo "✅ Distance analysis complete for ${LINEAGE}!"
echo "=========================================="
echo ""

# Display summaries
for DESIGN in "${DESIGNS[@]}"; do
    SUMMARY="${OUTPUT_DIR}/distance_summary_${DESIGN}.csv"
    if [ -f "$SUMMARY" ]; then
        echo "=== ${DESIGN} Summary ==="
        cat ${SUMMARY}
        echo ""
    fi
done
