#!/bin/bash
# Script 05: Calculate Evolutionary Distances
#
# Calculates distances between consensus designs and circulating strains
# Uses two methods:
# 1. Simple p-distance (observed differences)
# 2. LG+G4 evolutionary distance (from IQ-TREE ML estimation)

set -e

LINEAGE=$1
YEAR_START=${2:-2009}
YEAR_END=${3:-2025}

if [ -z "$LINEAGE" ]; then
    echo "Usage: bash scripts/05_calculate_distances.sh [LINEAGE] [START_YEAR] [END_YEAR]"
    echo "Example: bash scripts/05_calculate_distances.sh H1N1 2009 2025"
    exit 1
fi

echo "=========================================="
echo "Calculating Distances for ${LINEAGE}"
echo "Years: ${YEAR_START} to ${YEAR_END}"
echo "=========================================="

INPUT_DIR="results/per_year_clustered/${LINEAGE}"
OUTPUT_DIR="${INPUT_DIR}/distances"
mkdir -p ${OUTPUT_DIR}

# Summary file
SUMMARY="${OUTPUT_DIR}/distance_summary.csv"
echo "Year,N_Strains,Mean_P_Distance,Median_P_Distance,SD_P_Distance,Mean_ML_Distance,Median_ML_Distance,SD_ML_Distance" > ${SUMMARY}

for YEAR in $(seq ${YEAR_START} ${YEAR_END}); do
    
    CONSENSUS="${INPUT_DIR}/designs/${LINEAGE}_${YEAR}_consensus.fasta"
    ALIGNMENT="${INPUT_DIR}/alignments/${LINEAGE}_${YEAR}_aligned.fasta"
    MLDIST="${INPUT_DIR}/trees/${LINEAGE}_${YEAR}.mldist"
    
    # Skip if files don't exist
    if [ ! -f "$CONSENSUS" ] || [ ! -f "$ALIGNMENT" ] || [ ! -f "$MLDIST" ]; then
        echo "⏭️  ${YEAR}: Missing required files"
        continue
    fi
    
    echo ""
    echo "📅 Processing ${YEAR}..."
    
    OUTPUT_CSV="${OUTPUT_DIR}/${LINEAGE}_${YEAR}_distances.csv"
    
    # Calculate distances using Python
    python3 << EOF
import numpy as np
from Bio import AlignIO, SeqIO
import sys

# Read alignment
try:
    alignment = AlignIO.read("${ALIGNMENT}", "fasta")
except Exception as e:
    print(f"Error reading alignment: {e}")
    sys.exit(1)

# Find consensus sequence
consensus_seq = None
consensus_idx = None
for i, record in enumerate(alignment):
    if "Consensus" in record.id or "consensus" in record.id.lower():
        consensus_seq = str(record.seq)
        consensus_idx = i
        break

if consensus_seq is None:
    print("Error: Consensus sequence not found in alignment")
    sys.exit(1)

print(f"  Found consensus at position {consensus_idx}")

# Calculate p-distances to all other sequences
p_distances = []
strain_ids = []

for i, record in enumerate(alignment):
    if i == consensus_idx:
        continue
    
    seq = str(record.seq)
    
    # Remove gaps for p-distance calculation
    pairs = [(consensus_seq[j], seq[j]) for j in range(len(consensus_seq)) 
             if consensus_seq[j] != '-' and seq[j] != '-']
    
    if len(pairs) == 0:
        continue
    
    differences = sum(1 for a, b in pairs if a != b)
    p_dist = differences / len(pairs)
    
    p_distances.append(p_dist)
    strain_ids.append(record.id)

print(f"  Calculated p-distances for {len(p_distances)} strains")

# Read ML distances from .mldist file
try:
    with open("${MLDIST}", 'r') as f:
        lines = f.readlines()
    
    n_seqs = int(lines[0].strip())
    
    # Parse distance matrix
    ml_distances = []
    seq_names = []
    
    for i in range(1, n_seqs + 1):
        parts = lines[i].strip().split()
        seq_names.append(parts[0])
    
    # Extract distances to consensus
    ml_dist_to_consensus = []
    
    if consensus_idx < len(seq_names):
        for i in range(n_seqs):
            if i == consensus_idx:
                continue
            
            # ML distance matrix is lower triangular
            if i > consensus_idx:
                dist_line = lines[i + 1].strip().split()
                if len(dist_line) > consensus_idx + 1:
                    ml_dist_to_consensus.append(float(dist_line[consensus_idx + 1]))
            elif i < consensus_idx:
                dist_line = lines[consensus_idx + 1].strip().split()
                if len(dist_line) > i + 1:
                    ml_dist_to_consensus.append(float(dist_line[i + 1]))
    
    print(f"  Extracted {len(ml_dist_to_consensus)} ML distances")

except Exception as e:
    print(f"Warning: Could not parse ML distances: {e}")
    ml_dist_to_consensus = [np.nan] * len(p_distances)

# Ensure arrays are same length
min_len = min(len(p_distances), len(ml_dist_to_consensus), len(strain_ids))
p_distances = p_distances[:min_len]
ml_dist_to_consensus = ml_dist_to_consensus[:min_len]
strain_ids = strain_ids[:min_len]

# Write distances to CSV
with open("${OUTPUT_CSV}", 'w') as f:
    f.write("Strain_ID,P_Distance,ML_Distance_LG_G4,Difference\\n")
    for sid, pd, mld in zip(strain_ids, p_distances, ml_dist_to_consensus):
        diff = mld - pd if not np.isnan(mld) else np.nan
        f.write(f"{sid},{pd:.6f},{mld:.6f},{diff:.6f}\\n")

# Calculate summary statistics
p_mean = np.mean(p_distances)
p_median = np.median(p_distances)
p_std = np.std(p_distances)

ml_valid = [x for x in ml_dist_to_consensus if not np.isnan(x)]
if ml_valid:
    ml_mean = np.mean(ml_valid)
    ml_median = np.median(ml_valid)
    ml_std = np.std(ml_valid)
else:
    ml_mean = ml_median = ml_std = np.nan

print(f"  P-distance: {p_mean:.6f} ± {p_std:.6f} (median: {p_median:.6f})")
print(f"  ML distance: {ml_mean:.6f} ± {ml_std:.6f} (median: {ml_median:.6f})")

# Append to summary
with open("${SUMMARY}", 'a') as f:
    f.write(f"${YEAR},{len(p_distances)},{p_mean:.6f},{p_median:.6f},{p_std:.6f},")
    f.write(f"{ml_mean:.6f},{ml_median:.6f},{ml_std:.6f}\\n")

print(f"  ✅ Distances saved to ${OUTPUT_CSV}")

EOF

    if [ $? -eq 0 ]; then
        echo "✅ ${YEAR} complete"
    else
        echo "❌ ${YEAR} failed"
    fi
    
done

echo ""
echo "=========================================="
echo "✅ Distance analysis complete for ${LINEAGE}!"
echo "=========================================="
echo ""
echo "Summary:"
cat ${SUMMARY}
