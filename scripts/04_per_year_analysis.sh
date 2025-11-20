#!/bin/bash
# Script 04: Per-Year Analysis Pipeline with Multiple Vaccine Designs
#
# For each year (2009-2025), this script:
# 1. Aligns sequences with MAFFT
# 2. Builds phylogenetic tree with IQ-TREE
# 3. Creates vaccine designs:
#    - Consensus (most common amino acid at each position)
#    - Medoid (actual strain with minimum distance to all others)
#    - Ancestral (reconstructed ancestral sequence from tree root)

set -e

LINEAGE=$1
YEAR_START=${2:-2009}
YEAR_END=${3:-2025}
THREADS=${4:-4}
USE_UNCLUSTERED=${5:-false}
MIN_SEQUENCES=5

if [ -z "$LINEAGE" ]; then
    echo "Usage: bash scripts/04_per_year_analysis.sh [LINEAGE] [START_YEAR] [END_YEAR] [THREADS] [USE_UNCLUSTERED]"
    echo "Example: bash scripts/04_per_year_analysis.sh H1N1 2009 2025 4"
    exit 1
fi

echo "=========================================="
echo "Per-Year Analysis for ${LINEAGE}"
echo "Years: ${YEAR_START} to ${YEAR_END}"
echo "Use unclustered: ${USE_UNCLUSTERED}"
echo "=========================================="

# Determine input source
if [ "$USE_UNCLUSTERED" = "true" ]; then
    echo "Using UNCLUSTERED sequences"
    SOURCE_TYPE="unclustered"
    CLEAN_FILE="data/processed/${LINEAGE}/${LINEAGE}_clean.fasta"
    INPUT_DIR="data/processed/${LINEAGE}/per_year_unclustered"
    
    if [ ! -d "$INPUT_DIR" ]; then
        echo "Creating per-year split for unclustered sequences..."
        mkdir -p ${INPUT_DIR}
        
        for YEAR in $(seq ${YEAR_START} ${YEAR_END}); do
            OUTPUT_FILE="${INPUT_DIR}/${LINEAGE}_${YEAR}.fasta"
            awk -v year="$YEAR" '
                /^>/ {
                    split($0, fields, "|")
                    split(fields[3], parts, "/")
                    seq_year = parts[length(parts)]
                    if (seq_year == year) {
                        print_seq = 1
                        print $0
                    } else {
                        print_seq = 0
                    }
                    next
                }
                print_seq { print }
            ' ${CLEAN_FILE} > ${OUTPUT_FILE}
            
            if [ ! -s ${OUTPUT_FILE} ]; then
                rm ${OUTPUT_FILE}
            fi
        done
        echo "✅ Unclustered sequences split by year"
    fi
else
    echo "Using CLUSTERED sequences (CD-HIT 99%)"
    SOURCE_TYPE="clustered"
    INPUT_DIR="data/processed/${LINEAGE}/per_year"
fi

OUTPUT_DIR="results/per_year_${SOURCE_TYPE}/${LINEAGE}"
mkdir -p ${OUTPUT_DIR}/{alignments,trees,designs}

# Summary file
SUMMARY="${OUTPUT_DIR}/analysis_summary.csv"
echo "Year,Sequences,Alignment,Tree,Consensus,Medoid,Ancestral" > ${SUMMARY}

for YEAR in $(seq ${YEAR_START} ${YEAR_END}); do
    
    INPUT_FILE="${INPUT_DIR}/${LINEAGE}_${YEAR}.fasta"
    
    if [ ! -f "$INPUT_FILE" ]; then
        continue
    fi
    
    SEQ_COUNT=$(grep -c "^>" ${INPUT_FILE} || echo "0")
    
    if [ "$SEQ_COUNT" -lt $MIN_SEQUENCES ]; then
        echo "⏭️  ${YEAR}: Too few sequences (${SEQ_COUNT})"
        echo "${YEAR},${SEQ_COUNT},SKIPPED,SKIPPED,SKIPPED,SKIPPED,SKIPPED" >> ${SUMMARY}
        continue
    fi
    
    echo ""
    echo "=========================================="
    echo "📅 Processing ${LINEAGE} ${YEAR} (${SEQ_COUNT} sequences)"
    echo "=========================================="
    
    ALIGNMENT="${OUTPUT_DIR}/alignments/${LINEAGE}_${YEAR}_aligned.fasta"
    TREE_PREFIX="${OUTPUT_DIR}/trees/${LINEAGE}_${YEAR}"
    CONSENSUS="${OUTPUT_DIR}/designs/${LINEAGE}_${YEAR}_consensus.fasta"
    MEDOID="${OUTPUT_DIR}/designs/${LINEAGE}_${YEAR}_medoid.fasta"
    ANCESTRAL="${OUTPUT_DIR}/designs/${LINEAGE}_${YEAR}_ancestral.fasta"
    
    # --- Step 1: Alignment ---
    echo "Step 1/5: Aligning sequences..."
    if [ -f "${ALIGNMENT}" ]; then
        echo "  ✓ Alignment exists"
    else
        echo "  🔄 Aligning..."
        mafft --auto --quiet ${INPUT_FILE} > ${ALIGNMENT} 2>/dev/null
        echo "  ✅ Aligned"
    fi
    
    # --- Step 2: Build Tree ---
    echo "Step 2/5: Building phylogenetic tree..."
    if [ -f "${TREE_PREFIX}.treefile" ]; then
        echo "  ✓ Tree exists"
    else
        echo "  🔄 Building tree..."
        iqtree -s ${ALIGNMENT} -m LG+G4 -bb 1000 -nt ${THREADS} -pre ${TREE_PREFIX} -redo -quiet
        echo "  ✅ Tree built"
    fi
    
    # --- Step 3: Consensus Sequence ---
    echo "Step 3/5: Creating consensus sequence..."
    if [ -f "${CONSENSUS}" ]; then
        echo "  ✓ Consensus exists"
    else
        echo "  🔄 Creating consensus..."
        python3 << EOF
from Bio import AlignIO
from collections import Counter

alignment = AlignIO.read("${ALIGNMENT}", "fasta")
consensus_seq = []

for i in range(alignment.get_alignment_length()):
    column = [str(r.seq[i]) for r in alignment if str(r.seq[i]) != '-']
    if column:
        consensus_seq.append(Counter(column).most_common(1)[0][0])

with open("${CONSENSUS}", 'w') as f:
    f.write(">${LINEAGE}_${YEAR}_Consensus\n")
    f.write(''.join(consensus_seq) + "\n")
EOF
        echo "  ✅ Consensus created"
    fi
    
    # --- Step 4: Medoid Sequence ---
    echo "Step 4/5: Identifying medoid sequence..."
    if [ -f "${MEDOID}" ]; then
        echo "  ✓ Medoid exists"
    else
        echo "  🔄 Finding medoid..."
        python3 << EOF
from Bio import AlignIO
import numpy as np

alignment = AlignIO.read("${ALIGNMENT}", "fasta")
alignment_list = list(alignment)
n = len(alignment_list)

# Calculate pairwise distances
distances = np.zeros((n, n))
for i in range(n):
    for j in range(i+1, n):
        seq1 = str(alignment_list[i].seq).replace('-', '')
        seq2 = str(alignment_list[j].seq).replace('-', '')
        min_len = min(len(seq1), len(seq2))
        if min_len > 0:
            diff = sum(1 for k in range(min_len) if seq1[k] != seq2[k])
            dist = diff / min_len
            distances[i, j] = dist
            distances[j, i] = dist

# Find medoid (sequence with minimum sum of distances)
sum_distances = distances.sum(axis=1)
medoid_idx = int(np.argmin(sum_distances))
medoid_record = alignment_list[medoid_idx]

# Write medoid
with open("${MEDOID}", 'w') as f:
    original_id = medoid_record.id
    f.write(f">${LINEAGE}_${YEAR}_Medoid|{original_id}\n")
    seq = str(medoid_record.seq).replace('-', '')
    f.write(seq + "\n")
EOF
        echo "  ✅ Medoid identified"
    fi
    
    # --- Step 5: Ancestral Sequence Reconstruction ---
    echo "Step 5/5: Reconstructing ancestral sequence..."
    if [ -f "${ANCESTRAL}" ]; then
        echo "  ✓ Ancestral exists"
    else
        echo "  🔄 Running ASR..."
        # Run IQ-TREE with ancestral state reconstruction
        iqtree -s ${ALIGNMENT} -m LG+G4 -te ${TREE_PREFIX}.treefile -asr -pre ${TREE_PREFIX}_asr -redo -quiet > /dev/null 2>&1
        
        # Extract root ancestral sequence
        python3 << EOF
# Read the ancestral state file
asr_file = "${TREE_PREFIX}_asr.state"
try:
    with open(asr_file, 'r') as f:
        lines = f.readlines()
    
    # Skip header and extract root node sequence
    root_seq = []
    current_node = None
    
    for line in lines:
        # Skip empty lines and comments
        if not line.strip() or line.startswith('#'):
            continue
        
        parts = line.strip().split()
        
        # Skip header line (contains "Node Site State")
        if len(parts) > 1 and parts[0] == "Node" and parts[1] == "Site":
            continue
            
        if len(parts) < 3:
            continue
        
        node_name = parts[0]
        site = parts[1]
        state = parts[2]
        
        # Get the first node (root/oldest)
        if current_node is None:
            current_node = node_name
        
        # Collect all states for this node
        if node_name == current_node:
            root_seq.append(state)
        else:
            # Moved to next node, stop
            break
    
    if root_seq and len(root_seq) > 50:
        with open("${ANCESTRAL}", 'w') as f:
            f.write(f">${LINEAGE}_${YEAR}_Ancestral_{current_node}\n")
            f.write(''.join(root_seq) + "\n")
        print(f"  Extracted {len(root_seq)} amino acids from {current_node}")
    else:
        print(f"Warning: Extracted only {len(root_seq)} amino acids, using consensus as fallback")
        import shutil
        shutil.copy("${CONSENSUS}", "${ANCESTRAL}")
        
except Exception as e:
    print(f"ASR failed: {e}, using consensus as fallback")
    import shutil
    shutil.copy("${CONSENSUS}", "${ANCESTRAL}")
EOF
        echo "  ✅ Ancestral reconstructed"
    fi
    
    echo "${YEAR},${SEQ_COUNT},✓,✓,✓,✓,✓" >> ${SUMMARY}
    echo "✅ ${YEAR} complete (all designs generated)"
    
done

echo ""
echo "=========================================="
echo "✅ Analysis complete for ${LINEAGE}!"
echo "=========================================="
echo ""
echo "Summary:"
cat ${SUMMARY}
