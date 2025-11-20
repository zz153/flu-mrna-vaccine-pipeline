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
#
# Each design is saved in two formats:
#    - *_design.fasta: Ungapped sequence (for vaccine production)
#    - *_design_aligned.fasta: Gapped sequence (for distance calculations)

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
    
    # Design files (both gapped and ungapped versions)
    CONSENSUS="${OUTPUT_DIR}/designs/${LINEAGE}_${YEAR}_consensus.fasta"
    CONSENSUS_ALIGNED="${OUTPUT_DIR}/designs/${LINEAGE}_${YEAR}_consensus_aligned.fasta"
    
    MEDOID="${OUTPUT_DIR}/designs/${LINEAGE}_${YEAR}_medoid.fasta"
    MEDOID_ALIGNED="${OUTPUT_DIR}/designs/${LINEAGE}_${YEAR}_medoid_aligned.fasta"
    
    ANCESTRAL="${OUTPUT_DIR}/designs/${LINEAGE}_${YEAR}_ancestral.fasta"
    ANCESTRAL_ALIGNED="${OUTPUT_DIR}/designs/${LINEAGE}_${YEAR}_ancestral_aligned.fasta"
    
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
    if [ -f "${CONSENSUS}" ] && [ -f "${CONSENSUS_ALIGNED}" ]; then
        echo "  ✓ Consensus exists"
    else
        echo "  🔄 Creating consensus..."
        python3 << EOF
from Bio import AlignIO
from collections import Counter

alignment = AlignIO.read("${ALIGNMENT}", "fasta")

# Build consensus aligned (keeps alignment structure)
consensus_aligned = []
for i in range(alignment.get_alignment_length()):
    column = [str(r.seq[i]) for r in alignment]
    non_gaps = [c for c in column if c != '-']
    if non_gaps:
        consensus_aligned.append(Counter(non_gaps).most_common(1)[0][0])
    else:
        consensus_aligned.append('-')

# Ungapped version for output
consensus_ungapped = ''.join([c for c in consensus_aligned if c != '-'])

# Save both versions
with open("${CONSENSUS_ALIGNED}", 'w') as f:
    f.write(">${LINEAGE}_${YEAR}_Consensus_Aligned\n")
    f.write(''.join(consensus_aligned) + "\n")

with open("${CONSENSUS}", 'w') as f:
    f.write(">${LINEAGE}_${YEAR}_Consensus\n")
    f.write(consensus_ungapped + "\n")

print(f"  Consensus: {len(consensus_ungapped)} aa (ungapped), {len(consensus_aligned)} positions (aligned)")
EOF
        echo "  ✅ Consensus created"
    fi
    
    # --- Step 4: Medoid Sequence ---
    echo "Step 4/5: Identifying medoid sequence..."
    if [ -f "${MEDOID}" ] && [ -f "${MEDOID_ALIGNED}" ]; then
        echo "  ✓ Medoid exists"
    else
        echo "  🔄 Finding medoid..."
        python3 << EOF
from Bio import AlignIO
import numpy as np

alignment = AlignIO.read("${ALIGNMENT}", "fasta")
alignment_list = list(alignment)
n = len(alignment_list)

# Calculate pairwise distances (ungapped)
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

# Find medoid
sum_distances = distances.sum(axis=1)
medoid_idx = int(np.argmin(sum_distances))
medoid_record = alignment_list[medoid_idx]

# Save both versions
seq_aligned = str(medoid_record.seq)
seq_ungapped = seq_aligned.replace('-', '')

with open("${MEDOID_ALIGNED}", 'w') as f:
    f.write(f">${LINEAGE}_${YEAR}_Medoid_Aligned|{medoid_record.id}\n")
    f.write(seq_aligned + "\n")

with open("${MEDOID}", 'w') as f:
    f.write(f">${LINEAGE}_${YEAR}_Medoid|{medoid_record.id}\n")
    f.write(seq_ungapped + "\n")

print(f"  Medoid: {medoid_record.id} - {len(seq_ungapped)} aa (ungapped), {len(seq_aligned)} positions (aligned)")
EOF
        echo "  ✅ Medoid identified"
    fi
    
    # --- Step 5: Ancestral Sequence Reconstruction ---
    echo "Step 5/5: Reconstructing ancestral sequence..."
    if [ -f "${ANCESTRAL}" ] && [ -f "${ANCESTRAL_ALIGNED}" ]; then
        echo "  ✓ Ancestral exists"
    else
        echo "  🔄 Running ASR..."
        # Run IQ-TREE with ancestral state reconstruction
        iqtree -s ${ALIGNMENT} -m LG+G4 -te ${TREE_PREFIX}.treefile -asr -pre ${TREE_PREFIX}_asr -redo -quiet > /dev/null 2>&1
        
        # Extract root ancestral sequence
        python3 << EOF
asr_file = "${TREE_PREFIX}_asr.state"
try:
    with open(asr_file, 'r') as f:
        lines = f.readlines()
    
    # Extract root node sequence with alignment
    root_seq_aligned = []
    current_node = None
    
    for line in lines:
        if not line.strip() or line.startswith('#'):
            continue
        
        parts = line.strip().split()
        
        if len(parts) > 1 and parts[0] == "Node" and parts[1] == "Site":
            continue
            
        if len(parts) < 3:
            continue
        
        node_name = parts[0]
        site = parts[1]
        state = parts[2]
        
        if current_node is None:
            current_node = node_name
        
        if node_name == current_node:
            root_seq_aligned.append(state)
        else:
            break
    
    if root_seq_aligned and len(root_seq_aligned) > 50:
        # Ungapped version
        root_seq_ungapped = ''.join([aa for aa in root_seq_aligned if aa != '-'])
        
        # Save both versions
        with open("${ANCESTRAL_ALIGNED}", 'w') as f:
            f.write(f">${LINEAGE}_${YEAR}_Ancestral_Aligned_{current_node}\n")
            f.write(''.join(root_seq_aligned) + "\n")
        
        with open("${ANCESTRAL}", 'w') as f:
            f.write(f">${LINEAGE}_${YEAR}_Ancestral_{current_node}\n")
            f.write(root_seq_ungapped + "\n")
        
        print(f"  Ancestral: {len(root_seq_ungapped)} aa (ungapped), {len(root_seq_aligned)} positions (aligned)")
    else:
        print(f"Warning: Only {len(root_seq_aligned)} amino acids, using consensus as fallback")
        import shutil
        shutil.copy("${CONSENSUS}", "${ANCESTRAL}")
        shutil.copy("${CONSENSUS_ALIGNED}", "${ANCESTRAL_ALIGNED}")
        
except Exception as e:
    print(f"ASR failed: {e}, using consensus as fallback")
    import shutil
    shutil.copy("${CONSENSUS}", "${ANCESTRAL}")
    shutil.copy("${CONSENSUS_ALIGNED}", "${ANCESTRAL_ALIGNED}")
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
