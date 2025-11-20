#!/bin/bash
# Script 04: Per-Year Analysis Pipeline

set -e

LINEAGE=$1
YEAR_START=${2:-2009}
YEAR_END=${3:-2025}
THREADS=${4:-4}
MIN_SEQUENCES=5

if [ -z "$LINEAGE" ]; then
    echo "Usage: bash scripts/04_per_year_analysis.sh [LINEAGE] [START_YEAR] [END_YEAR] [THREADS]"
    echo "Example: bash scripts/04_per_year_analysis.sh H1N1 2009 2025 16"
    exit 1
fi

echo "=========================================="
echo "Per-Year Analysis for ${LINEAGE}"
echo "Years: ${YEAR_START} to ${YEAR_END}"
echo "=========================================="

INPUT_DIR="data/processed/${LINEAGE}/per_year"
OUTPUT_DIR="results/per_year_clustered/${LINEAGE}"
mkdir -p ${OUTPUT_DIR}/{alignments,trees,designs}

SUMMARY="${OUTPUT_DIR}/analysis_summary.csv"
echo "Year,Sequences,Status" > ${SUMMARY}

for YEAR in $(seq ${YEAR_START} ${YEAR_END}); do
    INPUT_FILE="${INPUT_DIR}/${LINEAGE}_${YEAR}.fasta"
    
    if [ ! -f "$INPUT_FILE" ]; then
        continue
    fi
    
    SEQ_COUNT=$(grep -c "^>" ${INPUT_FILE} || echo "0")
    
    if [ "$SEQ_COUNT" -lt $MIN_SEQUENCES ]; then
        echo "⏭️  ${YEAR}: Too few sequences (${SEQ_COUNT})"
        echo "${YEAR},${SEQ_COUNT},SKIPPED" >> ${SUMMARY}
        continue
    fi
    
    echo ""
    echo "📅 Processing ${YEAR} (${SEQ_COUNT} sequences)..."
    
    ALIGNMENT="${OUTPUT_DIR}/alignments/${LINEAGE}_${YEAR}_aligned.fasta"
    TREE_PREFIX="${OUTPUT_DIR}/trees/${LINEAGE}_${YEAR}"
    CONSENSUS="${OUTPUT_DIR}/designs/${LINEAGE}_${YEAR}_consensus.fasta"
    
    # Alignment
    if [ ! -f "${ALIGNMENT}" ]; then
        echo "  🔄 Aligning..."
        mafft --auto --quiet ${INPUT_FILE} > ${ALIGNMENT} 2>/dev/null
        echo "  ✅ Aligned"
    else
        echo "  ✓ Alignment exists"
    fi
    
    # Tree
    if [ ! -f "${TREE_PREFIX}.treefile" ]; then
        echo "  🔄 Building tree..."
        iqtree -s ${ALIGNMENT} -m LG+G4 -bb 1000 -nt ${THREADS} -pre ${TREE_PREFIX} -redo -quiet
        echo "  ✅ Tree built"
    else
        echo "  ✓ Tree exists"
    fi
    
    # Consensus
    if [ ! -f "${CONSENSUS}" ]; then
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
    else
        echo "  ✓ Consensus exists"
    fi
    
    echo "${YEAR},${SEQ_COUNT},COMPLETE" >> ${SUMMARY}
    echo "✅ ${YEAR} complete"
done

echo ""
echo "=========================================="
echo "✅ Analysis complete for ${LINEAGE}!"
echo "=========================================="
echo ""
echo "Summary:"
cat ${SUMMARY}
