#!/bin/bash
# Script 04: Per-Year Analysis Pipeline with Multiple Vaccine Designs
#
# For each year (2009-2025), this script:
# 1. Aligns sequences with MAFFT
# 2. Builds initial phylogenetic tree with IQ-TREE
# 3. Creates vaccine designs:
#    - Consensus (most common amino acid at each position)
#    - Medoid (actual strain with minimum distance to all others)
#    - Ancestral (reconstructed ancestral sequence from tree root)
#    - COBRA (2-round clustering approach for broad reactivity)
# 4. Re-aligns COBRA to match original alignment gap structure
# 5. Adds all designs to alignment
# 6. Rebuilds tree with designs included
# 7. Extracts ML distances for all designs
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

# Check for required tools
command -v cd-hit >/dev/null 2>&1 || { echo "❌ cd-hit not found. Required for COBRA generation."; exit 1; }

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
mkdir -p ${OUTPUT_DIR}/{alignments,trees,designs,logs}

# Summary file
SUMMARY="${OUTPUT_DIR}/analysis_summary.csv"
echo "Year,Sequences,Alignment,Tree,Consensus,Medoid,Ancestral,COBRA,Tree_with_Designs" > ${SUMMARY}

for YEAR in $(seq ${YEAR_START} ${YEAR_END}); do
    
    INPUT_FILE="${INPUT_DIR}/${LINEAGE}_${YEAR}.fasta"
    
    if [ ! -f "$INPUT_FILE" ]; then
        continue
    fi
    
    SEQ_COUNT=$(grep -c "^>" ${INPUT_FILE} || echo "0")
    
    if [ "$SEQ_COUNT" -lt $MIN_SEQUENCES ]; then
        echo "⏭️  ${YEAR}: Too few sequences (${SEQ_COUNT})"
        echo "${YEAR},${SEQ_COUNT},SKIPPED,SKIPPED,SKIPPED,SKIPPED,SKIPPED,SKIPPED,SKIPPED" >> ${SUMMARY}
        continue
    fi
    
    echo ""
    echo "=========================================="
    echo "📅 Processing ${LINEAGE} ${YEAR} (${SEQ_COUNT} sequences)"
    echo "=========================================="
    
    ALIGNMENT="${OUTPUT_DIR}/alignments/${LINEAGE}_${YEAR}_aligned.fasta"
    TREE_PREFIX="${OUTPUT_DIR}/trees/${LINEAGE}_${YEAR}"
    LOG_DIR="${OUTPUT_DIR}/logs/${YEAR}"
    mkdir -p ${LOG_DIR}
    
    # Design files (both gapped and ungapped versions)
    CONSENSUS="${OUTPUT_DIR}/designs/${LINEAGE}_${YEAR}_consensus.fasta"
    CONSENSUS_ALIGNED="${OUTPUT_DIR}/designs/${LINEAGE}_${YEAR}_consensus_aligned.fasta"
    
    MEDOID="${OUTPUT_DIR}/designs/${LINEAGE}_${YEAR}_medoid.fasta"
    MEDOID_ALIGNED="${OUTPUT_DIR}/designs/${LINEAGE}_${YEAR}_medoid_aligned.fasta"
    
    ANCESTRAL="${OUTPUT_DIR}/designs/${LINEAGE}_${YEAR}_ancestral.fasta"
    ANCESTRAL_ALIGNED="${OUTPUT_DIR}/designs/${LINEAGE}_${YEAR}_ancestral_aligned.fasta"
    
    COBRA="${OUTPUT_DIR}/designs/${LINEAGE}_${YEAR}_cobra.fasta"
    COBRA_ALIGNED="${OUTPUT_DIR}/designs/${LINEAGE}_${YEAR}_cobra_aligned.fasta"
    
    # --- Step 1: Alignment ---
    echo "Step 1/8: Aligning sequences..."
    if [ -f "${ALIGNMENT}" ]; then
        echo "  ✓ Alignment exists"
    else
        echo "  🔄 Aligning..."
        mafft --auto --quiet ${INPUT_FILE} > ${ALIGNMENT} 2>/dev/null
        echo "  ✅ Aligned"
    fi
    
    # --- Step 2: Build Initial Tree ---
    echo "Step 2/8: Building initial phylogenetic tree..."
    if [ -f "${TREE_PREFIX}.treefile" ]; then
        echo "  ✓ Tree exists"
    else
        echo "  🔄 Building tree..."
        iqtree -s ${ALIGNMENT} -m LG+G4 -bb 1000 -nt ${THREADS} -pre ${TREE_PREFIX} -redo -quiet
        echo "  ✅ Tree built"
    fi
    
    # --- Step 3: Consensus Sequence ---
    echo "Step 3/8: Creating consensus sequence..."
    if [ -f "${CONSENSUS}" ] && [ -f "${CONSENSUS_ALIGNED}" ]; then
        echo "  ✓ Consensus exists"
    else
        echo "  🔄 Creating consensus..."
        python3 << EOF
from Bio import AlignIO
from collections import Counter

alignment = AlignIO.read("${ALIGNMENT}", "fasta")

consensus_aligned = []
for i in range(alignment.get_alignment_length()):
    column = [str(r.seq[i]) for r in alignment]
    non_gaps = [c for c in column if c != '-']
    if non_gaps:
        consensus_aligned.append(Counter(non_gaps).most_common(1)[0][0])
    else:
        consensus_aligned.append('-')

consensus_ungapped = ''.join([c for c in consensus_aligned if c != '-'])

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
    echo "Step 4/8: Identifying medoid sequence..."
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

sum_distances = distances.sum(axis=1)
medoid_idx = int(np.argmin(sum_distances))
medoid_record = alignment_list[medoid_idx]

seq_aligned = str(medoid_record.seq)
seq_ungapped = seq_aligned.replace('-', '')

with open("${MEDOID_ALIGNED}", 'w') as f:
    f.write(f">${LINEAGE}_${YEAR}_Medoid_Aligned|{medoid_record.id}\n")
    f.write(seq_aligned + "\n")

with open("${MEDOID}", 'w') as f:
    f.write(f">${LINEAGE}_${YEAR}_Medoid|{medoid_record.id}\n")
    f.write(seq_ungapped + "\n")

# Save medoid index for later ML distance extraction
with open("${LOG_DIR}/medoid_idx.txt", 'w') as f:
    f.write(str(medoid_idx))

print(f"  Medoid: {medoid_record.id} - {len(seq_ungapped)} aa (ungapped)")
EOF
        echo "  ✅ Medoid identified"
    fi
    
    # --- Step 5: Ancestral Sequence Reconstruction ---
    echo "Step 5/8: Reconstructing ancestral sequence..."
    if [ -f "${ANCESTRAL}" ] && [ -f "${ANCESTRAL_ALIGNED}" ]; then
        echo "  ✓ Ancestral exists"
    else
        echo "  🔄 Running ASR..."
        iqtree -s ${ALIGNMENT} -m LG+G4 -te ${TREE_PREFIX}.treefile -asr -pre ${TREE_PREFIX}_asr -redo -quiet > /dev/null 2>&1
        
        python3 << EOF
asr_file = "${TREE_PREFIX}_asr.state"
try:
    with open(asr_file, 'r') as f:
        lines = f.readlines()
    
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
        state = parts[2]
        
        if current_node is None:
            current_node = node_name
        
        if node_name == current_node:
            root_seq_aligned.append(state)
        else:
            break
    
    if root_seq_aligned and len(root_seq_aligned) > 50:
        root_seq_ungapped = ''.join([aa for aa in root_seq_aligned if aa != '-'])
        
        with open("${ANCESTRAL_ALIGNED}", 'w') as f:
            f.write(f">${LINEAGE}_${YEAR}_Ancestral_Aligned_{current_node}\n")
            f.write(''.join(root_seq_aligned) + "\n")
        
        with open("${ANCESTRAL}", 'w') as f:
            f.write(f">${LINEAGE}_${YEAR}_Ancestral_{current_node}\n")
            f.write(root_seq_ungapped + "\n")
        
        print(f"  Ancestral: {len(root_seq_ungapped)} aa from {current_node}")
    else:
        print(f"Warning: ASR insufficient, using consensus")
        import shutil
        shutil.copy("${CONSENSUS}", "${ANCESTRAL}")
        shutil.copy("${CONSENSUS_ALIGNED}", "${ANCESTRAL_ALIGNED}")
        
except Exception as e:
    print(f"ASR failed: {e}, using consensus")
    import shutil
    shutil.copy("${CONSENSUS}", "${ANCESTRAL}")
    shutil.copy("${CONSENSUS_ALIGNED}", "${ANCESTRAL_ALIGNED}")
EOF
        echo "  ✅ Ancestral reconstructed"
    fi
    
    # --- Step 6: COBRA Sequence (2-Round Clustering) ---
    echo "Step 6/8: Generating COBRA sequence..."
    if [ -f "${COBRA}" ] && [ -f "${COBRA_ALIGNED}" ]; then
        echo "  ✓ COBRA exists"
    else
        echo "  🔄 Running 2-round COBRA clustering..."
        
        # COBRA Round 1: Cluster at 95%
        cd-hit -i ${INPUT_FILE} -o ${LOG_DIR}/cobra_r1.fasta \
            -c 0.95 -n 5 -M 0 -T ${THREADS} \
            > ${LOG_DIR}/cdhit_r1.log 2>&1 || {
            echo "  ⚠️  CD-HIT failed, using consensus as COBRA"
            cp ${CONSENSUS} ${COBRA}
            cp ${CONSENSUS_ALIGNED} ${COBRA_ALIGNED}
            sed -i "s/>.*/>$LINEAGE}_${YEAR}_COBRA/" ${COBRA}
            sed -i "s/>.*/>$LINEAGE}_${YEAR}_COBRA_Aligned/" ${COBRA_ALIGNED}
            echo "  ✅ COBRA (fallback to consensus)"
            echo "${YEAR},${SEQ_COUNT},✓,✓,✓,✓,✓,✓(fallback),SKIPPED" >> ${SUMMARY}
            continue
        }
        
        N_R1=$(grep -c '^>' ${LOG_DIR}/cobra_r1.fasta || echo 0)
        echo "    Round 1: ${N_R1} clusters"
        
        if [ "$N_R1" -lt 2 ]; then
            echo "  ⚠️  Insufficient diversity, using round 1 consensus"
            # Align and create consensus from round 1
            mafft --quiet ${LOG_DIR}/cobra_r1.fasta > ${LOG_DIR}/cobra_r1_aln.fasta 2>/dev/null
            
            python3 << EOF
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter

aln = AlignIO.read("${LOG_DIR}/cobra_r1_aln.fasta", "fasta")
cons = []
for i in range(aln.get_alignment_length()):
    col = [str(r.seq)[i] for r in aln]
    c = Counter(col)
    if "-" in c:
        c["-"] = int(c["-"] * 0.2)
    best = max(c, key=c.get)
    cons.append("X" if best == "-" else best)

cons_str = ''.join(cons)
cons_ungapped = cons_str.replace('X', '').replace('-', '')

# Save ungapped version
record_ungapped = SeqRecord(Seq(cons_ungapped), id="${LINEAGE}_${YEAR}_COBRA", description="")
SeqIO.write([record_ungapped], "${COBRA}", "fasta")
print(f"  COBRA ungapped: {len(cons_ungapped)} aa")
EOF
        else
            # COBRA Round 2: Cluster round 1 at 90%
            cd-hit -i ${LOG_DIR}/cobra_r1.fasta -o ${LOG_DIR}/cobra_r2.fasta \
                -c 0.90 -n 5 -M 0 -T ${THREADS} \
                > ${LOG_DIR}/cdhit_r2.log 2>&1 || {
                echo "  ⚠️  Round 2 failed, using round 1"
                mafft --quiet ${LOG_DIR}/cobra_r1.fasta > ${LOG_DIR}/cobra_r1_aln.fasta 2>/dev/null
                
                python3 << EOF
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter

aln = AlignIO.read("${LOG_DIR}/cobra_r1_aln.fasta", "fasta")
cons = []
for i in range(aln.get_alignment_length()):
    col = [str(r.seq)[i] for r in aln]
    c = Counter(col)
    if "-" in c:
        c["-"] = int(c["-"] * 0.2)
    best = max(c, key=c.get)
    cons.append("X" if best == "-" else best)

cons_str = ''.join(cons)
cons_ungapped = cons_str.replace('X', '').replace('-', '')

record_ungapped = SeqRecord(Seq(cons_ungapped), id="${LINEAGE}_${YEAR}_COBRA", description="")
SeqIO.write([record_ungapped], "${COBRA}", "fasta")
print(f"  COBRA ungapped: {len(cons_ungapped)} aa")
EOF
            }
            
            N_R2=$(grep -c '^>' ${LOG_DIR}/cobra_r2.fasta 2>/dev/null || echo 0)
            if [ "$N_R2" -gt 0 ]; then
                echo "    Round 2: ${N_R2} clusters"
                
                # Align round 2 and create final consensus
                mafft --quiet ${LOG_DIR}/cobra_r2.fasta > ${LOG_DIR}/cobra_r2_aln.fasta 2>/dev/null
                
                python3 << EOF
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter

aln = AlignIO.read("${LOG_DIR}/cobra_r2_aln.fasta", "fasta")
cons = []
for i in range(aln.get_alignment_length()):
    col = [str(r.seq)[i] for r in aln]
    c = Counter(col)
    if "-" in c:
        c["-"] = int(c["-"] * 0.2)
    best = max(c, key=c.get)
    cons.append("X" if best == "-" else best)

cons_str = ''.join(cons)
cons_ungapped = cons_str.replace('X', '').replace('-', '')

record_ungapped = SeqRecord(Seq(cons_ungapped), id="${LINEAGE}_${YEAR}_COBRA", description="")
SeqIO.write([record_ungapped], "${COBRA}", "fasta")
print(f"  COBRA ungapped: {len(cons_ungapped)} aa (2-round)")
EOF
            fi
        fi
        
        # --- CRITICAL: Re-align COBRA to match original alignment ---
        echo "  🔄 Re-aligning COBRA to match original alignment..."
        mafft --add ${COBRA} --keeplength --quiet ${ALIGNMENT} > ${LOG_DIR}/cobra_realigned.fasta 2>/dev/null
        
        # Extract the COBRA aligned sequence (it will be the last sequence in the output)
        python3 << EOF
from Bio import SeqIO

records = list(SeqIO.parse("${LOG_DIR}/cobra_realigned.fasta", "fasta"))
# COBRA should be the last record (added with --add)
cobra_record = records[-1]

# Update ID to match our naming
cobra_record.id = "${LINEAGE}_${YEAR}_COBRA_Aligned"
cobra_record.description = ""

# Save aligned version
SeqIO.write([cobra_record], "${COBRA_ALIGNED}", "fasta")
print(f"  COBRA aligned: {len(str(cobra_record.seq))} positions")
EOF
        
        echo "  ✅ COBRA created and aligned"
    fi
    
    # --- Step 7: Verify All Designs Have Same Length ---
    echo "Step 7/8: Verifying design alignment lengths..."
    python3 << EOF
from Bio import SeqIO

# Get alignment length from first sequence
first_record = next(SeqIO.parse("${ALIGNMENT}", "fasta"))
alignment_len = len(str(first_record.seq))

designs = ["consensus", "medoid", "ancestral", "cobra"]
all_good = True

for design in designs:
    design_file = "${OUTPUT_DIR}/designs/${LINEAGE}_${YEAR}_" + design + "_aligned.fasta"
    try:
        record = SeqIO.read(design_file, "fasta")
        design_len = len(str(record.seq))
        if design_len != alignment_len:
            print(f"  ❌ {design}: {design_len} (expected {alignment_len})")
            all_good = False
        else:
            print(f"  ✓ {design}: {design_len}")
    except:
        print(f"  ❌ {design}: file not found")
        all_good = False

if all_good:
    print("  ✅ All designs have correct alignment length")
else:
    print("  ❌ Alignment length mismatch detected!")
    import sys
    sys.exit(1)
EOF
        
    if [ $? -ne 0 ]; then
        echo "  ❌ Skipping tree building due to alignment errors"
        echo "${YEAR},${SEQ_COUNT},✓,✓,✓,✓,✓,✓,ERROR" >> ${SUMMARY}
        continue
    fi
    
    # --- Step 8: Rebuild Tree with All Designs ---
    echo "Step 8/8: Rebuilding tree with all designs..."
    
    ALIGNMENT_WITH_DESIGNS="${OUTPUT_DIR}/alignments/${LINEAGE}_${YEAR}_with_designs.fasta"
    TREE_WITH_DESIGNS="${OUTPUT_DIR}/trees/${LINEAGE}_${YEAR}_with_designs"
    
    if [ -f "${TREE_WITH_DESIGNS}.treefile" ]; then
        echo "  ✓ Tree with designs exists"
    else
        echo "  🔄 Adding designs to alignment..."
        
        # Combine original alignment with all designs (aligned versions)
        cat ${ALIGNMENT} \
            ${CONSENSUS_ALIGNED} \
            ${MEDOID_ALIGNED} \
            ${ANCESTRAL_ALIGNED} \
            ${COBRA_ALIGNED} \
            > ${ALIGNMENT_WITH_DESIGNS}
        
        echo "  🔄 Building tree with designs..."
        iqtree -s ${ALIGNMENT_WITH_DESIGNS} -m LG+G4 -bb 1000 -nt ${THREADS} \
            -pre ${TREE_WITH_DESIGNS} -redo -quiet
        
        echo "  ✅ Tree with designs built"
    fi
    
    echo "${YEAR},${SEQ_COUNT},✓,✓,✓,✓,✓,✓,✓" >> ${SUMMARY}
    echo "✅ ${YEAR} complete (all 4 designs + tree with designs)"
    
done

echo ""
echo "=========================================="
echo "✅ Analysis complete for ${LINEAGE}!"
echo "=========================================="
echo ""
echo "Summary:"
cat ${SUMMARY}
