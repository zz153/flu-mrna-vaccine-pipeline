#!/bin/bash -e
#SBATCH --job-name=yearly_cobra
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --output=logs/yearly_cobra.%j.out
#SBATCH --error=logs/yearly_cobra.%j.err

# =========================================================
#  Generic COBRA Builder (CD-HIT clustering approach)
#  Usage:
#     sbatch run_yearly_cobra.slurm H1N1 h1n1 [h1n1_per_year_fastas]
#  If INPUT_DIR provided, uses original unaligned fastas from there
# =========================================================

set -euo pipefail

# ---- Environment ----
ENV=/home/ranzo85p/conda_envs/bio
PY=$ENV/bin/python
MAFFT=$ENV/bin/mafft
CDHIT=$ENV/bin/cd-hit

PREFIX=${1:-H1N1}
OUT_DIR=${2:-h1n1}
INPUT_DIR=${3:-""}  # Optional: directory with original unaligned fastas
CPUS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p logs

echo "=== Tool versions ==="
$PY --version
$PY -c "import Bio; print('Biopython', Bio.__version__)"
$CDHIT -h 2>&1 | head -n 1 || echo "CD-HIT installed"
echo "=== Starting COBRA for all years in ${OUT_DIR} ==="

# =========================================================
#  Safe COBRA with fallback logic
# =========================================================
safe_cobra() {
    local input_fasta="$1"
    local year="$2"
    local designdir="$3"
    local logdir="$4"
    
    echo "[INFO] COBRA Round 1: Clustering at 95% identity"
    
    # Round 1: Cluster → Align → Consensus
    $CDHIT -i "$input_fasta" -o "$designdir/HA_r1.fasta" \
        -c 0.95 -n 5 -M 0 -T "$CPUS" \
        > "$logdir/cdhit_r1.log" 2>&1 || {
        echo "[ERROR] CD-HIT round 1 failed for ${year}"
        return 1
    }
    
    local n1=$(grep -c '^>' "$designdir/HA_r1.fasta" || echo 0)
    if (( n1 == 0 )); then
        echo "[ERROR] CD-HIT produced no clusters for ${year}"
        return 1
    fi
    
    echo "[INFO] Round 1: ${n1} cluster representatives"
    
    $MAFFT --auto "$designdir/HA_r1.fasta" \
        > "$designdir/HA_r1_aln.fasta" \
        2> "$logdir/mafft_r1.log" || {
        echo "[ERROR] MAFFT round 1 failed for ${year}"
        return 1
    }
    
    # Build consensus of round 1
    $PY <<EOF
from collections import Counter
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

try:
    aln = AlignIO.read("${designdir}/HA_r1_aln.fasta", "fasta")
except Exception as e:
    import sys
    sys.stderr.write(f"[ERROR] Could not read round 1 alignment: {e}\n")
    sys.exit(1)

cons = []
for i in range(aln.get_alignment_length()):
    col = [str(r.seq)[i] for r in aln]
    c = Counter(col)
    if "-" in c:
        c["-"] = int(c["-"] * 0.2)
    best = max(c, key=c.get)
    cons.append("X" if best == "-" else best)

record = SeqRecord(Seq("".join(cons)), id="HA_R1_CONS", description="")
SeqIO.write([record], "${designdir}/HA_r1_cons.fasta", "fasta")
print("[OK] Round 1 consensus created")
EOF
    
    if [[ ! -s "$designdir/HA_r1_cons.fasta" ]]; then
        echo "[ERROR] Round 1 consensus failed for ${year}"
        return 1
    fi
    
    # If only 1 cluster, use round 1 consensus as final
    if (( n1 < 2 )); then
        echo "[INFO] Only 1 cluster - using round 1 consensus as COBRA"
        cp "$designdir/HA_r1_cons.fasta" "$designdir/${PREFIX}_${year}_COBRA.fasta"
        sed -i "s/>HA_R1_CONS/>${PREFIX}_${year}_COBRA/" "$designdir/${PREFIX}_${year}_COBRA.fasta"
        return 0
    fi
    
    # Round 2: Cluster round 1 reps → Align → Final consensus
    echo "[INFO] COBRA Round 2: Clustering round 1 reps at 90% identity"
    
    $CDHIT -i "$designdir/HA_r1.fasta" -o "$designdir/HA_r2.fasta" \
        -c 0.90 -n 5 -M 0 -T "$CPUS" \
        > "$logdir/cdhit_r2.log" 2>&1 || {
       echo "[WARN] CD-HIT round 2 failed - using round 1 consensus"
        cp "$designdir/HA_r1_cons.fasta" "$designdir/${PREFIX}_${year}_COBRA.fasta"
        sed -i "s/>HA_R1_CONS/>${PREFIX}_${year}_COBRA/" "$designdir/${PREFIX}_${year}_COBRA.fasta"
        return 0
    }
    
    local n2=$(grep -c '^>' "$designdir/HA_r2.fasta" || echo 0)
    echo "[INFO] Round 2: ${n2} cluster representatives"
    
    if (( n2 == 0 )); then
        echo "[WARN] Round 2 produced no clusters - using round 1"
        cp "$designdir/HA_r1_cons.fasta" "$designdir/${PREFIX}_${year}_COBRA.fasta"
        sed -i "s/>HA_R1_CONS/>${PREFIX}_${year}_COBRA/" "$designdir/${PREFIX}_${year}_COBRA.fasta"
        return 0
    fi
    
    $MAFFT --auto "$designdir/HA_r2.fasta" \
        > "$designdir/HA_r2_aln.fasta" \
        2> "$logdir/mafft_r2.log" || {
        echo "[WARN] MAFFT round 2 failed - using round 1"
        cp "$designdir/HA_r1_cons.fasta" "$designdir/${PREFIX}_${year}_COBRA.fasta"
        sed -i "s/>HA_R1_CONS/>${PREFIX}_${year}_COBRA/" "$designdir/${PREFIX}_${year}_COBRA.fasta"
        return 0
    }
    
    # Final consensus
    $PY <<EOF
from collections import Counter
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

aln = AlignIO.read("${designdir}/HA_r2_aln.fasta", "fasta")
cons = []
for i in range(aln.get_alignment_length()):
    col = [str(r.seq)[i] for r in aln]
    c = Counter(col)
    if "-" in c:
        c["-"] = int(c["-"] * 0.2)
    best = max(c, key=c.get)
    cons.append("X" if best == "-" else best)

record = SeqRecord(Seq("".join(cons)), id="${PREFIX}_${year}_COBRA", description="")
SeqIO.write([record], "${designdir}/${PREFIX}_${year}_COBRA.fasta", "fasta")
print("[OK] Final COBRA consensus written")
EOF
}

# =========================================================
#  Process each year
# =========================================================
for year_dir in ${OUT_DIR}/*; do
    year=$(basename "$year_dir")
    
    # Priority: 1) INPUT_DIR if provided, 2) year_dir unaligned, 3) year_dir aligned
    if [[ -n "$INPUT_DIR" ]] && [[ -s "${INPUT_DIR}/${PREFIX}_${year}.fasta" ]]; then
        input_fasta="${INPUT_DIR}/${PREFIX}_${year}.fasta"
    elif [[ -s "${year_dir}/${PREFIX}_${year}.fasta" ]]; then
        input_fasta="${year_dir}/${PREFIX}_${year}.fasta"
        $PY <<EOF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

records = []
for rec in SeqIO.parse("${year_dir}/${PREFIX}_${year}_aligned.fasta", "fasta"):
    ungapped = str(rec.seq).replace("-", "")
    records.append(SeqRecord(Seq(ungapped), id=rec.id, description=rec.description))

SeqIO.write(records, "${input_fasta}", "fasta")
print(f"[OK] Ungapped {len(records)} sequences")
EOF
    else
        input_fasta=""
    fi
    
    designdir="${year_dir}/designs"
    logdir="${year_dir}/logs"
    
    mkdir -p "$designdir" "$logdir"
    
    echo "---------------------------------------------"
    echo "[CHECK] Year: ${year}"
    echo "[CHECK] Input: ${input_fasta}"
    
    if [[ ! -s "$input_fasta" ]]; then
        echo "[SKIP] No fasta file for ${year}"
        continue
    fi
    
    seqs=$(grep -c "^>" "$input_fasta" || echo 0)
    if (( seqs < 2 )); then
        echo "[SKIP] Only ${seqs} sequences for ${year}"
        continue
    fi
    
    echo "[INFO] Processing ${year} with ${seqs} sequences"
    
    if safe_cobra "$input_fasta" "$year" "$designdir" "$logdir"; then
        echo "[DONE] ${year} COBRA complete ✅"
    else
        echo "[ERROR] ${year} COBRA failed"
    fi
done

echo "=== All COBRA runs complete 🚀 ==="
