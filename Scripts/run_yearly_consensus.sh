#!/bin/bash
#SBATCH --job-name=yearly_consensus
#SBATCH --cpus-per-task=4
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --output=logs/yearly_consensus.%j.out
#SBATCH --error=logs/yearly_consensus.%j.err

# ===== USER CONFIG =====
ENV_PATH="/home/ranzo85p/conda_envs/bio"
PREFIX=${1:-"H1N1"}   # default to H1N1 if not provided
OUT_DIR=${2:-"h1n1"}  # default to h1n1 if not provided

# ===== SETUP =====
PY="${ENV_PATH}/bin/python"
mkdir -p logs

echo "=== Starting consensus generation for all years in ${OUT_DIR} ==="

for year_dir in ${OUT_DIR}/*; do
    year=$(basename "$year_dir")
    aln="${year_dir}/${PREFIX}_${year}_aligned.fasta"
    outdir="${year_dir}/designs"
    outfile="${outdir}/${PREFIX}_${year}_CONS.fasta"

    if [[ ! -f "$aln" ]]; then
        echo "[SKIP] No alignment found for ${year}"
        continue
    fi

    mkdir -p "$outdir"

    echo "[INFO] Generating consensus for ${year} ..."
    $PY <<EOF > "$outfile"
from collections import Counter
from Bio import AlignIO
aln = AlignIO.read("${aln}", "fasta")
cons = []
for i in range(aln.get_alignment_length()):
    col = [str(r.seq)[i] for r in aln]
    c = Counter(col)
    if "-" in c:
        c["-"] = int(c["-"] * 0.2)
    best = max(c, key=c.get)
    cons.append("X" if best == "-" else best)
print(">${PREFIX}_${year}_CONS")
print("".join(cons))
EOF

    echo "[DONE] Consensus saved to $outfile"
done

echo "=== Consensus generation complete 🚀 ==="
