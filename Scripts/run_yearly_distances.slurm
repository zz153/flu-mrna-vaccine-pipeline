#!/bin/bash -e
#SBATCH --job-name=yearly_distances
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --output=logs/yearly_distances.%j.out
#SBATCH --error=logs/yearly_distances.%j.err

# =========================================================
#  Calculate pairwise distances between designs & strains
#  Outputs: p-distance and LG model distance summaries
#  Usage:
#     sbatch run_yearly_distances.slurm H1N1 h1n1
# =========================================================

set -euo pipefail

# ---- Environment ----
ENV=/home/ranzo85p/conda_envs/bio
PY=$ENV/bin/python
MAFFT=$ENV/bin/mafft
IQT=$ENV/bin/iqtree2

PREFIX=${1:-H1N1}
OUT_DIR=${2:-h1n1}
CPUS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p logs

echo "=== Tool versions ==="
$PY --version
$PY -c "import Bio; print('Biopython', Bio.__version__)"
$IQT --version | head -1
echo "=== Starting distance calculations for ${OUT_DIR} ==="

# =========================================================
#  Process each year
# =========================================================
for year_dir in ${OUT_DIR}/*; do
    year=$(basename "$year_dir")
    aln="${year_dir}/${PREFIX}_${year}_aligned.fasta"
    designdir="${year_dir}/designs"
    distdir="${year_dir}/dist"
    logdir="${year_dir}/logs"
    
    mkdir -p "$distdir" "$logdir"
    
    echo "---------------------------------------------"
    echo "[CHECK] Year: ${year}"
    
    if [[ ! -s "$aln" ]]; then
        echo "[SKIP] No alignment for ${year}"
        continue
    fi
    
    # Check if we have design sequences
    design_count=$(ls "$designdir"/${PREFIX}_${year}_*.fasta 2>/dev/null | wc -l)
    if [[ $design_count -eq 0 ]]; then
        echo "[SKIP] No design sequences for ${year}"
        continue
    fi
    
    echo "[INFO] Found ${design_count} design sequences for ${year}"
    
    # =========================================================
    # Step 1: Combine designs + circulating strains
    # =========================================================
    combined="${year_dir}/${PREFIX}_${year}_designs_plus_strains.fasta"
    
    echo "[INFO] Combining designs and strains..."
    
    # First, create ungapped versions and combine
    $PY <<EOF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

all_records = []

# Add all design sequences (ungapped)
design_types = ['CONS', 'ASR', 'MEDOID', 'CENTROID', 'COBRA']
for design in design_types:
    design_file = f"${designdir}/${PREFIX}_${year}_{design}.fasta"
    try:
        for rec in SeqIO.parse(design_file, "fasta"):
            ungapped = str(rec.seq).replace("-", "")
            all_records.append(SeqRecord(Seq(ungapped), id=rec.id, description=rec.description))
    except FileNotFoundError:
        pass

# Add original circulating strains (ungapped)
for rec in SeqIO.parse("${aln}", "fasta"):
    ungapped = str(rec.seq).replace("-", "")
    all_records.append(SeqRecord(Seq(ungapped), id=rec.id, description=rec.description))

# Write all to combined file
SeqIO.write(all_records, "${combined}", "fasta")
print(f"[INFO] Combined {len(all_records)} sequences")
EOF
    
    combined_count=$(grep -c '^>' "$combined" || echo 0)
    echo "[INFO] Combined file has ${combined_count} sequences"
    
    # =========================================================
    # Step 2: Align combined sequences
    # =========================================================
    combined_aln="${year_dir}/${PREFIX}_${year}_combined_aligned.fasta"
    
    echo "[INFO] Aligning combined sequences..."
    $MAFFT --auto --thread "$CPUS" "$combined" > "$combined_aln" 2> "$logdir/mafft_combined.log"
    
    # =========================================================
    # Step 3: Calculate p-distances (Python)
    # =========================================================
    echo "[INFO] Calculating p-distances..."
    $PY <<EOF
import csv
import statistics as st
from Bio import AlignIO

aln = AlignIO.read("${combined_aln}", "fasta")
seq_ids = [rec.id for rec in aln]
seqs = [str(rec.seq) for rec in aln]

# Identify designs vs strains
design_keywords = ['CONS', 'ASR', 'MEDOID', 'CENTROID', 'COBRA']
is_design = lambda s: any(k in s for k in design_keywords)

design_indices = [i for i, sid in enumerate(seq_ids) if is_design(sid)]
strain_indices = [i for i, sid in enumerate(seq_ids) if not is_design(sid)]

print(f"[INFO] Designs: {len(design_indices)}, Strains: {len(strain_indices)}")

# p-distance: proportion of differences (ignoring gaps)
def p_distance(seq_a, seq_b):
    comparisons = 0
    mismatches = 0
    for a, b in zip(seq_a, seq_b):
        if a != '-' and b != '-':
            comparisons += 1
            if a != b:
                mismatches += 1
    return (mismatches / comparisons) if comparisons > 0 else float('nan')

# Calculate design → strain distances
results = []
for d_idx in design_indices:
    design_name = seq_ids[d_idx]
    design_seq = seqs[d_idx]
    
    distances = []
    for s_idx in strain_indices:
        strain_seq = seqs[s_idx]
        dist = p_distance(design_seq, strain_seq)
        if dist == dist:  # not NaN
            distances.append(dist)
    
    if distances:
        results.append({
            'design': design_name,
            'n_strains': len(distances),
            'mean_pdist': st.mean(distances),
            'median_pdist': st.median(distances),
            'min_pdist': min(distances),
            'max_pdist': max(distances),
            'stdev_pdist': st.stdev(distances) if len(distances) > 1 else 0.0
        })

# Write summary
if results:
    with open("${distdir}/p_distance_summary.tsv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(results)
    print(f"[OK] p-distance summary written to ${distdir}/p_distance_summary.tsv")
else:
    print("[WARN] No valid p-distances calculated")
EOF
    
    # =========================================================
    # Step 4: Calculate evolutionary model distances (IQ-TREE)
    # =========================================================
    echo "[INFO] Calculating evolutionary model distances with IQ-TREE..."
    
    # Use LG+G for proteins, calculate ML distances only (no tree needed for distance matrix)
    $IQT -s "$combined_aln" -m LG+G -nt "$CPUS" \
        -pre "$logdir/evoldist" --redo \
        > "$logdir/iqtree_evoldist.log" 2>&1 || {
        echo "[WARN] IQ-TREE distance calculation failed for ${year}"
        echo "[DONE] ${year} distance calculations complete (p-distance only) ✅"
        continue
    }
    
    # =========================================================
    # Step 5: Summarize evolutionary distances
    # =========================================================
    echo "[INFO] Summarizing evolutionary model distances..."
    $PY <<EOF
import csv
import numpy as np
import os

mldist_file = "${logdir}/evoldist.mldist"
if not os.path.exists(mldist_file):
    print(f"[WARN] ML distance file not found: {mldist_file}")
    exit(0)

# Parse IQ-TREE mldist format (lower triangular matrix)
with open(mldist_file) as f:
    lines = [l.strip().split() for l in f if l.strip()]

n = int(lines[0][0])
labels = [lines[i+1][0] for i in range(n)]

# Build full distance matrix
dist_matrix = np.zeros((n, n))
for i in range(n):
    row_data = lines[i+1][1:]  # skip the label
    for j, val in enumerate(row_data):
        if val != '-':
            dist_matrix[i][j] = float(val)
            dist_matrix[j][i] = float(val)  # symmetric

# Identify designs vs strains
design_keywords = ['CONS', 'ASR', 'MEDOID', 'CENTROID', 'COBRA']
is_design = lambda s: any(k in s for k in design_keywords)

design_indices = [i for i, lbl in enumerate(labels) if is_design(lbl)]
strain_indices = [i for i, lbl in enumerate(labels) if not is_design(lbl)]

print(f"[INFO] Designs: {len(design_indices)}, Strains: {len(strain_indices)}")

# Summarize design → strain distances
results = []
for d_idx in design_indices:
    design_name = labels[d_idx]
    
    distances = []
    for s_idx in strain_indices:
        dist = dist_matrix[d_idx][s_idx]
        if not np.isnan(dist) and dist > 0:
            distances.append(dist)
    
    if distances:
        results.append({
            'design': design_name,
            'n_strains': len(distances),
            'mean_evol_dist': float(np.mean(distances)),
            'median_evol_dist': float(np.median(distances)),
            'min_evol_dist': float(np.min(distances)),
            'max_evol_dist': float(np.max(distances)),
            'stdev_evol_dist': float(np.std(distances)) if len(distances) > 1 else 0.0
        })

# Write summary
if results:
    with open("${distdir}/evolutionary_distance_summary.tsv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(results)
    print(f"[OK] Evolutionary distance summary written to ${distdir}/evolutionary_distance_summary.tsv")
else:
    print("[WARN] No valid evolutionary distances calculated")
EOF
    
    echo "[DONE] ${year} distance calculations complete ✅"
done

echo "=== All distance calculations complete 🚀 ==="


