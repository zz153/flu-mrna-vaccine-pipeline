"""
Script: 04_find_conserved_stretches.py
Description:
    - Slides a window (default: 15 aa) across an MSA to identify and rank all possible conserved peptide stretches.
    - Calculates mean conservation for each window (fraction of most frequent amino acid per column, averaged).
    - Outputs a TSV file with all windows (start, end, mean conservation).
    - Outputs a FASTA file with the consensus sequence for each window.
Usage:
    python 05_find_conserved_stretches.py \
        --alignment ../data/processed/combined_HA_aligned.fa \
        --window 15 \
        --tsv ../results/tables/all_conserved_15aa_windows.tsv \
        --fasta ../results/peptides/conserved_15aa_consensus.fasta
"""

import argparse
from Bio import AlignIO, SeqIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

def compute_mean_conservation(aln, window_size):
    aln_length = aln.get_alignment_length()
    scores = []
    for i in range(aln_length - window_size + 1):
        window_conservation = []
        for j in range(i, i+window_size):
            column = aln[:, j]
            counts = {}
            for aa in column:
                if aa not in ['-', 'X', 'B', 'J', 'Z']:
                    counts[aa] = counts.get(aa, 0) + 1
            major_freq = max(counts.values()) / sum(counts.values()) if counts else 0
            window_conservation.append(major_freq)
        avg_cons = sum(window_conservation) / window_size
        scores.append((avg_cons, i+1, i+window_size))  # 1-based indexing
    return sorted(scores, reverse=True)

def main():
    parser = argparse.ArgumentParser(description="Identify and rank conserved windows in a protein MSA.")
    parser.add_argument("--alignment", "-a", default="../data/processed/combined_HA_aligned.fa", help="Input alignment FASTA")
    parser.add_argument("--window", "-w", type=int, default=15, help="Window size")
    parser.add_argument("--tsv", default="../results/tables/all_conserved_15aa_windows.tsv", help="Output TSV file")
    parser.add_argument("--fasta", default="../results/peptides/conserved_15aa_consensus.fasta", help="Output FASTA file")
    args = parser.parse_args()

    # Ensure output directories exist
    os.makedirs(os.path.dirname(args.tsv), exist_ok=True)
    os.makedirs(os.path.dirname(args.fasta), exist_ok=True)

    print(f"Loading alignment from {args.alignment}")
    aln = AlignIO.read(args.alignment, "fasta")
    print("Computing mean conservation for all windows...")
    scores = compute_mean_conservation(aln, args.window)

    # Write TSV of all windows
    print(f"Writing window stats to {args.tsv}")
    with open(args.tsv, "w") as f:
        f.write("start\tend\tmean_conservation\n")
        for avg_cons, start, end in scores:
            f.write(f"{start}\t{end}\t{avg_cons:.3f}\n")

    # Write consensus for each window
    print(f"Writing consensus sequences to {args.fasta}")
    cons_records = []
    for avg_cons, start, end in scores:
        sub_aln = aln[:, start-1:end]
        summary = AlignInfo.SummaryInfo(sub_aln)
        consensus = summary.dumb_consensus(threshold=0.7, ambiguous='X')
        rec = SeqRecord(Seq(str(consensus)), id=f"Pos{start}-{end}_cons{avg_cons:.3f}", description="")
        cons_records.append(rec)
    SeqIO.write(cons_records, args.fasta, "fasta")
    print("Done.")

if __name__ == "__main__":
    main()
