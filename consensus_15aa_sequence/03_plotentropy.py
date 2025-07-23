#!/usr/bin/env python3
"""
Script: 03_plot_entropy.py
Description: Plots Shannon entropy by alignment position for the HA MSA.
Usage: python 06_plot_entropy.py --input ../results/tables/HA_shannon_entropy.tsv --output ../results/figures/HA_entropy_plot.png
"""

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Plot Shannon entropy profile from TSV file.")
    parser.add_argument("--input", "-i", default="../results/tables/HA_shannon_entropy.tsv", help="Input entropy TSV file")
    parser.add_argument("--output", "-o", default="../results/figures/HA_entropy_plot.png", help="Output plot file (PNG)")
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    print(f"Loading entropy table from {args.input}")
    df = pd.read_csv(args.input, sep="\t")
    
    plt.figure(figsize=(12, 4))
    plt.plot(df["pos"], df["entropy"], lw=1, color="teal")
    plt.xlabel("Alignment Position (Amino Acid)")
    plt.ylabel("Shannon Entropy")
    plt.title("Conservation Landscape of HA (Shannon Entropy)")
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Plot saved to {args.output}")
    plt.show()

if __name__ == "__main__":
    main()
