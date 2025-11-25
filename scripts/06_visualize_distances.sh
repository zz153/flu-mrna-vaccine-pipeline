#!/bin/bash
# Script 06: Visualize Vaccine Design Performance
#
# Creates visualizations comparing vaccine designs

set -e

LINEAGE=$1
SOURCE_TYPE=${2:-clustered}

if [ -z "$LINEAGE" ]; then
    echo "Usage: bash scripts/06_visualize_distances.sh [LINEAGE] [SOURCE_TYPE]"
    echo "Example: bash scripts/06_visualize_distances.sh H1N1 clustered"
    exit 1
fi

echo "=========================================="
echo "Creating Visualizations for ${LINEAGE}"
echo "Source: ${SOURCE_TYPE}"
echo "=========================================="

INPUT_DIR="results/per_year_${SOURCE_TYPE}/${LINEAGE}"
OUTPUT_DIR="${INPUT_DIR}/figures"
mkdir -p ${OUTPUT_DIR}

# Check if distance files exist
if [ ! -d "${INPUT_DIR}/distances" ]; then
    echo "❌ Error: Distance directory not found!"
    exit 1
fi

N_FILES=$(find ${INPUT_DIR}/distances -name "*_distances.csv" | wc -l)
echo "Found ${N_FILES} distance files"

if [ "$N_FILES" -eq 0 ]; then
    echo "❌ Error: No distance files found!"
    exit 1
fi

echo ""
echo "Creating visualizations..."

# Generate plots using Python - pass variables as arguments
python3 - "$LINEAGE" "$SOURCE_TYPE" "$INPUT_DIR" "$OUTPUT_DIR" << 'PYEOF'
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import sys

# Get arguments
lineage = sys.argv[1]
source_type = sys.argv[2]
input_dir = Path(sys.argv[3])
output_dir = Path(sys.argv[4])
distances_dir = input_dir / "distances"

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10

# Read summary files - NOW INCLUDING COBRA!
designs = ['consensus', 'medoid', 'ancestral', 'cobra']
summary_data = {}

for design in designs:
    summary_file = distances_dir / f"distance_summary_{design}.csv"
    if summary_file.exists():
        df = pd.read_csv(summary_file)
        summary_data[design] = df
        print(f"  Loaded {design}: {len(df)} years")
    else:
        print(f"  ⚠️  Missing {design} summary")

if not summary_data:
    print("❌ No summary files found!")
    sys.exit(1)

# ============================================
# Plot 1: Year-by-Year Mean Distance Comparison
# ============================================
print("\n📊 Creating Plot 1: Year-by-year comparison...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))

# Define colors for each design
colors = {
    'consensus': '#3498db',
    'medoid': '#2ecc71', 
    'ancestral': '#e74c3c',
    'cobra': '#f39c12'
}

# Plot mean p-distance
for design, df in summary_data.items():
    ax1.plot(df['Year'], df['Mean_ML_Distance'], 
             marker='o', linewidth=2.5, markersize=7, 
             label=design.capitalize(), alpha=0.85,
             color=colors.get(design, '#95a5a6'))
    ax1.fill_between(df['Year'], 
                      df['Mean_ML_Distance'] - df['SD_ML_Distance'],
                      df['Mean_ML_Distance'] + df['SD_ML_Distance'],
                      alpha=0.15, color=colors.get(design, '#95a5a6'))

ax1.set_xlabel('Year', fontsize=13, fontweight='bold')
ax1.set_ylabel('Mean P-Distance', fontsize=13, fontweight='bold')
ax1.set_title(f'{lineage} Vaccine Design Performance (P-Distance)', 
              fontsize=15, fontweight='bold')
ax1.legend(loc='best', fontsize=12, framealpha=0.9)
ax1.grid(True, alpha=0.3)

# Plot ML distance (for ALL designs now!)
for design, df in summary_data.items():
    valid_ml = df[~df['Mean_ML_Distance'].isna()]
    if len(valid_ml) > 0:
        ax2.plot(valid_ml['Year'], valid_ml['Mean_ML_Distance'],
                marker='s', linewidth=2.5, markersize=7,
                label=f'{design.capitalize()} (ML)', alpha=0.85,
                color=colors.get(design, '#95a5a6'))
        ax2.fill_between(valid_ml['Year'],
                         valid_ml['Mean_ML_Distance'] - valid_ml['SD_ML_Distance'],
                         valid_ml['Mean_ML_Distance'] + valid_ml['SD_ML_Distance'],
                         alpha=0.15, color=colors.get(design, '#95a5a6'))

ax2.set_xlabel('Year', fontsize=13, fontweight='bold')
ax2.set_ylabel('Mean ML Distance (LG+G4)', fontsize=13, fontweight='bold')
ax2.set_title(f'{lineage} ML Distance Performance (All Designs)',
              fontsize=15, fontweight='bold')
ax2.legend(loc='best', fontsize=11, framealpha=0.9)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
output_file = output_dir / f"{lineage}_yearly_comparison.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()
print(f"  ✅ Saved: {output_file}")

# ============================================
# # # Plot 2: Heatmap of Mean Distances (CLEANER!)
# # # ============================================
# # print("\n📊 Creating Plot 2: Distance heatmap...")
# # 
# # # Prepare data for heatmap
# # heatmap_data = []
# # years = sorted(set().union(*[set(df['Year']) for df in summary_data.values()]))
# # 
# # for year in years:
# #     row = {'Year': year}
# #     for design, df in summary_data.items():
# #         year_data = df[df['Year'] == year]
# #         if len(year_data) > 0:
# #             row[design.capitalize()] = year_data['Mean_ML_Distance'].values[0]
# #         else:
# #             row[design.capitalize()] = np.nan
# #     heatmap_data.append(row)
# # 
# # heatmap_df = pd.DataFrame(heatmap_data).set_index('Year')
# # 
# # # Create cleaner heatmap WITHOUT numbers
# # fig, ax = plt.subplots(figsize=(14, 8))
# # sns.heatmap(heatmap_df.T, annot=False, cmap='RdYlGn_r',
# #             cbar_kws={'label': 'Mean P-Distance'},
# #             linewidths=1, linecolor='white', ax=ax,
# #             vmin=0, vmax=0.5)  # Fixed scale for consistency
# # 
# # ax.set_title(f'{lineage} Mean P-Distance by Design and Year',
# #              fontsize=16, fontweight='bold', pad=20)
# # ax.set_xlabel('Year', fontsize=13, fontweight='bold')
# # ax.set_ylabel('Design Type', fontsize=13, fontweight='bold')
# # ax.tick_params(axis='both', labelsize=11)
# # 
# # plt.tight_layout()
# # output_file = output_dir / f"{lineage}_distance_heatmap.png"
# # plt.savefig(output_file, dpi=300, bbox_inches='tight')
# # plt.close()
# # print(f"  ✅ Saved: {output_file}")
# # 
# ============================================
# Plot 3: Overall Performance Comparison
# ============================================
print("\n📊 Creating Plot 3: Overall comparison...")

fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Plot A: Mean distance comparison
ax = axes[0, 0]
means = {design: df['Mean_ML_Distance'].mean() for design, df in summary_data.items()}
design_colors = [colors.get(d, '#95a5a6') for d in means.keys()]
bars = ax.bar(means.keys(), means.values(), color=design_colors, alpha=0.8, edgecolor='black', linewidth=1.5)
ax.set_ylabel('Overall Mean P-Distance', fontsize=12, fontweight='bold')
ax.set_title('Average Performance Across All Years', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')

for bar in bars:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height,
            f'{height:.4f}', ha='center', va='bottom', fontweight='bold', fontsize=10)

# Plot B: Median distance comparison
ax = axes[0, 1]
medians = {design: df['Median_ML_Distance'].median() for design, df in summary_data.items()}
bars = ax.bar(medians.keys(), medians.values(), color=design_colors, alpha=0.8, edgecolor='black', linewidth=1.5)
ax.set_ylabel('Overall Median P-Distance', fontsize=12, fontweight='bold')
ax.set_title('Median Performance Across All Years', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')

for bar in bars:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height,
            f'{height:.4f}', ha='center', va='bottom', fontweight='bold', fontsize=10)

# Plot C: Number of strains per year
ax = axes[1, 0]
design = list(summary_data.keys())[0]
df = summary_data[design]
ax.bar(df['Year'], df['N_Strains'], color='#95a5a6', alpha=0.8, edgecolor='black', linewidth=0.5)
ax.set_xlabel('Year', fontsize=12, fontweight='bold')
ax.set_ylabel('Number of Strains', fontsize=12, fontweight='bold')
ax.set_title('Strains Analyzed Per Year', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')

# Plot D: Variability comparison
ax = axes[1, 1]
variability = {design: df['SD_ML_Distance'].mean() for design, df in summary_data.items()}
bars = ax.bar(variability.keys(), variability.values(), color=design_colors, alpha=0.8, edgecolor='black', linewidth=1.5)
ax.set_ylabel('Average Standard Deviation', fontsize=12, fontweight='bold')
ax.set_title('Distance Variability Across Years', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')

for bar in bars:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height,
            f'{height:.4f}', ha='center', va='bottom', fontweight='bold', fontsize=10)

plt.tight_layout()
output_file = output_dir / f"{lineage}_overall_comparison.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()
print(f"  ✅ Saved: {output_file}")

# ============================================
# Generate Summary Statistics
# ============================================
print("\n📊 Generating summary statistics...")

summary_stats = []
for design, df in summary_data.items():
    stats = {
        'Design': design.capitalize(),
        'Mean_Distance': df['Mean_ML_Distance'].mean(),
        'Median_Distance': df['Median_ML_Distance'].median(),
        'Min_Distance': df['Mean_ML_Distance'].min(),
        'Max_Distance': df['Mean_ML_Distance'].max(),
        'Avg_SD': df['SD_ML_Distance'].mean(),
        'Years_Analyzed': len(df)
    }
    summary_stats.append(stats)

summary_df = pd.DataFrame(summary_stats)
summary_file = output_dir / f"{lineage}_summary_statistics.csv"
summary_df.to_csv(summary_file, index=False)
print(f"  ✅ Saved: {summary_file}")

print("\n" + "="*50)
print("Summary Statistics:")
print("="*50)
print(summary_df.to_string(index=False))
print("="*50)
PYEOF

echo ""
echo "=========================================="
echo "✅ Visualizations complete for ${LINEAGE}!"
echo "=========================================="
echo ""
echo "Output directory: ${OUTPUT_DIR}"
ls -lh ${OUTPUT_DIR}
