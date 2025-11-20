# Results Organization

This directory contains two parallel analyses:

## Clustered Analysis (Primary)

**Location:** `results/clustered/`

**Data:** Sequences clustered at 99% identity using CD-HIT
- H1N1: 535 representative sequences
- H3N2: 806 representative sequences  
- VicB: 89 representative sequences
- **Total: 1,430 sequences**

**Purpose:** 
- Main analysis pipeline
- Publication-ready figures
- Faster computational time
- Reduces sampling bias

---

## Unclustered Analysis (Validation)

**Location:** `results/unclustered/`

**Data:** All filtered, high-quality sequences (no clustering)
- H1N1: 12,835 sequences
- H3N2: 15,595 sequences
- VicB: 4,883 sequences  
- **Total: 33,313 sequences**

**Purpose:**
- Validate clustering didn't lose important diversity
- Comprehensive diversity metrics
- Supplementary analysis
- Quality assurance

---

## Comparison

Both analyses follow the same pipeline:
1. Multiple sequence alignment (MAFFT)
2. Phylogenetic tree building (IQ-TREE)
3. Vaccine design (Consensus, Medoid, ASR, COBRA)
4. Distance analysis
5. Statistical visualization

The clustered analysis is **~20-30x faster** while maintaining biological insights.
