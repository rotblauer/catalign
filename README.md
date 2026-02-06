# Catalign

**Where sequences find each other like cats find the warmest spot.**

Catalign is a bioinformatics sequence alignment tool that mimics natural molecular forces for aligning DNA sequences. Instead of arbitrary scoring heuristics, Catalign models alignment as an energy minimisation problem — long-range attraction finds approximate matching regions, short-range forces refine base-level alignment, and the final result is the most energetically stable configuration.

---

## Quick Demo with Real Data

The fastest way to see Catalign in action with real genomic data:

```bash
# Install with visualization support
pip install catalign[viz]

# Run the mitochondrial genome demo
python examples/mito_demo.py
```

This downloads human mitochondrial sequences from NCBI and generates:
- Interactive dot plots showing sequence similarity
- Energy and identity heatmaps
- CIGAR string visualizations
- Quality assessment reports

See `examples/README.md` for more demo options.

## Why Catalign?

Current sequence aligners (minimap2, BWA-MEM, LAST, …) are engineering marvels, but their scoring models diverge from the biology they serve:

| Problem with current aligners | Catalign's natural approach |
|---|---|
| **Seed-and-extend** is rigid — seeds are either exact or ignored | **Long-range attraction** via minimizer sketches provides a smooth, distance-weighted signal that degrades gracefully |
| **Gap penalties** (affine, convex) are mathematical conveniences, not biology | **Energy-based gaps** model strand separation: opening a gap costs real energy, extending it costs less — just like pulling apart a DNA duplex |
| **Poor handling of structural variants** — large indels fall outside band, SVs break aligners | **Multi-scale chaining** naturally handles large rearrangements: anchor chains can span structural variants |
| **Repetitive regions** cause spurious mappings and MQ0 scores | **Energy minimisation** across the full alignment landscape finds the globally optimal placement, not just the first adequate seed chain |
| **No multi-scale quality assessment** — MAPQ is a single number that hides local problems | **Four-level quality evaluation** (base → block → region → genome) exposes exactly where and why an alignment is uncertain |

## Architecture

```
┌──────────────────────────────────────────────────┐
│                  Catalign Pipeline                │
│                                                  │
│  ┌────────────┐   Long-range     ┌────────────┐ │
│  │  Minimizer  │  attraction     │   Anchor    │ │
│  │  Sketching  │ ──────────────► │  Chaining   │ │
│  └────────────┘   (k-mer match)  └─────┬──────┘ │
│                                        │         │
│                                        ▼         │
│                                  ┌────────────┐  │
│                Short-range       │   Banded    │  │
│                forces            │  DP Align   │  │
│                                  └─────┬──────┘  │
│                                        │         │
│                                        ▼         │
│                                  ┌────────────┐  │
│                Energy            │  Multi-Scale│  │
│                minimisation      │  Quality    │  │
│                                  └────────────┘  │
└──────────────────────────────────────────────────┘
```

1. **Sketching** (`catalign.sketch`) — Extract minimizer profiles from query and target. These compressed representations enable efficient long-range comparison.
2. **Anchor finding & chaining** (`catalign.chain`) — Shared minimizers become anchors; dynamic programming chains compatible anchors into co-linear groups.
3. **Banded alignment** (`catalign.align`) — Between anchors, a banded Smith-Waterman-style DP fills in base-level alignment using the energy model.
4. **Quality evaluation** (`catalign.quality`) — The finished alignment is assessed at four scales: base, block, region, and genome.

## Installation

```bash
pip install catalign
```

For development with all extras:

```bash
git clone https://github.com/your-org/catalign.git
cd catalign
pip install -e ".[dev]"
```

Optional dependency groups:
- `test` — pytest, coverage
- `viz` — plotly, streamlit (for interactive visualization)
- `bench` — pytest-benchmark, memory-profiler
- `all` — all of the above

## Quick Start

### Command Line

```bash
# Align two FASTA files
catalign align query.fa target.fa

# PAF-style output
catalign align query.fa target.fa --output paf

# Custom energy parameters
catalign align query.fa target.fa --match-energy -2.5 --gap-open 6.0

# Quality evaluation
catalign quality query.fa target.fa
```

### Python API

```python
from catalign import CatalignAligner, EnergyModel, evaluate_quality

# Configure energy model
em = EnergyModel(match_energy=-2.0, mismatch_energy=3.0,
                 gap_open_energy=5.0, gap_extend_energy=1.0)

# Align
aligner = CatalignAligner(energy_model=em, k=15, w=50)
aln = aligner.align(query_seq, target_seq)

print(f"CIGAR: {aln.cigar}")
print(f"Energy: {aln.energy_score}")

# Multi-scale quality
qual = evaluate_quality(aln, query_seq, target_seq)
print(f"Identity: {qual.overall_identity:.2%}")
print(f"Quality score: {qual.quality_score:.1f}/100")
```

## Multi-Scale Quality Evaluation

Catalign evaluates alignment quality at four levels:

### Base Level
Per-position assessment: is this a match, mismatch, insertion, or deletion? Each base carries an energy value reflecting confidence.

### Block Level
Contiguous aligned segments are grouped into blocks. Each block reports identity percentage, length, and cumulative energy — revealing whether a region is solidly matched or weakly aligned.

### Region Level
Blocks are aggregated into larger regions, reporting structural concordance, query/target coverage, and the number of aligned blocks. This is where structural variants become visible.

### Genome Level
A `QualityReport` aggregates all alignments for a genome-wide view: mean identity, mean quality score, total energy, and overall coverage.

## Test Cases

### Running tests

```bash
pip install -e ".[test]"
pytest tests/ -v
```

### Generating Synthetic Test Data

Catalign includes a comprehensive synthetic test data generator for validation:

```bash
python scripts/generate_test_data.py
```

This creates test sequences in `tests/resources/` with:
- **Identical sequences** — baseline for 100% identity alignment
- **SNP mutations** — single and multiple nucleotide polymorphisms
- **Indels** — insertions and deletions of various sizes
- **Structural variants** — inversions, tandem duplications, large deletions
- **Repetitive sequences** — tandem repeats and interspersed elements
- **Benchmark sequences** — larger sequences for performance testing

Each test case includes a ground truth JSON file with expected alignment properties.

### Running Benchmarks

```bash
# Run benchmark suite against ground truth
catalign benchmark --resources-dir tests/resources

# Generate metrics report
catalign metrics query.fa target.fa --json
```

The test suite includes synthetic sequences covering:
- Identical sequence alignment
- Single-base mismatches
- Small insertions and deletions
- Repetitive sequences
- Longer sequences with shared cores

### Real-world test data

For production-scale testing, we recommend:

- **T2T CHM13 v2.0** — the first complete human genome assembly from the Telomere-to-Telomere Consortium. Download from [NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCA_009914755.4/).
- **HPRC assemblies** — haplotype-resolved assemblies (hap1 vs hap2) from the Human Pangenome Reference Consortium, ideal for testing structural variant handling and haplotype-aware alignment. Available at [https://humanpangenome.org](https://humanpangenome.org).

Example with real data:

```bash
# Align hap1 vs hap2 for a single chromosome
catalign align hap1_chr1.fa hap2_chr1.fa --kmer-size 19 --window-size 100 --output paf
```

## Computational Efficiency

- **Minimizer sketching** reduces the sequence representation by ~50× compared to all k-mers, enabling fast long-range comparison.
- **Anchor chaining** uses a look-back window in the DP to keep chaining O(n·w) instead of O(n²).
- **Banded alignment** constrains the DP matrix to a diagonal band, reducing base-level alignment from O(nm) to O(n·bandwidth).
- **NumPy arrays** are used for base encoding and energy matrices for vectorised computation.

For whole-genome alignment, future versions will add:
- Tiled/blocked processing for memory efficiency
- Multiprocessing for independent chromosome arms
- Optional GPU acceleration via CuPy

## Roadmap

- [x] **Visualisation** — quality heat-maps, dot-plots, energy landscape plots, and interactive dashboard
- [ ] **GPU acceleration** — CuPy backend for banded DP on large sequences
- [ ] **Population-scale alignment** — align many samples against a pangenome graph
- [ ] **RNA-seq support** — splice-aware energy model
- [ ] **Structural variant calling** — leverage anchor chain breaks as SV evidence

## Interactive Visualization Dashboard

Catalign includes a comprehensive visualization suite for alignment analysis. Launch the interactive dashboard:

```bash
catalign viz
# Or specify port
catalign viz --port 8080
```

### Python Visualization API

```python
from catalign import CatalignAligner, evaluate_quality
from catalign.viz import (
    create_dotplot,
    create_energy_heatmap,
    create_identity_heatmap,
    visualize_cigar,
    create_alignment_view,
)

# Align sequences
aligner = CatalignAligner()
aln = aligner.align(query_seq, target_seq)

# Dot plot
dp = create_dotplot(query_seq, target_seq, k=11)
fig = dp.to_figure()
fig.write_html("dotplot.html")

# Energy landscape
energy_fig = create_energy_heatmap(aln, query_seq, target_seq)

# CIGAR visualization
cigar_viz = visualize_cigar(aln.cigar)
cigar_viz.to_figure().show()

# Text alignment view
view = create_alignment_view(aln, query_seq, target_seq)
print(view.to_text())
```

### Visualization Types

| Visualization | Description |
|--------------|-------------|
| **Dot Plot** | K-mer match visualization showing sequence similarity patterns |
| **Energy Heatmap** | Sliding window energy across the alignment |
| **Identity Heatmap** | Local identity percentage along the alignment |
| **CIGAR View** | Colored block representation of alignment operations |
| **Quality Heatmap** | Multi-track view of operation types and energies |
| **Text Alignment** | Traditional side-by-side alignment with match indicators |

## License

MIT — see [LICENSE](LICENSE).
