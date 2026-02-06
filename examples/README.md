# Catalign Examples

This directory contains example scripts demonstrating Catalign's capabilities with real-world genomic data.

## Quick Start

### 1. Mitochondrial Genome Demo (Recommended)

The fastest way to see Catalign in action with real data:

```bash
# Install visualization dependencies
pip install catalign[viz]

# Run the demo
python examples/mito_demo.py
```

This downloads real human mitochondrial sequences (~16kb each) from NCBI and:
- Aligns a sample mitochondrial genome against the reference (rCRS)
- Generates interactive visualizations
- Produces quality metrics and reports

Output is saved to `mito_demo_output/`.

### 2. Real-World Demo with Simulated T2T Data

For a more comprehensive demo with larger sequences:

```bash
python examples/real_world_demo.py --quick
```

Options:
- `--quick` - Use simulated data (no downloads)
- `--output-dir <path>` - Specify output directory
- `--query <file> --target <file>` - Use your own FASTA files

## Generated Visualizations

After running either demo, you'll find these interactive HTML files:

| File | Description |
|------|-------------|
| `dotplot.html` | K-mer dot plot showing sequence similarity patterns |
| `energy_landscape.html` | Energy profile across the alignment |
| `identity_heatmap.html` | Local sequence identity in sliding windows |
| `quality_heatmap.html` | Multi-track quality assessment |
| `cigar_view.html` | CIGAR string visualization |
| `alignment.html` | Base-by-base alignment view |
| `*_report.html` | Combined dashboard (real_world_demo only) |

## Using Your Own Data

### Command Line

```bash
# Align two FASTA files
catalign align query.fa target.fa --output text

# Generate metrics
catalign metrics query.fa target.fa --json

# Launch interactive dashboard
catalign viz
```

### Python API

```python
from catalign import CatalignAligner, evaluate_quality, read_fasta
from catalign.metrics import compute_metrics
from catalign.viz import create_dotplot, create_energy_heatmap

# Load sequences
query = list(read_fasta("query.fa"))[0][1]
target = list(read_fasta("target.fa"))[0][1]

# Align
aligner = CatalignAligner(k=15, w=50)
alignment = aligner.align(query, target)

# Evaluate
quality = evaluate_quality(alignment, query, target)
metrics = compute_metrics(alignment, query, target, quality)

print(f"Identity: {metrics.identity:.2%}")
print(f"Quality Score: {metrics.quality_score:.1f}/100")

# Visualize
dp = create_dotplot(query, target, k=11)
dp.to_figure().write_html("dotplot.html")
```

## Data Sources

### Mitochondrial Demo
- **Reference**: NC_012920.1 (revised Cambridge Reference Sequence, rCRS)
- **Sample**: KC911424.1 (Human mitochondrial complete genome)

### T2T Reference
- **T2T-CHM13 v2.0**: Complete telomere-to-telomere human reference
- URL: https://github.com/marbl/CHM13

## Performance Notes

| Sequence Length | Approximate Time | Memory |
|----------------|------------------|--------|
| 16 kb (mito) | < 1 second | ~50 MB |
| 100 kb | 2-5 seconds | ~200 MB |
| 1 Mb | 30-60 seconds | ~1 GB |
| 10 Mb | 5-10 minutes | ~5 GB |

For very large sequences, consider:
- Increasing k-mer size (`k=21` or higher)
- Using larger minimizer windows (`w=100` or higher)
- Processing in chunks

## Troubleshooting

### "plotly not installed"
```bash
pip install catalign[viz]
```

### "streamlit not installed" (for dashboard)
```bash
pip install catalign[viz]
```

### Download failures
If NCBI downloads fail, you can manually download the sequences:
1. Go to https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1
2. Click "Send to" → "File" → "FASTA"
3. Save as `mito_reference.fa`

## Example Output

After running `mito_demo.py`, you should see output like:

```
===============================================
ALIGNMENT RESULTS
===============================================
  Alignment length: 16,571 bp
  Identity:         99.85%
  Matches:          16,546
  Mismatches:       25
  Insertions:       0
  Deletions:        0
  Quality Score:    99.8/100
  Energy:           -33,017.00
```

The visualizations will show:
- **Dot plot**: Clean diagonal line (high similarity)
- **Identity heatmap**: Mostly green with occasional yellow spots (variants)
- **Energy landscape**: Low (favorable) energy throughout
