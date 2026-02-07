# Catalign Viewer (caliview)

A high-performance genome alignment viewer optimized for multi-scale alignment visualization.

## Features

- **Multi-scale visualization**: View alignments at base, block, region, and chromosome scales
- **Custom `.cali` format**: Pre-computed multi-scale metrics for fast browsing
- **BAM/CRAM support**: Import standard alignment formats
- **Static binary**: Easy deployment with no dependencies
- **GPU-accelerated rendering**: Smooth navigation of large genomes

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                     Caliview Architecture                    │
├─────────────────────────────────────────────────────────────┤
│                                                              │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐  │
│  │   Data Layer  │    │  Index Layer │    │  View Layer   │  │
│  ├──────────────┤    ├──────────────┤    ├──────────────┤  │
│  │ • BAM/CRAM   │───►│ • .cali file │───►│ • Tracks     │  │
│  │ • FASTA      │    │ • R-tree idx │    │ • Heatmaps   │  │
│  │ • BED/GFF    │    │ • Tile cache │    │ • Coverage   │  │
│  │ • VCF        │    │              │    │ • Metrics    │  │
│  └──────────────┘    └──────────────┘    └──────────────┘  │
│                                                              │
│  ┌──────────────────────────────────────────────────────┐   │
│  │                    Rendering Engine                    │   │
│  │  • GPU-accelerated (wgpu)                             │   │
│  │  • Level-of-detail management                         │   │
│  │  • Async tile loading                                 │   │
│  └──────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────┘
```

## File Formats

### `.cali` - Catalign Index Format

Binary format storing pre-computed multi-scale metrics:

```
Header (64 bytes):
  - Magic: "CALI" (4 bytes)
  - Version: u32
  - Flags: u32
  - Num chromosomes: u32
  - Num samples: u32
  - Index offset: u64
  - Reserved: 40 bytes

Chromosome Table:
  - Name (32 bytes, null-padded)
  - Length: u64
  - Data offset: u64
  - Data length: u64

Per-chromosome data (tiled):
  - Tile size: typically 1kb, 10kb, 100kb, 1Mb
  - Per-tile metrics:
    - Identity (f32)
    - Coverage (f32)
    - Energy (f32)
    - Gap rate (f32)
    - Quality score (f32)
    - Variant count (u32)
    - Block boundaries (variable)

Index:
  - R-tree for fast region queries
```

## Installation

### Pre-built binaries

```bash
# macOS
curl -L https://github.com/catalign/caliview/releases/latest/download/caliview-macos -o caliview
chmod +x caliview

# Linux
curl -L https://github.com/catalign/caliview/releases/latest/download/caliview-linux -o caliview
chmod +x caliview
```

### Build from source

```bash
cd catalign/viewer
cargo build --release
```

## Usage

```bash
# View a .cali file
caliview view alignment.cali

# Convert BAM to .cali
caliview convert input.bam --reference ref.fa -o output.cali

# Generate from Catalign alignment
catalign align query.fa target.fa --output-cali alignment.cali
```

## Python API

```python
from catalign.viewer import CaliFile, generate_cali

# Generate .cali from alignment
generate_cali(
    alignment_bam="aligned.bam",
    reference="reference.fa",
    output="metrics.cali",
    tile_sizes=[1000, 10000, 100000, 1000000],
)

# Read .cali file
cali = CaliFile("metrics.cali")
metrics = cali.get_region("chr1", 1000000, 2000000, resolution=1000)
```
