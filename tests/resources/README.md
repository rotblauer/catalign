# Test Resources
3. Register the test case in `tests/test_benchmarks.py`
2. Add ground truth JSON with expected alignment properties
1. Create sequence files in the appropriate subdirectory

## Adding New Test Cases

```
}
    "description": "Single SNP at position 10"
    "mutation_positions": [10],
    "expected_identity": 0.95,
    "expected_cigar": "10M1X10M",
    "target_file": "mutations/snp_single_target.fa",
    "query_file": "mutations/snp_single.fa",
{
```json

Each test case has a corresponding JSON file with expected alignment properties:

## Ground Truth Format

All synthetic data is generated with fixed random seeds for reproducibility.

```
python scripts/generate_test_data.py
```bash

To regenerate synthetic test data:

## Generating Test Data

```
└── benchmarks/         # Performance benchmark sequences
├── ground_truth/       # Expected alignment results (JSON sidecar files)
│   └── repetitive/     # Repetitive sequence patterns
│   ├── structural/     # Structural variant test cases
│   ├── indels/         # Sequences with insertions/deletions
│   ├── mutations/      # Sequences with known mutations
│   ├── identical/      # Identical sequence pairs
├── synthetic/          # Programmatically generated test sequences
resources/
```

## Directory Structure

This directory contains test data for Catalign alignment validation.

