#!/usr/bin/env python3
"""Test script for SV detection."""

import sys
sys.path.insert(0, '.')

output_lines = []

def log(msg):
    output_lines.append(msg)
    print(msg)

from catalign.metrics import detect_structural_variants, detect_svs_by_sequence_comparison, StructuralVariant

class MockAlignment:
    def __init__(self, pairs):
        self.aligned_pairs = pairs

# Test 1: Perfect alignment with known large INS and DEL
log("Test 1: Large INS (400bp) and DEL (5000bp)")
log("="*60)

pairs = []

# Matches 0-30000
for i in range(30000):
    pairs.append((i, i, 'M'))

# Insertion of 400bp (in query, not in target)
for i in range(400):
    pairs.append((30000 + i, None, 'I'))

# Matches 30000-45000 (but query is offset by 400)
for i in range(15000):
    pairs.append((30400 + i, 30000 + i, 'M'))

# Deletion of 5000bp (in target, not in query)
for i in range(5000):
    pairs.append((None, 45000 + i, 'D'))

# Matches 50000-59280 (query offset by 400-5000 = -4600)
for i in range(9280):
    pairs.append((45400 + i, 50000 + i, 'M'))

log(f"Total pairs: {len(pairs)}")

alignment = MockAlignment(pairs)

query_seq = 'A' * 54680  # 59280 + 400 - 5000
target_seq = 'A' * 59280

log(f"Query length: {len(query_seq)}")
log(f"Target length: {len(target_seq)}")
log(f"Difference: {len(query_seq) - len(target_seq)}")

log("\n--- Method 1: From Alignment ---")
svs = detect_structural_variants(alignment, query_seq, target_seq, min_sv_size=50)
log(f"Found {len(svs)} SVs:")
for sv in svs:
    log(f"  {sv}")
    log(f"    Details: {sv.to_dict()}")

log("\n--- Method 2: Direct Sequence Comparison ---")
svs2 = detect_svs_by_sequence_comparison(query_seq, target_seq, min_sv_size=50)
log(f"Found {len(svs2)} SVs:")
for sv in svs2:
    log(f"  {sv}")
    log(f"    Details: {sv.to_dict()}")

log("\n" + "="*60)
log("Expected: 1 INS (~400bp at pos 30000), 1 DEL (~5000bp at pos 45000)")

# Write output to file
with open("sv_test_results.txt", "w") as f:
    f.write("\n".join(output_lines))
