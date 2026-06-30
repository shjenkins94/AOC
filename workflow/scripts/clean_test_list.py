#!/usr/bin/env python3
"""
clean_test_list.py

For now this converts the sequence names in the test list into the normalized FASTA-safe format used by the AOC pipeline and outputs:

  1) mapping CSV
  2) Test-set text file (one sequence per line)

Normalization rule:
  - replace non-alphanumeric characters with "_"
  - collapse multiple "_"
  - strip leading/trailing "_"

Inputs:
  --test-headers e.g. BDNF.test_headers.csv
  --out-map      e.g. BDNF.sequence_header_map.csv
  --out-test     e.g. BDNF.test_sequences.txt

"""

import argparse
import re
import pandas as pd
from pathlib import Path
from collections import Counter

# ---------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------


def sanitize_header(h: str) -> str:
    """Convert NCBI-ish headers into single-token HyPhy-safe IDs:
    - Replace '.' and any non-alnum with '_'
    - Collapse multiple underscores
    - Strip leading/trailing underscores
    Example:
        'NM_001709.5 Homo sapiens brain derived...' ->
        'NM_001709_5_Homo_sapiens_brain_derived...'
    """
    h = (h or "").strip().lstrip(">")
    h = re.sub(r"[^A-Za-z0-9]+", "_", h)  # spaces, punctuation -> _
    h = re.sub(r"_+", "_", h)  # collapse
    h = h.strip("_")
    return h


# ---------------------------------------------------------------------
# Build mapping
# ---------------------------------------------------------------------
def build_mapping(headers_txt: Path):

    print(f"[INFO] test_headers={headers_txt}\n")

    with headers_txt.open() as f:
        raw_headers = [line.strip() for line in f if line.strip()]

    clean_headers = [sanitize_header(h) for h in raw_headers]
    # ensure uniqueness after sanitization (very important)
    counts = Counter()
    mapped_unique = []
    dup_count = 0
    for m in clean_headers:
        counts[m] += 1
        if counts[m] == 1:
            mapped_unique.append(m)
        else:
            mapped_unique.append(f"{m}_dup{counts[m]}")
            dup_count += 1
    unique_headers = pd.Series(mapped_unique)

    map_df = pd.DataFrame(
        {"original_header": raw_headers, "mapped_header": unique_headers}
    )

    print(f"[INFO] total_rows={len(map_df)}\n")
    print(f"[INFO] duplicate headers={dup_count}\n")
    print(f"[INFO] first_mappings:\n")
    for i in range(min(5, len(map_df))):
        print(
            f"  - {map_df.iloc[i]['original_header']} -> {map_df.iloc[i]['mapped_header']}\n"
        )

    return map_df


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Build normalized header mapping CSV and Test-set list (FASTA-free)."
    )
    ap.add_argument("--test", required=True)
    ap.add_argument("--out-map", required=True)
    ap.add_argument("--out-test", required=True)
    args = ap.parse_args()

    raw_headers = Path(args.test)
    out_map = Path(args.out_map)
    out_test = Path(args.out_test)

    mapped_headers = build_mapping(raw_headers)

    mapped_headers.to_csv(out_map, index=False)
    mapped_headers["mapped_header"].to_csv(out_test, index=False, header=False)

    print(f"[OK] Mapping CSV written: {out_map}\n")
    print(f"[OK] Test list written:   {out_test}\n")


if __name__ == "__main__":
    main()
