#!/usr/bin/env python3
"""
clean_test_list.py

This script takes in an fasta MSA file and optionally a test headers txt file and outputs:

    1. a mapping CSV of original headers to sanitized headers (empty if no test headers provided)
    2. a file with the sanitized test headers (empty if no test headers provided)
    3. a file with the test headers that are present in the MSA (empty if no test headers provided or none are present in the MSA)

Normalization rule:
  - replace non-alphanumeric characters with "_"
  - collapse multiple "_"
  - strip leading/trailing "_"

Inputs:
  --msa          e.g. BDNF.fasta
  --test-headers e.g. BDNF.test_headers.csv
  --out-map      e.g. BDNF.sequence_header_map.csv
  --out-test     e.g. BDNF.test_sequences.txt
  --out-test-in-msa e.g. BDNF.test_sequences_in_msa.txt
"""

import argparse
import re
import pandas as pd
from pathlib import Path
from Bio import SeqIO
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
    ap.add_argument("--msa", required=True)
    ap.add_argument("--test-headers", required=False)
    ap.add_argument("--out-map", required=True)
    ap.add_argument("--out-test", required=True)
    ap.add_argument("--out-test-in-msa", required=True)
    args = ap.parse_args()

    msa = Path(args.msa)
    raw_headers = Path(args.test_headers) if args.test_headers else None
    out_map = Path(args.out_map)
    out_test = Path(args.out_test)
    out_test_in_msa = Path(args.out_test_in_msa)

    if raw_headers is not None:
        mapped_headers = build_mapping(raw_headers)
        records = list(SeqIO.parse(msa, "fasta"))
        if not records:
            raise ValueError(f"[ERROR] Empty MSA: {msa}")
        msa_headers = {record.id for record in records}

        print(f"[INFO] MSA headers: {msa_headers}\n")
        print(f"[INFO] Test headers: {mapped_headers['mapped_header']}\n")
        
        present_in_msa = mapped_headers[mapped_headers["mapped_header"].isin(msa_headers)]

        mapped_headers.to_csv(out_map, index=False)
        print(f"[OK] Mapping CSV written: {out_map}\n")
        mapped_headers["mapped_header"].to_csv(out_test, index=False, header=False)
        print(f"[OK] Test list written:   {out_test}\n")
        present_in_msa["mapped_header"].to_csv(out_test_in_msa, index=False, header=False)
        print(f"[OK] Test list in MSA written:   {out_test_in_msa}\n")
    else:
        # if no test headers provided, write empty files
        print(f"[INFO] No test headers provided; writing empty output files.\n")
        out_map.touch()
        out_test.touch()
        out_test_in_msa.touch()


if __name__ == "__main__":
    main()
