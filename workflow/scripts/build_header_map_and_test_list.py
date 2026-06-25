#!/usr/bin/env python3
"""
Essentialy the same as what build_header_map_and_test_list rule did, 
but rewritten to be a snakemake script. Stuff that checks stuff
that Snakemake already checks is removed.
"""
import sys
import os, re
import pandas as pd
from collections import Counter

import snakemake

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

with open(snakemake.log[0], "w") as L:
    sys.stderr =sys.stdout = L

    print(f"[INFO] sample={snakemake.wildcards['sample']}\n")
    print(f"[INFO] labels={snakemake.input['labels']}\n")

    df = pd.read_csv(snakemake.input['labels'])
    # normalize column names
    df.columns = [c.strip().lower() for c in df.columns]
    # heuristic column detection
    label_candidates = [
        "label",
        "group",
        "set",
        "class",
        "foreground",
        "branchset",
    ]
    header_candidates = [
        "header",
        "sequence_header",
        "seq_header",
        "id",
        "name",
        "taxon",
        "tip",
        "sequence",
    ]
    label_col = next((c for c in label_candidates if c in df.columns), None)
    header_col = next((c for c in header_candidates if c in df.columns), None)
    if label_col is None or header_col is None:
        if df.shape[1] == 2:
            label_col = df.columns[0]
            header_col = df.columns[1]
        else:
            raise ValueError(
                f"[build_header_map_and_test_list] Could not infer label/header columns.\n"
                f"Columns found: {list(df.columns)}\n"
                f"Expected something like LABEL + HEADER (or a 2-column CSV)."
            )
    labels_raw = df[label_col].astype(str).str.strip()
    headers_raw = df[header_col].astype(str).str.strip().str.lstrip(">")
    # sanitize -> mapped headers
    mapped = headers_raw.map(sanitize_header)
    # ensure uniqueness after sanitization (very important)
    counts = Counter()
    mapped_unique = []
    for m in mapped.tolist():
        counts[m] += 1
        if counts[m] == 1:
            mapped_unique.append(m)
        else:
            mapped_unique.append(f"{m}_dup{counts[m]}")
    mapped_unique = pd.Series(mapped_unique)
    # define what counts as Test
    labels_norm = labels_raw.str.lower()
    is_test = labels_norm.isin(
        {"test", "foreground", "fg", "case"}
    ) | labels_norm.str.contains(r"\btest\b", regex=True)
    test_headers_mapped = mapped_unique[is_test].dropna().tolist()
    # write map
    map_df = pd.DataFrame(
        {"original_header": headers_raw, "mapped_header": mapped_unique}
    )
    map_df.to_csv(snakemake.output["map_csv"], index=False)
    # write test list (MAPPED headers, one per line)
    with open(snakemake.output["test_txt"], "w") as out:
        out.write(
            "\n".join(test_headers_mapped)
            + ("\n" if test_headers_mapped else "")
        )
    print(f"[INFO] label_col={label_col}, header_col={header_col}\n")
    print(f"[INFO] total_rows={len(df)}\n")
    print(f"[INFO] test_rows={len(test_headers_mapped)}\n")
    print(f"[INFO] first_mappings:\n")
    for i in range(min(5, len(map_df))):
        print(
            f"  - {map_df.iloc[i]['original_header']} -> {map_df.iloc[i]['mapped_header']}\n"
        )




    # strict checks
    if os.path.getsize(output.map_csv) == 0:
        raise RuntimeError(
            f"[build_header_map_and_test_list] map_csv empty: {output.map_csv}"
        )
    if (
        snakemake.config.get("require_test_sequences", False)
        and len(test_headers_mapped) == 0
    ):
        raise RuntimeError(
            f"[build_header_map_and_test_list] No Test sequences found for sample={snakemake.wildcards['sample']} "
            f"but require_test_sequences=True. See log: {snakemake.log[0]}"
        )


