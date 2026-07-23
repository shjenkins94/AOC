#!/usr/bin/env bash
set -e

echo "Running AOC installation test..."

# Clean previous test
rm -rf tests/tmp_output_legacy
rm -rf tests/tmp_output_sequences

# Run minimal workflow legacy mode
snakemake \
  --cores 1 \
  --snakefile workflow/Snakefile \
  --config samples_csv=tests/data/mini_samples_only.csv outdir=tests/tmp_output_legacy

# run minimal workflow sequences mode
snakemake \
  --cores 1 \
  --snakefile workflow/Snakefile \
  --config samples_csv=tests/data/mini_samples.csv sequences_csv=tests/data/mini_sequences.csv outdir=tests/tmp_output_sequences

# Basic checks
#if [ ! -f tests/tmp_output/summary/run_manifest.csv ]; then
#  echo "Installation test failed."
#  exit 1
#fi

echo "Installation test passed."
