#!/usr/bin/env bash
set -e

echo "Running AOC installation test..."

# Clean previous test
rm -rf tests/results
rm -rf tests/logs

# Run minimal workflow
snakemake \
  --cores 1 \
  --directory tests

# Basic checks
#if [ ! -f tests/tmp_output/summary/run_manifest.csv ]; then
#  echo "Installation test failed."
#  exit 1
#fi

echo "Installation test passed."
