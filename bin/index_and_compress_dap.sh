#!/bin/bash

set -euo pipefail

DAP_FILE=/vast/blangme2/sjhwang/data/bacteria_rc_index/results/plain/full_dap.txt
FAI_FILE=/vast/blangme2/sjhwang/data/bacteria_rc_index/query/ecoli_.fa.fai
OUTPUT_BEDFILE=bacteria_overlap_mem.bed
################################################################################

DAP_TO_MS_BED_PY=/scratch4/blangme2/sjhwang/omem/src/dap_to_ms_bed.py
PARQUET_COMPRESS_BED_PY=/scratch4/blangme2/sjhwang/omem/bin/parquet_compress_bed.py
INDEX_RECORDS=$(cat $FAI_FILE | cut -f1 | paste -sd " ")

# Extract overlap MEM intervals
echo -e "\nConverting document array profile to overlap order MEM intervals."
/usr/bin/time --format='Run stats: user=%U, system=%S, elapsed=%e, CPU=%P, MemMax=%M' \
  $DAP_TO_MS_BED_PY \
    --mems \
    --overlap \
    --fai $FAI_FILE \
    --dap $DAP_FILE \
    --query $INDEX_RECORDS \
    > $OUTPUT_BEDFILE
    # --sort_lcps \     # for order MEMs

# Compress and index
echo -e "\nCompressing bed file."
/usr/bin/time --format='Run stats: user=%U, system=%S, elapsed=%e, CPU=%P, MemMax=%M' \
  $PARQUET_COMPRESS_BED_PY -f $OUTPUT_BEDFILE

