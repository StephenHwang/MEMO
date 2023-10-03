#!/bin/bash

set -euo pipefail

OUTPUT_BEDFILE=xxx
THREADS=1

# Sort, compress, and index
echo "Sorting, compressing, and indexing interval file."
sort -k1,1V -k2,2n -o $OUTPUT_BEDFILE $OUTPUT_BEDFILE
bgzip -f --threads $THREADS $OUTPUT_BEDFILE
tabix -p bed $OUTPUT_BEDFILE.gz

echo -e "Output index file: $OUTPUT_BEDFILE.gz"
