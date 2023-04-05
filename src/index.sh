#!/bin/bash
#  Create overlap MEM interval index from document array.

set -euxo pipefail

BED_FILE=e_coli_ordered_mems.bed
FAI_FILE=input.fna.fai
DAP_FILE=full_dap.txt
INDEX_RECORDS="NZ_CP015023.1 NZ_CP015022.1"

# Extact overlap MEM intervals
echo -e "STARTING: convert dap to bed\n"
/scratch4/mschatz1/sjhwang/scripts/dap_to_ms_bed.py \
  --mems \
  --overlap \
  --fai $FAI_FILE \
  --dap $DAP_FILE \
  --query $INDEX_RECORDS \
  > $BED_FILE

# Sort, compress, and index
echo -e "Sorting, compressing, and indexing bed file\n"
sort -k1,1V -k2,2n -o $BED_FILE $BED_FILE
# TODO: adjust sort... can we only sort by
rm -f $BED_FILE.gz
bgzip $BED_FILE
tabix -p bed $BED_FILE.gz

echo "DONE!"
