#!/bin/bash

set -euo pipefail

# OUTPUT_BEDFILE=e_coli_ordered_mems.bed
# FAI_FILE=input.fna.fai
# DAP_FILE=full_dap.txt
# INDEX_RECORDS="NZ_CP015023.1 NZ_CP015022.1"

FAI_FILE=
DAP_FILE=
INDEX_RECORDS=
OUTPUT_DIR='/.'
OUTPUT_BEDFILE=
SHOW_PROGRESS='false'

usage() {
echo \
"Usage: ./omem index [options]

Create a compressed and indexed overlap MEM interval bed file from a full document array.

Basic options:
  -f FILE              document list fai file
  -d FILE              full document array
  -r RECORDS           fasta records to query
  -p                   show progress

Output options:
  -o FILE              output directory ['/.']
  -b FILE              output file name
"
exit 0
}

if [ "$#" -eq 0 ] || [ "$1" = "-h" ]; then
    usage
fi

# parse flags
while getopts "l:d:r:o:b:p" OPTION
do
    case $OPTION in
        l )
            FAI_FILE=$OPTARG
            ;;
        d )
            DAP_FILE=$OPTARG
            ;;
        r )
            INDEX_RECORDS=$OPTARG
            ;;
        o )
            OUTPUT_DIR=$OPTARG
            ;;
        b )
            OUTPUT_BEDFILE=$OPTARG
            ;;
        p )
            SHOW_PROGRESS='true'
            ;;
        * )
            usage
            ;;
    esac
done

# Extract overlap MEM intervals
if [ "$SHOW_PROGRESS" = "true" ]; then
  echo "Converting document array profile to overlap order MEM intervals."
fi
dap_to_ms_bed.py \
  --mems \
  --overlap \
  --fai $FAI_FILE \
  --dap $DAP_FILE \
  --query $INDEX_RECORDS \
  > $OUT_DIR/$OUTPUT_BEDFILE

# Sort, compress, and index
if [ "$SHOW_PROGRESS" = "true" ]; then
  echo "Sorting, compressing, and indexing interval file."
fi
sort -k1,1V -k2,2n -o $OUTPUT_DIR/$OUTPUT_BEDFILE $OUTPUT_DIR/$OUTPUT_BEDFILE
bgzip -f $OUTPUT_DIR/$OUTPUT_BEDFILE
tabix -p bed $OUTPUT_DIR/$OUTPUT_BEDFILE.gz

if [ "$SHOW_PROGRESS" = "true" ]; then
  echo "Done indexing!"
fi
