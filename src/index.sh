#!/bin/bash

set -euo pipefail

FAI_FILE=
DAP_FILE=
INDEX_RECORDS=
OUTPUT_DIR='.'
OUTPUT_BEDFILE=
THREADS=1
SHOW_PROGRESS='false'

usage() {
echo \
"Usage: ./omem index [options]

Create a compressed and indexed overlap MEM interval bed file from a full document array.

Basic options:
  -f FILE              document list fai file
  -d FILE              full document array
  -r RECORDS           fasta records to query
  -t INT               number threads [1]
  -p                   show progress

Output options:
  -o FILE              output directory ['.']
  -b FILE              output file name
"
exit 0
}

if [ "$#" -eq 0 ] || [ "$1" = "-h" ]; then
    usage
fi

# parse flags
while getopts "f:d:r:o:b:t:p" OPTION
do
    case $OPTION in
        f )
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
        t )
            THREADS=$OPTARG
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

./dap_to_ms_bed.py \
  --mems \
  --overlap \
  --fai $FAI_FILE \
  --dap $DAP_FILE \
  --query $INDEX_RECORDS \
  > $OUTPUT_DIR/$OUTPUT_BEDFILE

# Sort, compress, and index
if [ "$SHOW_PROGRESS" = "true" ]; then
  echo "Sorting, compressing, and indexing interval file."
fi
sort -k1,1V -k2,2n -o $OUTPUT_DIR/$OUTPUT_BEDFILE $OUTPUT_DIR/$OUTPUT_BEDFILE
bgzip -f --threads $THREADS $OUTPUT_DIR/$OUTPUT_BEDFILE
tabix -p bed $OUTPUT_DIR/$OUTPUT_BEDFILE.gz

if [ "$SHOW_PROGRESS" = "true" ]; then
  echo "Done indexing!"
fi
