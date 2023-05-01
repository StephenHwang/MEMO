#!/bin/bash

set -euo pipefail

FAI_FILE=
DAP_FILE=
INDEX_RECORDS=
OUTPUT_DIR='.'
OUTPUT_BEDFILE=
THREADS=1

usage() {
echo \
"
omem index - index overlap order MEMs from a document array profile
Usage: omem index [options]

Basic options:
  -f [FILE]              document list fai file
  -d [FILE]              full document array
  -r [RECORDS]           fasta records to query
  -o [FILE]              output directory ['.']
  -b [FILE]              output file prefix
  -t [INT]               number threads [1]
"
exit 0
}

if [ "$#" -eq 0 ] || [ "$1" = "-h" ]; then
    usage
fi

# parse flags
while getopts "f:d:r:o:b:t:" OPTION
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
            OUTPUT_BEDFILE=$OPTARG.bed
            ;;
        t )
            THREADS=$OPTARG
            ;;
        * )
            usage
            ;;
    esac
done

# Extract overlap MEM intervals
echo -e "\nConverting document array profile to overlap order MEM intervals."
./dap_to_ms_bed.py \
  --mems \
  --overlap \
  --fai $FAI_FILE \
  --dap $DAP_FILE \
  --query $INDEX_RECORDS \
  > $OUTPUT_DIR/$OUTPUT_BEDFILE

# Sort, compress, and index
echo "Sorting, compressing, and indexing interval file."
sort -k1,1V -k2,2n -o $OUTPUT_DIR/$OUTPUT_BEDFILE $OUTPUT_DIR/$OUTPUT_BEDFILE
bgzip -f --threads $THREADS $OUTPUT_DIR/$OUTPUT_BEDFILE
tabix -p bed $OUTPUT_DIR/$OUTPUT_BEDFILE.gz

echo -e "Output index file: $OUTPUT_DIR/$OUTPUT_BEDFILE.gz"
