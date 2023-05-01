#!/bin/bash

set -euo pipefail

# basic params
BED_FILE=
QUERY_REGION=
OUTPUT_DIR='.'

usage() {
echo \
"
omem extract - extract region chr:start-end from overlap order MEM index
Usage: omem extract [options]

Options:
  -b [FILE]              compressed, indexed bed file of overlap order MEMs
  -r [CHR:START-END]     target query region to extract from indexed omem bed
  -o [FILE]              output directory ['/.']
"
exit 0
}

if [ "$#" -eq 0 ] || [ "$1" = "-h" ]; then
    usage
fi

# parse flags
while getopts "b:r:o:" OPTION
do
    case $OPTION in
        b )
            BED_FILE=$OPTARG
            ;;
        r )
            QUERY_REGION=$OPTARG
            ;;
        o )
            OUTPUT_DIR=$OPTARG
            ;;
        * )
            usage
            ;;
    esac
done

# Parse QUERY_REGION into record, start, end
QUERY_CHR=$(echo $QUERY_REGION | cut -d':' -f1)
QUERY_START_END=$(echo $QUERY_REGION | cut -d':' -f2)
QUERY_START=$(echo $QUERY_START_END | cut -d'-' -f1)
QUERY_END=$(echo $QUERY_START_END | cut -d'-' -f2)

# Start query
echo -e "\nQuerying $QUERY_CHR : $QUERY_START - $QUERY_END"
OUT_FILE=$OUTPUT_DIR/omem_olaps_${QUERY_CHR}_${QUERY_START}_${QUERY_END}.bed

# Make query bed
echo -e "$QUERY_CHR\t$QUERY_START\t$QUERY_END" > $OUTPUT_DIR/query.bed

# extract window
tabix $OUTPUT_DIR/$BED_FILE $QUERY_CHR:$QUERY_START-$QUERY_END | \
  bedtools intersect -sorted -wa -f 1 -a "stdin" -b $OUTPUT_DIR/query.bed \
  > $OUT_FILE
rm $OUTPUT_DIR/query.bed

echo -e "Output order MEM overlaps file: $OUT_FILE"
