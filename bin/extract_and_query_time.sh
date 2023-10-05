#!/bin/bash

set -euo pipefail

# basic params
K=12
BED_FILE=
NUM_ORDER_MEMS=
QUERY_REGION=
OUTPUT_DIR='.'
SORT_OUTPUT_BED='false'

usage() {
echo \
"
omem extract and query - extract region chr:start-end from overlap order MEM index and cast k-shadows
Usage: ./extract_and_query.sh [options]

Options:
  -b [FILE]              compressed, indexed bed file of overlap order MEMs
  -r [CHR:START-END]     target query region to extract from indexed omem bed
  -k [INT]               k-mer size [12]
  -n [INT]               number of other (non-pivot) documents in pangenome
  -o [FILE]              output directory ['/.']
  -s                     sort bed output
"
exit 0
}

if [ "$#" -eq 0 ] || [ "$1" = "-h" ]; then
    usage
fi

# parse flags
while getopts "b:r:k:n:o:s" OPTION
do
    case $OPTION in
        b )
            BED_FILE=$OPTARG
            ;;
        r )
            QUERY_REGION=$OPTARG
            ;;
        k )
            K=$OPTARG
            ;;
        n )
            NUM_ORDER_MEMS=$OPTARG
            ;;
        o )
            OUTPUT_DIR=$OPTARG
            ;;
        s )
            SORT_OUTPUT_BED='true'
            ;;
        * )
            usage
            ;;
    esac
done

################################################################################
# subtract the size of the k-mer from the end and re-assign windows
echo -e "\nCasting shadows"
cast_awk_shadows() {
  awk -v k=$1 '{
    d=$3-(k-1);
    if (d < $2) {
      $3=$2;
      $2=d*(d>0)
    }
    else {
      next
    }
    }1' OFS='\t'
}
export -f cast_awk_shadows
################################################################################

# Parse QUERY_REGION into record, start, end
QUERY_CHR=$(echo $QUERY_REGION | cut -d':' -f1)
QUERY_START_END=$(echo $QUERY_REGION | cut -d':' -f2)
QUERY_START=$(echo $QUERY_START_END | cut -d'-' -f1)
QUERY_END=$(echo $QUERY_START_END | cut -d'-' -f2)

# Start query
echo -e "\nQuerying $QUERY_CHR : $QUERY_START - $QUERY_END"
OUT_FILE=$OUTPUT_DIR/omem_olaps_${QUERY_CHR}_${QUERY_START}_${QUERY_END}.bed

# Make query bed of desired regions
echo -e "$QUERY_CHR\t$QUERY_START\t$QUERY_END" > $OUTPUT_DIR/query.bed

# timing tabix extract
start_extract=`date +%s%N`
tabix $BED_FILE $QUERY_CHR:$QUERY_START-$QUERY_END | \
  bedtools intersect -sorted -wa -f 1 -a "stdin" -b $OUTPUT_DIR/query.bed > 'out.tmp'
end_extract=`date +%s%N`
echo Extract time was `expr $end_extract - $start_extract` nanoseconds.

# extract window and cast shadows
start_cast=`date +%s%N`
cat 'out.tmp' | cast_awk_shadows $K > $OUT_FILE.tmp
end_cast=`date +%s%N`
echo Cast time was `expr $end_cast - $start_cast` nanoseconds.

# merge intervals for each order
echo "Merging windows"
> $OUT_FILE     # rm file before re-populating
start_merge=`date +%s%N`
parallel -j $NUM_ORDER_MEMS \
  "grep -w {}$ $OUT_FILE.tmp |
  bedtools merge -c 4 -o first" \
  ::: $(seq 1 $NUM_ORDER_MEMS) \
  >> $OUT_FILE
end_merge=`date +%s%N`
echo Merge time was `expr $end_merge - $start_merge` nanoseconds.

rm $OUT_FILE.tmp

# Sort the final bed file
if [ "$SORT_OUTPUT_BED" = true ] ; then
  echo "Sorting final bed"
  start_sort=`date +%s%N`
  sort -k1,1V -k2,2n -o $OUT_FILE $OUT_FILE
  end_sort=`date +%s%N`
  echo Sort time was `expr $end_sort - $start_sort` nanoseconds.
fi

echo "Output bed file: $OUT_FILE"
