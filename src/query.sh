#!/bin/bash

set -euo pipefail

# basic params
K=12
NUM_ORDER_MEMS=
QUERY_REGION_BED=
OUTPUT_DIR='.'
SORT_OUTPUT_BED='false'

usage() {
echo \
"
omem query - query k-mers on order MEMs of input bed file region
Usage: ./omem query [options]

Basic options:
  -k [INT]               k-mer size [12]
  -n [INT]               number of other (non-pivot) documents in pangenome
  -r [FILE]              bed file of target query region
  -o [FILE]              output directory ['/.']
  -s                     sort bed output
"
exit 0
}

if [ "$#" -eq 0 ] || [ "$1" = "-h" ]; then
    usage
fi

# parse flags
while getopts "k:n:b:r:o:s" OPTION
do
    case $OPTION in
        k )
            K=$OPTARG
            ;;
        n )
            NUM_ORDER_MEMS=$OPTARG
            ;;
        r )
            QUERY_REGION_BED=$OPTARG
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

# Output files
OUT_FILE=$OUTPUT_DIR/omem_$K\mer.bed

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

# cast shadows parallelized across each order
parallel -j $NUM_ORDER_MEMS \
  "grep {}$ $QUERY_REGION_BED | cast_awk_shadows $K " \
  ::: $(seq 1 $NUM_ORDER_MEMS) \
  > $OUT_FILE.tmp

# merge intervals for each order
echo "Merging windows"
> $OUT_FILE     # rm file before re-populating
parallel -j $NUM_ORDER_MEMS \
  "grep {}$ $OUT_FILE.tmp |
  bedtools merge -c 4 -o first" \
  ::: $(seq 1 $NUM_ORDER_MEMS) \
  >> $OUT_FILE
rm $OUT_FILE.tmp

# Sort the final bed file
if [ "$SORT_OUTPUT_BED" = true ] ; then
  echo "Sorting final bed"
  sort -k1,1V -k2,2n -o $OUT_FILE $OUT_FILE
fi

echo "Output bed file: $OUT_FILE"
