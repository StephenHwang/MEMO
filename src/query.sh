#!/bin/bash

set -euo pipefail

# basic params
K=12
BED_FILE=
NUM_ORDER_MEMS=
QUERY_REGION=

# output
OUTPUT_DIR='\.'

# other
SAVE_INTERMEDIATES='false'
SHOW_PROGRESS='false'
SANITY_CHECK='false'

# NUM_ORDER_MEMS=4
# BED_FILE=e_coli_ordered_mems.bed
# QUERY_CHR=NZ_CP015023.1
# QUERY_START=0
# QUERY_END=5506800


usage() {
echo \
"Usage: ./omem query [options]

Query k-mer conservation across specified genome coordinates.

Output: bed file of interval position found in order-number other documents in
        the pangenome

Basic options:
  -k INT               k-mer size [12]
  -b FILE              compressed, indexed bed file of overlap order MEMs
  -n INT               number of other (non-pivot) documents in pangenome
  -r CHR:start-end     target query region
  -s                   save intermediate files
  -p                   show progress

Output options:
  -o FILE              output directory ['/.']
"
exit 0
}

if [ "$#" -eq 0 ] || [ "$1" = "-h" ]; then
    usage
fi

# parse flags
while getopts "k:b:n:r:o:sp" OPTION
do
    case $OPTION in
        k )
            K=$OPTARG
            ;;
        b )
            BED_FILE=$OPTARG
            ;;
        n )
            NUM_ORDER_MEMS=$OPTARG
            ;;
        r )
            QUERY_REGION=$OPTARG
            ;;
        o )
            OUTPUT_DIR=$OPTARG
            ;;
        s )
            SAVE_INTERMEDIATES='true'
            ;;
        p )
            SHOW_PROGRESS='true'
            ;;
        * )
            usage
            ;;
    esac
done

# Parse QUERY_REGION into record, start, end
QUERY_START_END=$(echo $QUERY_REGION | cut -d':' -f2)
QUERY_CHR=$(echo $QUERY_REGION | cut -d':' -f1)
QUERY_START=$(echo $QUERY_START_END | cut -d'-' -f1)
QUERY_END=$(echo $QUERY_START_END | cut -d'-' -f2)

# Output files
OUT_FILE=$OUTPUT_DIR/omem_$K\mer.bed
OUT_SUMMARY_FILE=$OUTPUT_DIR/omem_$K\mer.stats

################################################################################
# Start query
echo -e "STARTING: Querying $QUERY_CHR : $QUERY_START - $QUERY_END"
echo "Output bed file: $OUT_FILE"
echo "Output summary file: $OUT_SUMMARY_FILE"

# Make query bed
echo -e "$QUERY_CHR\t$QUERY_START\t$QUERY_END" > query.bed     # whole region

# extract window
if [ "$SHOW_PROGRESS" = "true" ]; then
  echo "Extracting window"
fi
tabix $BED_FILE.gz $QUERY_CHR:$QUERY_START-$QUERY_END | \
  bedtools intersect -sorted -wa -f 1 -a "stdin" -b query.bed \
  > $OUT_FILE

# Save intermediate files
if [ "$SAVE_INTERMEDIATES" = true ] ; then
  cp $OUT_FILE $OUTPUT_DIR/order_mem_overlaps.bed
fi

################################################################################
# subtract the size of the k-mer from the end and re-assign windows
awk -v k=$K '{
  d=$3-(k-1);
  if (d < $2) {
    $3=$2;
    $2=d*(d>0)
  }
  else {
    next
  }
  }1' OFS='\t' $OUT_FILE \
  > $OUT_FILE.tmp
# }1' OFS='\t' order_mem_overlaps.bed \    # use pre-extracted region


################################################################################
# merge the for each
> $OUT_FILE     # clean out file before repopulating
for MEM_IDX in $(seq 1 $NUM_ORDER_MEMS)
do
  awk -v mem_idx=$MEM_IDX '$4==mem_idx' OFS='\t' $OUT_FILE.tmp | \
    bedtools merge | \
    sed "s/$/\t$MEM_IDX/" \
    >> $OUT_FILE
done

# Save intermediate files
if [ "$SAVE_INTERMEDIATES" = true ] ; then
  cp $OUT_FILE.tmp $OUTPUT_DIR/pre_merged_ordered_mems.bed
fi
rm $OUT_FILE.tmp

# sort the final bed file
sort -k1,1V -k2,2n -o $OUT_FILE $OUT_FILE

################################################################################
# Output some summary stats
echo -e "order\tnum_records\tbp_over_window\tper_coverage" > $OUT_SUMMARY_FILE
for MEM_IDX in $(seq 1 $NUM_ORDER_MEMS)
do
  awk -v mem_idx=$MEM_IDX '$4==mem_idx' OFS='\t' $OUT_FILE | \
    bedtools coverage -sorted -a query.bed -b "stdin" | \
    awk -v mem_idx=$MEM_IDX {'printf ("%s\t%s\t%s/%s\t%s\n", mem_idx, $4, $5, $6, $7)'} \
    >> $OUT_SUMMARY_FILE
done

################################################################################

if [ "$SHOW_PROGRESS" = "true" ]; then
  echo "DONE!"
fi

