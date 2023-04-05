#!/bin/bash
# Query region:
#   Given a region, query all kmers of specified K in that window
#
#   Output: bed file (inclusive end) with intervals where a kmer starting there is
#           unique for that order MEM

set -euo pipefail
# set -euxo pipefail

# query
BED_FILE=e_coli_ordered_mems.bed
NUM_ORDER_MEMS=4
QUERY_CHR=NZ_CP015023.1
QUERY_START=0
QUERY_END=5506800
SAVE_INTERMEDIATES=true
SANITY_CHECK=true

# output
#K=$1
K=12
OUT_DIR=out
OUT_FILE=$OUT_DIR/order_mem_Xs_k$K.bed
OUT_SUMMARY_FILE=$OUT_DIR/order_mem_Xs_k$K.stats
echo "Output bed file: $OUT_FILE"
echo "Output summary file: $OUT_SUMMARY_FILE"


################################################################################
# Query to mems
echo -e "STARTING: Querying $QUERY_CHR : $QUERY_START - $QUERY_END"

################################################################################
# Make query bed
echo -e "$QUERY_CHR\t$QUERY_START\t$QUERY_END" > query.bed     # whole region

################################################################################
echo -e "          Extracting window"
# extract window
tabix $BED_FILE.gz $QUERY_CHR:$QUERY_START-$QUERY_END | \
  bedtools intersect -sorted -wa -f 1 -a "stdin" -b query.bed \
  > $OUT_FILE

# Save intermediate files
if [ "$SAVE_INTERMEDIATES" = true ] ; then
  cp $OUT_FILE $OUT_DIR/order_mem_overlaps.bed
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


# Save intermediate files
if [ "$SAVE_INTERMEDIATES" = true ] ; then
  cp $OUT_FILE.tmp $OUT_DIR/pre_merged_ordered_mems.bed
fi

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
rm $OUT_FILE.tmp

################################################################################
# sort the final bed file
sort -k1,1V -k2,2n -o $OUT_FILE $OUT_FILE

################################################################################
# Output some summary stats
# - Output windows mark Xs - positions where a kmer starting there
#                              is unique for that order MEM
# - Interpretation:
#   - for an order:
#     - % of pivot genome that is unique for that order statistic
#   - sum of orders:
#     - (100% - %) = % of pivot genome that is not unique
#     - individual docs don't sum to all due to overlap of Xs between orders

echo -e "order\tnum_records\tbp_over_window\tper_coverage" > $OUT_SUMMARY_FILE

# individual order MEMs
for MEM_IDX in $(seq 1 $NUM_ORDER_MEMS)
do
  awk -v mem_idx=$MEM_IDX '$4==mem_idx' OFS='\t' $OUT_FILE | \
    bedtools coverage -sorted -a query.bed -b "stdin" | \
    awk -v mem_idx=$MEM_IDX {'printf ("%s\t%s\t%s/%s\t%s\n", mem_idx, $4, $5, $6, $7)'} \
    >> $OUT_SUMMARY_FILE
done


################################################################################
# Assert (all_order_Xs == last_order_mem)
if [ "$SANITY_CHECK" = true ] ; then
  ALL_ORDER_X=$(bedtools coverage -sorted -a query.bed -b $OUT_FILE | cut -f5)
  HIGHEST_ORDER_X=$(awk -v mem_idx=$NUM_ORDER_MEMS '$4==mem_idx' OFS='\t' $OUT_FILE | \
                     bedtools coverage -sorted -a query.bed -b "stdin" | \
                     cut -f5)
  if [ $ALL_ORDER_X != $HIGHEST_ORDER_X ]; then
    echo "ERROR: Lower order overlap MEMs are not contained in higher order MEMs."
    exit 1
  else
    echo "Passes order overlap MEM containment test."
  fi
fi


################################################################################
echo "DONE!"
