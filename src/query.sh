#!/bin/bash

set -euo pipefail

# basic params
K=12
BED_FILE=
NUM_ORDER_MEMS=
QUERY_REGION=
OUTPUT_DIR='.'

# other
SAVE_INTERMEDIATES='false'
SHOW_PROGRESS='false'
PRINT_SUMMARY_STATS='false'
SANITY_CHECK='false'


usage() {
echo \
"Usage: ./omem query [options]

Query k-mer conservation across specified genome coordinates.

Output: bed file of interval position found in order-number other documents in
        the pangenome

Basic options:
  -k INT               k-mer size [12]
  -n INT               number of other (non-pivot) documents in pangenome
  -b FILE              compressed, indexed bed file of overlap order MEMs
  -r CHR:start-end     target query region
  -s                   save intermediate files
  -t                   output summary stats
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
while getopts "k:n:b:r:o:stp" OPTION
do
    case $OPTION in
        k )
            K=$OPTARG
            ;;
        n )
            NUM_ORDER_MEMS=$OPTARG
            ;;
        b )
            BED_FILE=$OPTARG
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
        s )
            PRINT_SUMMARY_STATS='true'
            ;;
        p )
            SHOW_PROGRESS='true'
            ;;
        * )
            usage
            ;;
    esac
done

# Output files
OUT_FILE=$OUTPUT_DIR/omem_$K\mer.bed
OUT_SUMMARY_FILE=$OUTPUT_DIR/omem_$K\mer.stats
echo -e "STARTING:"
echo "Output bed file: $OUT_FILE"
echo "Output summary file: $OUT_SUMMARY_FILE"

if [[ $QUERY_REGION == *bed ]];      # direct file
then
  echo -e "Querying all of file: $QUERY_REGION"
  cp $OUTPUT_DIR/$QUERY_REGION $OUT_FILE
else                                    # must extract region
  # Parse QUERY_REGION into record, start, end
  QUERY_START_END=$(echo $QUERY_REGION | cut -d':' -f2)
  QUERY_CHR=$(echo $QUERY_REGION | cut -d':' -f1)
  QUERY_START=$(echo $QUERY_START_END | cut -d'-' -f1)
  QUERY_END=$(echo $QUERY_START_END | cut -d'-' -f2)

  ################################################################################
  # Start query
  echo -e "Querying $QUERY_CHR : $QUERY_START - $QUERY_END"

  # Make query bed
  echo -e "$QUERY_CHR\t$QUERY_START\t$QUERY_END" > $OUTPUT_DIR/query.bed     # whole region

  # extract window
  if [ "$SHOW_PROGRESS" = "true" ]; then
    echo "Extracting window"
  fi
  tabix $OUTPUT_DIR/$BED_FILE $QUERY_CHR:$QUERY_START-$QUERY_END | \
    bedtools intersect -sorted -wa -f 1 -a "stdin" -b $OUTPUT_DIR/query.bed \
    > $OUT_FILE

  # Save intermediate files
  if [ "$SAVE_INTERMEDIATES" = true ] ; then
    cp $OUT_FILE $OUTPUT_DIR/omem_overlaps.bed
  fi

fi


################################################################################
# subtract the size of the k-mer from the end and re-assign windows
if [ "$SHOW_PROGRESS" = "true" ]; then
  echo "Casting shadows"
fi
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


################################################################################
# merge the for each
if [ "$SHOW_PROGRESS" = "true" ]; then
  echo "Merging windows"
fi
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
if [ "$SHOW_PROGRESS" = "true" ]; then
  echo "Sorting final bed"
fi
sort -k1,1V -k2,2n -o $OUT_FILE $OUT_FILE

################################################################################
# Output some summary stats
if [ "$PRINT_SUMMARY_STATS" = true ] ; then
  if [ "$SHOW_PROGRESS" = "true" ]; then
    echo "Calculating summary stats"
  fi
  echo -e "order\tnum_records\tbp_over_window\tper_coverage" > $OUT_SUMMARY_FILE
  for MEM_IDX in $(seq 1 $NUM_ORDER_MEMS)
  do
    awk -v mem_idx=$MEM_IDX '$4==mem_idx' OFS='\t' $OUT_FILE | \
      bedtools coverage -sorted -a query.bed -b "stdin" | \
      awk -v mem_idx=$MEM_IDX {'printf ("%s\t%s\t%s/%s\t%s\n", mem_idx, $4, $5, $6, $7)'} \
      >> $OUT_SUMMARY_FILE
  done
fi

################################################################################

if [ "$SHOW_PROGRESS" = "true" ]; then
  echo "DONE!"
fi

