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
SORT_OUTPUT_BED='false'
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
  -s                   sort bed output
  -i                   save intermediate files
  -m                   output summary stats
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
while getopts "k:n:b:r:o:t:simp" OPTION
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
            SORT_OUTPUT_BED='true'
            ;;
        i )
            SAVE_INTERMEDIATES='true'
            ;;
        m )
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
if [ "$SHOW_PROGRESS" = "true" ]; then
  echo "Output summary file: $OUT_SUMMARY_FILE"
fi


if [[ $QUERY_REGION == *bed ]];      # direct file
then
  echo -e "Querying file: $QUERY_REGION"
  QUERY_BED_FILE=$OUTPUT_DIR/$QUERY_REGION
else                                    # must extract region
  # Parse QUERY_REGION into record, start, end
  QUERY_CHR=$(echo $QUERY_REGION | cut -d':' -f1)
  QUERY_START_END=$(echo $QUERY_REGION | cut -d':' -f2)
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
  QUERY_BED_FILE=$OUT_FILE

  # Save intermediate files
  if [ "$SAVE_INTERMEDIATES" = true ] ; then
    if [ "$SHOW_PROGRESS" = "true" ]; then
      echo -e "Saving omem overlaps to: $OUTPUT_DIR/omem_overlaps.bed"
    fi
    rm -f $OUTPUT_DIR/omem_overlaps.bed
    cp $OUT_FILE $OUTPUT_DIR/omem_overlaps.bed
  fi
fi


################################################################################
# subtract the size of the k-mer from the end and re-assign windows
if [ "$SHOW_PROGRESS" = "true" ]; then
  echo "Casting shadows"
fi

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
parallel -j $NUM_ORDER_MEMS \
  "grep {}$ $QUERY_BED_FILE | cast_awk_shadows $K " \
  ::: $(seq 1 $NUM_ORDER_MEMS) \
  > $OUT_FILE.tmp

################################################################################
# merge the for each
if [ "$SHOW_PROGRESS" = "true" ]; then
  echo "Merging windows"
fi
> $OUT_FILE     # rm file before re-populating
parallel -j $NUM_ORDER_MEMS \
  "grep {}$ $OUT_FILE.tmp |
  bedtools merge -c 4 -o first" \
  ::: $(seq 1 $NUM_ORDER_MEMS) \
  >> $OUT_FILE

# Save intermediate files
if [ "$SAVE_INTERMEDIATES" = true ] ; then
  cp $OUT_FILE.tmp $OUTPUT_DIR/pre_merged_ordered_mems.bed
fi
rm $OUT_FILE.tmp

# Sort the final bed file
if [ "$SORT_OUTPUT_BED" = true ] ; then
  if [ "$SHOW_PROGRESS" = "true" ]; then
    echo "Sorting final bed"
  fi
  sort -k1,1V -k2,2n -o $OUT_FILE $OUT_FILE
fi

################################################################################
# Output some summary stats
if [ "$PRINT_SUMMARY_STATS" = true ] ; then
  if [ "$SHOW_PROGRESS" = "true" ]; then
    echo "Calculating summary stats"
  fi
  echo -e "order\tnum_records\tbp_over_window\tper_coverage" > $OUT_SUMMARY_FILE
  for MEM_IDX in $(seq 1 $NUM_ORDER_MEMS)
  do
    cat $OUT_FILE | grep "$MEM_IDX$" | \
      bedtools coverage -sorted -a query.bed -b "stdin" | \
      awk -v mem_idx=$MEM_IDX {'printf ("%s\t%s\t%s/%s\t%s\n", mem_idx, $4, $5, $6, $7)'} \
      >> $OUT_SUMMARY_FILE
  done
fi

################################################################################

if [ "$SHOW_PROGRESS" = "true" ]; then
  echo "DONE!"
fi

