#!/bin/bash

set -euo pipefail

GENOME_LIST=
OUTPUT_DIR='.'
OUTPUT_PREFIX=
MEMBERSHIP_INDEX='false'

usage() {
echo \
"
omem index - index overlap order MEMs from a document array profile
Usage: omem index [options]

Basic options:
  -g [FILE]              document list
  -o [FILE]              output directory ['.']
  -p [FILE]              output file prefix
  -m                     make membership index
"
exit 0
}

if [ "$#" -eq 0 ] || [ "$1" = "-h" ]; then
    usage
fi

# parse flags
while getopts "g:o:p:m" OPTION
do
    case $OPTION in
        g )
            GENOME_LIST=$OPTARG
            ;;
        o )
            OUTPUT_DIR=$OPTARG
            ;;
        p )
            OUTPUT_PREFIX=$OPTARG
            ;;
        m )
            MEMBERSHIP_INDEX='true'
            ;;
        * )
            usage
            ;;
    esac
done

# prep pivot
PIVOT=$(head -n 1 $GENOME_LIST)
samtools faidx $PIVOT
mv $PIVOT.fai $OUTPUT_DIR

tail -n+2 $GENOME_LIST | while read FILE
do
  FILE_BASE=$(basename $FILE)
  # pre-process input fasta files
  seqtk seq -S $FILE > $OUTPUT_DIR/$FILE_BASE.w_rc
  samtools faidx -i $OUTPUT_DIR/$FILE_BASE.w_rc $(cat $OUTPUT_DIR/$FILE_BASE.w_rc | grep '^>' | tr -d '>') >> $OUTPUT_DIR/$FILE_BASE.w_rc
  sed -i -e '/^>/ !s/$/\$/g' $OUTPUT_DIR/$FILE_BASE.w_rc

  # Build MONI index and find MS to pivot
  echo "Build $FILE_BASE MONI index"
  moni build \
    -r $OUTPUT_DIR/$FILE_BASE.w_rc -f \
    -o $OUTPUT_DIR/$FILE_BASE.w_rc
  echo "Finding MS between pivot and $FILE_BASE"
  moni ms \
    -i $OUTPUT_DIR/$FILE_BASE.w_rc \
    -p $PIVOT \
    -o $OUTPUT_DIR/$FILE_BASE.w_rc

  # make vertical
  cat $OUTPUT_DIR/$FILE_BASE.w_rc.lengths | grep -v '^>' | tr ' ' '\n' | grep . > $OUTPUT_DIR/$FILE_BASE.w_rc.lengths.vert
done

# Make DAP
paste -d ' ' $(tail -n+2 $GENOME_LIST | xargs -I {} basename {} | sed 's/$/.w_rc.lengths.vert/' | sed "s|^|$OUTPUT_DIR\/|") | nl -v0 -w1 -s' ' > $OUTPUT_DIR/dap.txt

# From MSs to BED files
if [ "$MEMBERSHIP_INDEX" = true ] ; then
  echo "Making membership index"
  ./dap_to_ms_bed.py \
    --mem \
    --overlap \
    --fai $OUTPUT_DIR/$(basename $PIVOT.fai) \
    --dap $OUTPUT_DIR/dap.txt \
    > $OUTPUT_DIR/$OUTPUT_PREFIX.bed
else
  echo "Making conservation index"
  ./dap_to_ms_bed.py \
    --mem \
    --order \
    --overlap \
    --fai $OUTPUT_DIR/$(basename $PIVOT.fai) \
    --dap $OUTPUT_DIR/dap.txt \
    > $OUTPUT_DIR/$OUTPUT_PREFIX.bed
fi

# Compressing BED index
echo "Compressing index."
./parquet_compress_bed.py \
  -f $OUTPUT_DIR/$OUTPUT_PREFIX.bed \
  -o $OUTPUT_DIR/$OUTPUT_PREFIX.parquet

echo "DONE"
