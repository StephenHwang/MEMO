#!/bin/bash

set -euo pipefail

GENOME_LIST=
OUTPUT_DIR='.'
OUTPUT_PREFIX=
THREADS=1
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
  -t [INT]               number threads [1]
  -m                     make membership index
"
exit 0
}

if [ "$#" -eq 0 ] || [ "$1" = "-h" ]; then
    usage
fi

# parse flags
while getopts "g:o:p:t:m" OPTION
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
        t )
            THREADS=$OPTARG
            ;;
        m )
            MEMBERSHIP_INDEX='true'
            ;;
        * )
            usage
            ;;
    esac
done

echo ""
echo $GENOME_LIST
echo $OUTPUT_DIR
echo $OUTPUT_PREFIX
echo $THREADS
echo $MEMBERSHIP_INDEX
echo ""

MONI=/Users/stephenhwang/Documents/projects/langmead_lab/moni/build/moni

# Pre-process fasta files
#   for every non-pivot genome:
#     reverse complement
#     add spacer
#   for query
#     faidx

# prep pivot
PIVOT=$(head -n 1 $GENOME_LIST)
samtools faidx $PIVOT

tail -n+2 $GENOME_LIST | while read FILE
do
  echo $FILE
  seqtk seq -S $FILE > $FILE.w_rc
  samtools faidx -i $FILE.w_rc $(cat $FILE.w_rc | grep '^>' | tr -d '>') >> $FILE.w_rc
  sed -i '' -e '/^>/ !s/$/\$/g' $FILE.w_rc
  echo "Build $FILE MONI index"
  $MONI build \
    -r $FILE.w_rc -f

  # -r "${header_path}" -f \
  # -o "${index_dir}/${header}"

done




#
#echo "Finding ${query_fasta} MS on ${header} index"
#moni ms \
#  -t $THREADS \
#  -i "${index_dir}/${header}" \
#  -p ${query_fasta} \
#  -o "${index_dir}/${header}"
#
## convert MONI.lengths to vertical format and then to DAP
#cat $MONI_LENGTHS_FILE | grep -v '^>' | tr ' ' '\n' | grep . > $OUTDIR/$(basename $LENGTHS_FILE)
# }}}




#
## then paste into DAP
#paste -d ' ' $(ls *.fasta | sort) | nl -v0 -w1 -s' ' > dap.txt
#
## From MSs to BED files
#if [ "$MEMBERSHIP_INDEX" = true ] ; then
#  echo "Making membership index"
#  dap_to_ms_bed.py \
#    --mem \
#    --overlap \
#    --fai $PIVOT.fai \
#    --dap dap.txt \
#    > $OUTPUT_DIR/$OUTPUT_PREFIX.bed
#else
#  echo "Making conservation index"
#  dap_to_ms_bed.py \
#    --mem \
#    --order \
#    --overlap \
#    --fai $PIVOT.fai \
#    --dap dap.txt \
#    > $OUTPUT_DIR/$OUTPUT_PREFIX.bed
#fi
#
## Compressing BED index
#echo "Compressing index."
#parquet_compress_bed.py \
#  -f $OUTPUT_DIR/$OUTPUT_PREFIX.bed \
#  -o $OUTPUT_DIR/$OUTPUT_PREFIX.parquet
#
#echo "Index creation done."
