#!/bin/bash

set -euo pipefail

# basic params
LOG_FILE=
INDEX_DIR=
OUTDIR=

usage() {
echo \
"
Options:
  -f [LOG_FILE]          slurm log file
  -l [LENGTH_DIR]        MONI length directory
  -i [INDEX_DIR]         MONI index directory
  -o [OUTDIR]            output directory
"
exit 0
}

if [ "$#" -eq 0 ] || [ "$1" = "-h" ]; then
    usage
fi

# parse flags
while getopts "f:l:i:o:" OPTION
do
    case $OPTION in
        f )
            LOG_FILE=$OPTARG
            ;;
        l )
            LENGTH_DIR=$OPTARG
            ;;
        i )
            INDEX_DIR=$OPTARG
            ;;
        o )
            OUTDIR=$OPTARG
            ;;
        * )
            usage
            ;;
    esac
done


HEADER=$(head -n1 $LOG_FILE | awk '{print $NF}')
LENGTHS_SIZE=$(ls -l --block-size=M $LENGTH_DIR/$HEADER.lengths | cut -d' ' -f5)
INDEX_SIZE=$(ls -l --block-size=M $INDEX_DIR/$HEADER* | grep "thrbv\|slp" | awk '{print $5}' | tr -d 'M' | paste -sd+ | bc)
echo -n "| $HEADER | $LENGTHS_SIZE |" > $OUTDIR/$HEADER.txt
cat $LOG_FILE | grep "^\s" | tail -n23 | grep "Percent\|Elapsed\|Maximum" | awk '{print $NF}' | tr '\n' '|' >> $OUTDIR/$HEADER.txt
echo -n "${INDEX_SIZE}M |" >> $OUTDIR/$HEADER.txt
cat $LOG_FILE | grep "^\s" | head -n23 | grep "Percent\|Elapsed\|Maximum" | awk '{print $NF}' | tr '\n' '|' >> $OUTDIR/$HEADER.txt
echo "" >> $OUTDIR/$HEADER.txt

