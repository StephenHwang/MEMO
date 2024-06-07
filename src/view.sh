#!/bin/bash

set -euo pipefail

# basic params
INPUT_FILE=
OUTPUT_FILE=
NUM_GENOMES=
NUM_BINS=500
DPI=600


usage() {
echo \
"
MEMO query - query k-mer membership or conservation on pivot genome region
Usage: ./memo query [options]

Basic options:
  -i [FILE]              input conservation.out from MEMO conservation query
  -o [FILE]              output plot.png
  -n [INT]               total number of documents in pangenome (include the pivot)
  -b [INT]               number of genomic bins to visualize conservation [500]
  -d [INT]               plot DPI [600]
"
exit 0
}

if [ "$#" -eq 0 ] || [ "$1" = "-h" ]; then
    usage
fi

# parse flags
while getopts "i:o:n:b:d:" OPTION
do
    case $OPTION in
        i )
            INPUT_FILE=$OPTARG
            ;;
        o )
            OUTPUT_FILE=$OPTARG
            ;;
        n )
            NUM_GENOMES=$OPTARG
            ;;
        b )
            NUM_BINS=$OPTARG
            ;;
        d )
            DPI=$OPTARG
            ;;
        * )
            usage
            ;;
    esac
done

# Get script dir
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

echo "MEMO - plotting sequence conservation"
$SCRIPT_DIR/plot_conservation.py \
  -i $INPUT_FILE \
  -o $OUTPUT_FILE \
  -n $NUM_GENOMES \
  -b $NUM_BINS \
  -d $DPI
