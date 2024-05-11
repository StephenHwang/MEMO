#!/bin/bash

set -euo pipefail

# basic params
PARQUET_FILE=
K=31
NUM_GENOMES=
QUERY_REGION=
OUTPUT_FILE=
MEMBERSHIP_QUERY='false'


usage() {
echo \
"
MEMO query - query k-mer membership or conservation on pivot genome region
Usage: ./memo query [options]

Basic options:
  -b [FILE]              parquet membership or conservation MEMO index
  -k [INT]               k-mer size [31]
  -n [INT]               total number of documents in pangenome (include the pivot)
  -r [CHR:START-END]     query region (0-indexed, half open '[)' coordinates)
  -o [FILE]              output file
  -m                     perform the membership query (instead of conservation query)
"
exit 0
}

if [ "$#" -eq 0 ] || [ "$1" = "-h" ]; then
    usage
fi

# parse flags
while getopts "b:k:n:r:o:m" OPTION
do
    case $OPTION in
        b )
            PARQUET_FILE=$OPTARG
            ;;
        k )
            K=$OPTARG
            ;;
        n )
            NUM_GENOMES=$OPTARG
            ;;
        r )
            QUERY_REGION=$OPTARG
            ;;
        o )
            OUTPUT_FILE=$OPTARG
            ;;
        m )
            MEMBERSHIP_QUERY='true'
            ;;
        * )
            usage
            ;;
    esac
done

if [ "$MEMBERSHIP_QUERY" = true ] ; then
  echo "MEMO - membership query"
  ./memo_query.py \
    -m \
    -b $PARQUET_FILE \
    -r $QUERY_REGION \
    -k $K \
    -n $NUM_GENOMES \
    -o $OUTPUT_FILE
else
  echo "MEMO - conservation query"
  ./memo_query.py \
    -b $PARQUET_FILE \
    -r $QUERY_REGION \
    -k $K \
    -n $NUM_GENOMES \
    -o $OUTPUT_FILE
fi
