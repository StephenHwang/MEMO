#!/bin/bash

set -euo pipefail

# basic params
LENGTHS_FILE=
INDEX_DIR=
OUTDIR=

usage() {
echo \
"
Options:
  -f [LENGTHS_FILE]      lengths file
  -o [OUTDIR]            output directory
"
exit 0
}

if [ "$#" -eq 0 ] || [ "$1" = "-h" ]; then
    usage
fi

# parse flags
while getopts "f:o:" OPTION
do
    case $OPTION in
        f )
            LENGTHS_FILE=$OPTARG
            ;;
        o )
            OUTDIR=$OPTARG
            ;;
        * )
            usage
            ;;
    esac
done

# echo $LENGTHS_FILE
# echo $OUTDIR/$(basename $LENGTHS_FILE)
cat $LENGTHS_FILE | grep -v '^>' | tr ' ' '\n' | grep . > $OUTDIR/$(basename $LENGTHS_FILE)

