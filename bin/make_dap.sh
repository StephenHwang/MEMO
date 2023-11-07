#!/bin/bash

set -euo pipefail

# basic params
LENGTHS_DIR=
OUTDIR=

usage() {
echo \
"
Options:
  -l [LENGTHS_DIR]       lengths directory
  -o [OUTDIR]            output directory
"
exit 0
}

if [ "$#" -eq 0 ] || [ "$1" = "-h" ]; then
    usage
fi

# parse flags
while getopts "l:o:" OPTION
do
    case $OPTION in
        l )
            LENGTHS_DIR=$OPTARG
            ;;
        o )
            OUTDIR=$OPTARG
            ;;
        * )
            usage
            ;;
    esac
done

ls $LENGTHS_DIR/*.lengths | xargs -I file /scratch4/blangme2/sjhwang/omem/bin/verticalize_lengths.sh -f file -o $OUTDIR

echo "Done convert to vertical"
ls $OUTDIR

echo "Start make dap"
paste -d ' ' $(ls $OUTDIR/*.lengths | sort) | nl -v0 -w1 -s' ' > full_dap.txt


