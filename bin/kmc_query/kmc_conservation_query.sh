#!/bin/bash

set -euo pipefail

# basic params
# ./kmc_conservation_query.sh -d /home/shwang45/vast_sjhwang/kmc_query/file -f test.fa -r record1:8-200 -o out.txt
KMC_QUERY_SCRIPT=/scratch4/blangme2/sjhwang/omem/bin/kmc_query/kmc_query.out
KMC_DATABASE_PREFIX=
FASTA_FILE=
QUERY_REGION=
OUT_FILE=

usage() {
echo \
"
kmc conservation query - query k-mer pangenome sequence conservation in chr:start-end from a KMC database
Usage: ./kmc_conservation_query [options]

Options:
  -d [FILE]              kmc index prefix
  -f [FILE]              fasta file
  -r [CHR:START-END]     target query region to extract from indexed omem bed
  -o [FILE]              output file name ['/.']

Note: Region is indexed as [START-END) by 0-index.
  ie. to query a whole contig, use: 0 to samtools faidx size
"
exit 0
}

if [ "$#" -eq 0 ] || [ "$1" = "-h" ]; then
    usage
fi

# parse flags
while getopts "d:f:r:o:" OPTION
do
    case $OPTION in
        d )
            KMC_DATABASE_PREFIX=$OPTARG
            ;;
        f )
            FASTA_FILE=$OPTARG
            ;;
        r )
            QUERY_REGION=$OPTARG
            ;;
        o )
            OUT_FILE=$OPTARG
            ;;
        * )
            usage
            ;;
    esac
done

# Make region file
echo -e "$QUERY_REGION" > _region.txt
samtools faidx $FASTA_FILE
samtools faidx $FASTA_FILE -r _region.txt > _region.fa

# Query KMC database
$KMC_QUERY_SCRIPT \
  --d $KMC_DATABASE_PREFIX \
  --f _region.fa \
  > $OUT_FILE

# Clean temporary files
rm _region.txt _region.fa

echo "Done"
