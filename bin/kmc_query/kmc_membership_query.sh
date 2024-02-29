#!/bin/bash

set -euo pipefail

# basic params
# ./kmc_membership_query.sh -m genome_to_kmc_database_map.txt -f test.fa -r chr:1-10 -o out.txt
KMC_QUERY_SCRIPT=/scratch4/blangme2/sjhwang/omem/bin/kmc_query/kmc_query.out
KMC_DATABASE_MAPPING=
FASTA_FILE=
QUERY_REGION=
OUT_FILE=

usage() {
echo \
"
kmc membership query - query k-mer pangenome sequence membership in chr:start-end from a list of KMC databases
Usage: ./kmc_membership_query [options]

Options:
  -m [FILE]              file containing paths to each kmc index prefix, one path per line
  -f [FILE]              fasta file
  -r [CHR:START-END]     target query region
  -o [FILE]              output file

Note:
  - Region is 1-indexed, inclusive [START-END] (to match the MEMO region START%+1, END+K-1)
  - Uses samtools to extract the query region
"
exit 0
}

if [ "$#" -eq 0 ] || [ "$1" = "-h" ]; then
    usage
fi

# parse flags
while getopts "m:f:r:o:" OPTION
do
    case $OPTION in
        m )
            KMC_DATABASE_MAPPING=$OPTARG
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

# Use a while loop to read and query each KMC database
COUNT=0
echo -n '' > _order_file.txt
while IFS=$' ' read -r KMC_DATABASE_PATH; do
    echo "Querying genome $COUNT at: $KMC_DATABASE_PATH"
    echo ${COUNT}_counts.txt >> _order_file.txt
    # Query KMC database
    $KMC_QUERY_SCRIPT \
        --d $KMC_DATABASE_PATH \
        --f _region.fa \
        > ${COUNT}_counts.txt
    COUNT=$((COUNT+1))
done < $KMC_DATABASE_MAPPING
paste -d ' ' $(cat _order_file.txt) > $OUT_FILE

# Clean temporary files
cat _order_file.txt | xargs rm
rm _region.txt _region.fa _order_file.txt

echo "Done"
