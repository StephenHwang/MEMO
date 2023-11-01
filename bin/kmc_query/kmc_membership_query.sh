#!/bin/bash

set -euo pipefail

# basic params
# ./kmc_membership_query.sh -m /home/shwang45/vast_sjhwang/kmc_query/kmc_database_to_hprc_genome.txt -f test.fa -r record1:8-200 -o out.txt
KMC_QUERY_SCRIPT=/scratch4/blangme2/sjhwang/omem/bin/kmc_query/kmc_query.out
KMC_DATABASE_MAPPING=
FASTA_FILE=
QUERY_REGION=
OUT_FILE=

usage() {
echo \
"
kmc conservation query - query k-mer pangenome sequence conservation in chr:start-end from a KMC database
Usage: ./kmc_conservation_query [options]

Options:
  -m [FILE]              file mapping of genome to kmc index prefix
  -f [FILE]              fasta file
  -r [CHR:START-END]     target query region to extract from indexed omem bed
  -o [FILE]              output file name ['/.']
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
while IFS=$' ' read -r genome kmc_database_path; do
    echo "Querying genome $genome at: $kmc_database_path"
    # Query KMC database
    $KMC_QUERY_SCRIPT \
        --d $kmc_database_path \
        --f _region.fa \
        > ${genome}_counts.txt
done < $KMC_DATABASE_MAPPING

paste -d ' ' $(cat $KMC_DATABASE_MAPPING | cut -d' ' -f1 | sed 's/$/_counts.txt/g') > $OUT_FILE

# Clean temporary files
rm _region.txt _region.fa *_counts.txt

echo "Done"
