#!/bin/bash
#  Create overlap MEM interval index from document array.

# set -euxo pipefail
set -euo pipefail

# list of flags
K=12


## Index
#BED_FILE=
#FAI_FILE=
#DAP_FILE=
#NUM_ORDER_MEMS=
#INDEX_RECORDS=
## BED_FILE=e_coli_ordered_mems.bed
## FAI_FILE=input.fna.fai
## DAP_FILE=full_dap.txt
## NUM_ORDER_MEMS=4
## INDEX_RECORDS="NZ_CP015023.1 NZ_CP015022.1"

usage() {
echo "omem: order mem toolkit

usage: ./omem <commmand> [options]

main omem commands:
  -- index         overlap order MEM index construction from a full document array
  -- query         query overlap order MEMs of a specified region

For more commands, type omem -h
"
exit 0
}

if [ "$#" -eq 0 ] || [ "$1" = "-h" ]; then
    usage
fi



# parse flags
while getopts "k:" OPTION
do
    case $OPTION in
        k )
            K=$OPTARG
            ;;
        n )
            n="n"
            ;;
        * )
            usage
            ;;
    esac
done

echo "K=$K"











#BED_FILE=e_coli_ordered_mems.bed
#FAI_FILE=input.fna.fai
#DAP_FILE=full_dap.txt
#INDEX_RECORDS="NZ_CP015023.1 NZ_CP015022.1"
#
## Extact overlap MEM intervals
#echo -e "STARTING: convert dap to bed\n"
#/scratch4/mschatz1/sjhwang/scripts/dap_to_ms_bed.py \
#  --mems \
#  --overlap \
#  --fai $FAI_FILE \
#  --dap $DAP_FILE \
#  --query $INDEX_RECORDS \
#  > $BED_FILE
#
## Sort, compress, and index
#echo -e "Sorting, compressing, and indexing bed file\n"
#sort -k1,1V -k2,2n -o $BED_FILE $BED_FILE
## TODO: adjust sort... can we only sort by
#rm -f $BED_FILE.gz
#bgzip $BED_FILE
#tabix -p bed $BED_FILE.gz
#
#echo "DONE!"
