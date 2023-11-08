#!/usr/bin/env python3
#
# Run: query k against moni.lengths file
#  python parse_lengths.py \
#    -f moni.lengths \
#    -k 31 \
#    > out.txt

import argparse

def file_reader(path):
    ''' Read file from path line-by-line. '''
    with open(path, 'r') as inFile:
        for line in inFile:
            yield line.strip()

################################################################################

def parse_arguments():
    """ Parse and return the command-line arguments. """
    parser = argparse.ArgumentParser(description="Reads fasta file from stdin. Output sterilzized sequence with optional rc.")
    parser.add_argument('-f', '--file', action="store", dest='file', help='Print reverse of sequence.')
    parser.add_argument('-k', '--kmer', action="store", dest='k', help='k-mer size.')
    args = parser.parse_args()
    return args

def main(args):
    ''' Sterize fasta file.

    Input:
        - stdin: fasta
    Output:
        - stdout: sterilized fasta file of records and, optionally, their
                  reverse complement
    '''
    path = args.file
    k = int(args.k)
    print('Querying k:', k)

    reader = file_reader(path)
    vals = []
    for _ in reader:
        vals.append(_)

    for length in vals:
        if int(length) >= k:
            print(1)
        else:
            print(0)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
