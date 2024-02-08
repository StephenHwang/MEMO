#!/usr/bin/env python3
#
# Name: concat_pos_to_chr_pos.py
# Description: This script converts concatenated position to genome (chr) position.
#
# Date: Jan 8, 2024
#
# Run:
#   ./concat_pos_to_chr_pos.py --fai input.fna.fai --pos X

import argparse
import os

def read_file(path):
    ''' Read file from path line-by-line. '''
    with open(path, 'r') as inFile:
        for line in inFile:
            yield line.strip()

def parse_fai(fai_path):
    ''' Return fai file as a list of intervals. '''
    fai_stream = read_file(fai_path)
    intervals = []
    csum = 0
    for row in fai_stream:
        header, length, *_ = row.split() # get header and length, discard rest of fai
        intervals.append((header, csum, csum := csum+int(length)))
    return intervals

def pos_to_record(pos, record_intervals):
    ''' Return document header and document start position for given position
    relative to the concatenated length of the fasta file. '''
    for record, start, end in record_intervals:
        if start <= pos < end:
            return record, start
    if pos >= end:
        raise Exception('Position beyond all intervals; ensure your fai file is from fasta of initial query.')

def get_exact_pos(pos, record_intervals):
    record, offset = pos_to_record(pos, record_intervals)
    return record, pos - offset

################################################################################

def parse_arguments():
    """ Parse and return the command-line arguments. """
    parser = argparse.ArgumentParser(description="Takes in .fai and full document array profile and converts to bed-style MEM intervals to stdout.")
    parser.add_argument('--fai', dest='fai_path', help='path to fai file', required=True)
    parser.add_argument('--pos', dest='pos', help='path to full document profile', required=True)
    args = parser.parse_args()
    return args


def main(args):
    fai_path = args.fai_path
    pos = int(args.pos)
    record_intervals = parse_fai(fai_path)
    chrm, pos = get_exact_pos(pos, record_intervals)
    print(chrm, pos)

if __name__ == "__main__":
    args = parse_arguments()
    main(args)
