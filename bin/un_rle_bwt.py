#!/usr/bin/env python3
#
# Name: un_RLE_bwt.py
# Description: This script reads a binary file into ints, N-bytes at a time
#
# Date: Jan 23, 2024
#
# Run:
#   ./un_RLE_bwt.py --heads heads.txt --lens lens.txt > bwt.txt

import argparse

def lens_reader(path):
    ''' Read len from path line-by-line. '''
    with open(path, 'r') as f:
        for length in f:
            yield length.strip()

def heads_reader(path):
    ''' Read head from path char-by-char. '''
    with open(path, 'r') as f:
        while (head := f.read(1)):
            yield head.strip()

def un_rle_bwt(heads_path, lens_path):
    ''' Print the un-RLE bwt. '''
    for head, length in zip(heads_reader(heads_path), lens_reader(lens_path)):
        print(head * int(length), end='')


################################################################################

def parse_arguments():
    ''' Parse and return the command-line arguments. '''
    parser = argparse.ArgumentParser(description="Converts binary file of ints with byte size to stdout.")
    parser.add_argument('--heads', dest='heads_path', help='path to heads text file', required=True)
    parser.add_argument('--lens', dest='lens_path', help='path to lens text file', required=True)
    args = parser.parse_args()
    return args

def main(args):
    heads_path = args.heads_path
    lens_path = args.lens_path
    un_rle_bwt(heads_path, lens_path)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
