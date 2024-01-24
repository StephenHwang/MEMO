#!/usr/bin/env python3
#
# Name: read_N_byte_ints.py
# Description: This script reads a binary file into ints, N-bytes at a time
#
# Date: Jan 23, 2024
#
# Run:
#   ./read_N_byte_ints.py --file binary_file --int_byte_size N > out.txt

import argparse

def read_N_byte_ints(file_path, int_byte_size):
    ''' Decode text from parse and dict. '''
    with open(file_path, "rb") as f:
        while (n_bytes := f.read(int_byte_size)):
            print(int.from_bytes(n_bytes, 'little'))


################################################################################

def parse_arguments():
    ''' Parse and return the command-line arguments. '''
    parser = argparse.ArgumentParser(description="Converts binary file of ints with byte size to stdout.")
    parser.add_argument('--file', dest='file', help='path to binary file', required=True)
    parser.add_argument('--isize', dest='int_byte_size', help='int byte size', required=True)
    args = parser.parse_args()
    return args

def main(args):
    file = args.file
    int_byte_size = int(args.int_byte_size)
    read_N_byte_ints(file, int_byte_size)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
