#!/usr/bin/env python3
#
# Name: decode_parse.py
# Description: This script decodes the text ce from its dict.txt and .parse
#
# Date: Jan 23, 2024
#
# Run:
#   ./decode_parse.py --dict dict.txt --parse _.parse --window 10 > out.txt

import argparse
import os
import numpy as np

def reconstruct_text(dict_path, parse_path, wsize):
    ''' Decode text from parse and dict. '''
    # read parse
    parse = np.fromfile(parse_path, dtype='uint32')

    # read dict - removing byte chars (keeping '\x01' for splitting); can view chars by set(dict_txt)
    dict_txt = open(dict_path, 'r').read()
    dict_txt = dict_txt.replace('\x00','').replace('\x02','')
    dict_list = dict_txt.split('\x01')

    # reconstruct and print text
    print(dict_list[parse[0] - 1], end='')
    for phrase_idx in parse[1:]:
        print(dict_list[phrase_idx - 1][wsize:], end='')
    print()


################################################################################

def parse_arguments():
    """ Parse and return the command-line arguments. """
    parser = argparse.ArgumentParser(description="Takes in .fai and full document array profile and converts to bed-style MEM intervals to stdout.")
    parser.add_argument('--dict', dest='dict_path', help='path to dict.txt file; made from gzip -cdfq _.dict > dict.txt', required=True)
    parser.add_argument('--parse', dest='parse_path', help='path to .parse file', required=True)
    parser.add_argument('--wsize', dest='wsize', default=10, help='[10] parse window size', required=False)
    args = parser.parse_args()
    return args

def main(args):
    dict_path = args.dict_path
    parse_path = args.parse_path
    wsize = int(args.wsize)
    reconstruct_text(dict_path, parse_path, wsize)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
