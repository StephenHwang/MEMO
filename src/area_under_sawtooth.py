#!/usr/bin/env python3
#
# Name: area_under_sawtooth.py
# Description: This script finds the area under sawtooth of overapping MEM interval triangles.
# Date: Jun 15, 2023
#
# Run:
#   ./area_under_sawtooth.py \
#       -f file.bed
#       -d 4 \
#       -l 5506808
#
# Example data:
#   file = /home/stephen/Documents/projects/langmead_lab/omem/data/triangle_intervals/e_coli_ordered_mems.bed
#   num_docs = 4
#   genome_len = 5506808


import numpy as np
import argparse
import os


def fileReader(path):
    ''' Read file from path line-by-line. '''
    with open(path, 'r') as inFile:
        for line in inFile:
            yield line.strip()

def calculate_area(path, num_docs, genome_len):
    ''' Calculate area under sawtooth (genome intervals). '''
    height_matrix = np.zeros((num_docs, genome_len))
    for line in fileReader(path):
        record, *fields = line.strip().split()
        start_idx, end_idx, order = map(int, fields)
        interval_array = np.array(range(end_idx - start_idx, 0, -1))
        height_matrix[order-1, start_idx:end_idx] = np.maximum(height_matrix[order-1, start_idx:end_idx], interval_array)
    return np.sum(height_matrix, axis=1) - (np.count_nonzero(height_matrix, axis=1) / 2)

def print_vals(areas):
    ''' Print formatted areas '''
    for order, area in enumerate(areas):
        print('Order-' + str(order + 1) + ':', area)


################################################################################

def parse_arguments():
    """ Parse and return the command-line arguments. """
    parser = argparse.ArgumentParser(description="Reads bed file from stdin with arguments for number of documents in document array profile and genome len. Outputs k* of each order-MEM to stdout.")
    parser.add_argument('-f', '--file', dest='path', help='bed file', required=True)
    parser.add_argument('-d', '--ndocs', dest='num_docs', help='number of documents', required=True)
    parser.add_argument('-l', '--genome_len', dest='genome_len', help='genome len', required=True)
    args = parser.parse_args()
    return args

def main(args):
    path = args.path
    num_docs = int(args.num_docs)
    genome_len = int(args.genome_len)
    areas = calculate_area(path, num_docs, genome_len)
    print_vals(areas)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
