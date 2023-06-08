#!/usr/bin/env python3
#
# Name: find_k_star.py
# Description: This script finds k star given the summarized K x Xs drawn per
#              order-MEM file
# Date: May 9, 2023
#
# Run:
#   ./find_k_star.py -f summary.txt -d 4 > file.out
#
# Example data:
#   summary = /home/stephen/Documents/projects/langmead_lab/omem/workflow/results/summary.txt
#   num_docs = 4

import numpy as np
import argparse
import os

def find_k_star(path, num_docs):
    ''' Find k* for each order MEM. '''
    x_array = np.fromfile(path, dtype=int, sep='\t')
    x_mat = np.reshape(x_array, (-1, num_docs+1))
    K_s = x_mat[:,0]
    for doc_idx in range(1, num_docs+1):
        sub_mat = x_mat[:,doc_idx]
        print('order-'+str(doc_idx)+' k*:', K_s[np.argmax(sub_mat[1:] - sub_mat[:-1]) + 1])

def parse_arguments():
    """ Parse and return the command-line arguments. """
    parser = argparse.ArgumentParser(description="Takes in summary.txt and number of documents in document array profile and outputs k* of each order-MEM to stdout.")
    parser.add_argument('-f', '--file', dest='summary_path', help='path to X summary file', required=True)
    parser.add_argument('-d', '--ndocs', dest='num_docs', help='number of documents', required=True)
    args = parser.parse_args()
    return args

def check_args(args):
    """ Check that the command-line arguments are valid. """
    if not os.path.isfile(args.summary_path):
        print("Error: the fai file does not exist.")
        exit(1)

def main(args):
    path = args.summary_path
    num_docs = int(args.num_docs)
    find_k_star(path, num_docs)

if __name__ == "__main__":
    args = parse_arguments()
    check_args(args)
    main(args)
