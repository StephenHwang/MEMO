#!/usr/bin/env python
#
# Run: Extract bed-intervals overlapping query region
#  ./memo_query.py \
#    -b bed.parquet \
#    -r chr:start-end \
#    -k kmer_size \
#    -n num_docs_in_pangenome \
#    -o out.txt
#
# optional: -m flag for membership query

import numpy as np
import pyarrow.dataset as ds
import argparse
import os
from numba import jit

def filter_pq(in_file, query_record, query_start, query_end):
    ''' Filter parquet bed file into k-mer shadow casted intervals. '''
    pq_ds = ds.dataset(in_file, format="parquet")
    qs_in_f_filter = ((ds.field('f0') == query_record) &
                      (ds.field('f1') <= query_start) &
                      (ds.field('f2') >  query_start))
    f1_in_qs_filter = ((ds.field('f0') == query_record) &
                       (ds.field('f1') >  query_start) &
                       (ds.field('f1') <  query_end))
    qs_in_f_arr = np.array(pq_ds.to_table(
                           filter=qs_in_f_filter,
                           columns=['f1', 'f2', 'f3']
                          ).to_pandas(), np.uint)
    f1_in_qs_arr = np.array(pq_ds.to_table(
                            filter=f1_in_qs_filter,
                            columns=['f1', 'f2', 'f3']
                           ).to_pandas(), np.uint)
    return np.concatenate([qs_in_f_arr, f1_in_qs_arr], axis=0)

def save_to_bed_file(np_arr, query_record, save_file):
    ''' Save np array to output save file. '''
    np.savetxt(save_file, np_arr, fmt=query_record+'\t%1i\t%1i\t%1i')

def memo_init(mem_arr, k, true_start, true_end, num_docs, membership_query):
    ''' Initialize MEM intervals and result array. '''
    true_len = true_end - true_start
    mem_arr = mem_arr.astype('int64')
    mem_arr[:, 0:2] -= true_start  # re-center mem intervals (start, ends) to query region
    mem_arr[:, 1] -= k - 1         # then "shadow cast" k
    mem_arr[:, 0:2] = mem_arr[:, 0:2].clip(min=0, max=true_len)
    mem_arr = mem_arr[mem_arr[:,1] < mem_arr[:,0]]   # subset to valid intervals
    if membership_query: # membership
        rec = np.ones([true_len, num_docs], dtype='bool')
    else:                # conservation
        rec = np.zeros([true_len, num_docs+1], dtype='bool')
        rec[:, num_docs] = True
    return mem_arr, rec

@jit
def memo_query(mem_arr, rec, membership_query):
    ''' Peform MEMO membership or conservation query by recording k-mer casting for each doc/order. '''
    set_bit = False if membership_query else True
    for start, casted_end, order in mem_arr:
        rec[casted_end:start, order] = set_bit
    return rec

def print_res(rec, out_file, membership_query):
    ''' Output result to file. '''
    if membership_query:
        np.savetxt(out_file, rec.astype('byte'), delimiter=' ', fmt='%i')
    else:
        rec = np.argmax(rec, axis=1)
        print(*rec, sep='\n', file=open(out_file, 'w'))


################################################################################

def parse_arguments():
    """ Parse and return the command-line arguments. """
    parser = argparse.ArgumentParser(description="Extract and query overlap MEMs for k-mer presence/absence.")
    parser.add_argument('-b', '--pq_bed_file', dest='in_file', help='parquet bed file', required=True)
    parser.add_argument('-o', '--out_file', dest='out_file', help='output file', required=True)
    parser.add_argument('-n', '--ndocs', dest='num_docs', help='total number of genomes in the pangenome', required=True)
    parser.add_argument('-k', '--kmer_size', dest='k', help='k-mer size', required=True)
    parser.add_argument('-r', '--genome_region', dest='genome_region', help='genome region, formatted as chr:start-end', required=True)
    parser.add_argument('-m', '--membership_query', dest='membership_query', action='store_true', default=False,
                        help='Perform membership query instead of conservation query')
    args = parser.parse_args()
    return args

def main(args):
    in_file = args.in_file
    out_file = args.out_file
    membership_query = args.membership_query
    num_docs = int(args.num_docs)
    k = int(args.k)
    genome_region = args.genome_region
    query_record, start_end = genome_region.split(':')
    query_start, query_end = map(int, start_end.split('-'))

    # extract bed region
    genome_mems_arr = filter_pq(in_file, query_record, query_start, query_end + k)

    # query
    mem_arr, rec = memo_init(genome_mems_arr, k, query_start, query_end, num_docs, membership_query)
    rec = memo_query(mem_arr, rec, membership_query)
    print_res(rec, out_file, membership_query)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
