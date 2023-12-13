#!/usr/bin/env python3
#
# Run: Extract bed-intervals overlapping query region
#  ./parquet_extract_bed.py \
#    -b bed.parquet \
#    -r chr:start-end \
#    -o out.bed

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import argparse
import os


# def filter_pq(in_file, query_record, query_start, query_end):
def filter_pq_og(in_file, query_record, query_start, query_end):
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

def filter_pq_simple(in_file, query_record, query_start, query_end):
# def filter_pq(in_file, query_record, query_start, query_end):
    ''' Filter parquet bed file into k-mer shadow casted intervals. '''
    pq_ds = ds.dataset(in_file, format="parquet")
    chr_filter = (ds.field('f0') == query_record)
    simple_filter = (query_start <= (ds.field('f1'))) & ((ds.field('f1') < query_end))
    joint_filter = chr_filter & simple_filter
    return np.array(pq_ds.to_table(
             filter=joint_filter,
             columns=['f1', 'f2', 'f3']
            ).to_pandas(), np.uint)


# def filter_pq_(in_file, query_record, query_start, query_end):
def filter_pq(in_file, query_record, query_start, query_end):
    ''' Filter parquet bed file into k-mer shadow casted intervals. '''
    pq_ds = ds.dataset(in_file, format="parquet")
    chr_filter = (ds.field('f0') == query_record)
    qs_in_f_filter = ((ds.field('f1') <= query_start) & (ds.field('f2') >  query_start))
    f1_in_qs_filter = ((ds.field('f1') >  query_start) & (ds.field('f1') <  query_end))
    joint_filter = chr_filter & (qs_in_f_filter | f1_in_qs_filter)
    return pq_ds.to_table(
             filter=joint_filter,
             columns=['f1', 'f2', 'f3'],
             use_threads=True
            )
    # return np.array(pq_ds.to_table(
             # filter=joint_filter,
             # columns=['f1', 'f2', 'f3'],
             # use_threads=True
            # ).to_pandas(), np.uint)


def save_to_file(arr, query_record, out_file):
    np.savetxt(out_file, arr, delimiter="\t", fmt=query_record+'\t%1i\t%1i\t%1i')


################################################################################

def parse_arguments():
    """ Parse and return the command-line arguments. """
    parser = argparse.ArgumentParser(description="Extract and query overlap MEMs for k-mer presence/absence.")
    parser.add_argument('-b', '--pq_bed_file', dest='in_file', help='parquet bed file', required=True)
    parser.add_argument('-o', '--out_file', dest='out_file', help='output bed file', required=True)
    parser.add_argument('-r', '--genome_region', dest='genome_region', help='genome region, formatted as record:start-end', required=True)
    args = parser.parse_args()
    return args

def main(args):
    in_file = args.in_file
    out_file = args.out_file
    genome_region = args.genome_region
    query_record, start_end = genome_region.split(':')
    query_start, query_end = map(int, start_end.split('-'))
    res = filter_pq(in_file, query_record, query_start, query_end)
    # save_to_file(res, query_record, out_file)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
