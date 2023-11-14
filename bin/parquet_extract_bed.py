#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import argparse
import os


def _filter_pq(in_file, query_record, query_start, query_end):
    ''' Filter parquet bed file into k-mer shadow casted intervals. '''
    e_coli_mem_bed_pq_ds = ds.dataset(in_file, format="parquet")
    interval_filter_1 = (
        (ds.field('f0') == query_record) &
        (ds.field('f1') < query_start) &
        (ds.field('f2') > query_start)
    )
    interval_filter_2 = (
        (ds.field('f0') == query_record) &
        (ds.field('f1') >= query_start) &
        (ds.field('f2') <= query_end)
    )
    interval_filter_3 = (
        (ds.field('f0') == query_record) &
        (ds.field('f1') < query_end) &
        (ds.field('f2') > query_end) &
        (ds.field('f1') >= query_start)
    )
    genome_omems_arr_1 =  np.array(
        e_coli_mem_bed_pq_ds.to_table(
            filter=interval_filter_1,
            columns=['f1', 'f2', 'f3']
        ).to_pandas(), np.uint)
    genome_omems_arr_2 =  np.array(
        e_coli_mem_bed_pq_ds.to_table(
            filter=interval_filter_2,
            columns=['f1', 'f2', 'f3']
        ).to_pandas(), np.uint)
    genome_omems_arr_3 =  np.array(
        e_coli_mem_bed_pq_ds.to_table(
            filter=interval_filter_3,
            columns=['f1', 'f2', 'f3']
        ).to_pandas(), np.uint)
    genome_omems_arr = np.concatenate([genome_omems_arr_1, genome_omems_arr_2, genome_omems_arr_3], axis=0)
    return genome_omems_arr


def filter_pq(in_file, query_record, query_start, query_end):
    ''' Filter parquet bed file into k-mer shadow casted intervals. '''
    e_coli_mem_bed_pq_ds = ds.dataset(in_file, format="parquet")
    interval_filter_1 = (
        (ds.field('f0') == query_record) &
        (ds.field('f1') <= query_start) &
        (ds.field('f2') >  query_start)
    )
    interval_filter_2 = (
        (ds.field('f0') == query_record) &
        (ds.field('f1') >  query_start) &
        (ds.field('f1') <  query_end)
    )
    genome_omems_arr_1 =  np.array(
        e_coli_mem_bed_pq_ds.to_table(
            filter=interval_filter_1,
            columns=['f1', 'f2', 'f3']
        ).to_pandas(), np.uint)
    genome_omems_arr_2 =  np.array(
        e_coli_mem_bed_pq_ds.to_table(
            filter=interval_filter_2,
            columns=['f1', 'f2', 'f3']
        ).to_pandas(), np.uint)
    genome_omems_arr = np.concatenate([genome_omems_arr_1, genome_omems_arr_2], axis=0)
    return genome_omems_arr



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
    save_to_file(res, query_record, out_file)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
