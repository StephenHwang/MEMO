#!/usr/bin/env python3

import numpy as np
import pyarrow.dataset as ds
import argparse
import os

def load_pq(in_file, query_record, query_start, query_end, k):
    ''' Filter parquet bed file into k-mer shadow casted intervals. '''
    # extract and filter for desired region
    e_coli_mem_bed_pq_ds = ds.dataset(in_file, format="parquet")
    interval_filter = (
        (ds.field('f0') == query_record) &
        ((query_start <= ds.field('f1')) & (ds.field('f1') <= query_end)) |
        ((query_start <= ds.field('f2')) & (ds.field('f2') <= query_end))
    )
    genome_omems_arr =  np.array(
        e_coli_mem_bed_pq_ds.to_table(
            filter=interval_filter,
            columns=['f1', 'f2', 'f3']
        ).to_pandas(), np.uint)

    # subset for candidate rows and shadow cast
    genome_omems_arr_subet = genome_omems_arr[genome_omems_arr[:,1] >= k]
    diff = genome_omems_arr_subet[:,1] - k
    genome_omems_arr_subet[:,1] = genome_omems_arr_subet[:,0]
    genome_omems_arr_subet[:,0] = diff
    return genome_omems_arr_subet[genome_omems_arr_subet[:,0] < genome_omems_arr_subet[:,1], :]

def merge_intervals(result, num_docs):
    ''' Merge intervals. '''
    merged_intervals = []
    for doc_idx in range(1, num_docs+1):
        # subset for recs of current doc
        doc_idx_intervals = result[result[:,2] == doc_idx, 0:2]
        index = 0
        for i in range(1, len(doc_idx_intervals)):
            if (doc_idx_intervals[index][1] >= doc_idx_intervals[i][0]):
                doc_idx_intervals[index][1] = max(doc_idx_intervals[index][1],
                                                  doc_idx_intervals[i][1])
            else:
                index += 1
                doc_idx_intervals[index] = doc_idx_intervals[i]
        doc_idx_intervals = doc_idx_intervals[:index+1, ]
        merged_intervals.append( np.insert(doc_idx_intervals, doc_idx_intervals.shape[1], doc_idx, axis=1) )
    # concat and sort
    merged_intervals_array = np.concatenate(merged_intervals)
    merged_sorted_array = merged_intervals_array[np.lexsort((merged_intervals_array[:,1], merged_intervals_array[:,0]))]
    return merged_sorted_array

def save_to_file(arr, query_record, out_file):
    np.savetxt(out_file, arr, delimiter="\t", fmt=query_record+'\t%1i\t%1i\t%1i')

################################################################################

def parse_arguments():
    """ Parse and return the command-line arguments. """
    parser = argparse.ArgumentParser(description="Extract and query overlap MEMs for k-mer presence/absence.")
    parser.add_argument('-b', '--bed_file', dest='in_file', help='parquet bed file', required=True)
    parser.add_argument('-o', '--out_file', dest='out_file', help='output bed file', required=True)
    parser.add_argument('-n', '--ndocs', dest='num_docs', help='number of documents', required=True)
    parser.add_argument('-k', '--kmer_size', dest='k', help='k-mer size', required=True)
    parser.add_argument('-r', '--genome_region', dest='genome_region', help='genome region, formatted as record:start-end', required=True)
    args = parser.parse_args()
    return args

def main(args):
    in_file = args.in_file
    out_file = args.out_file
    num_docs = int(args.num_docs)
    k = int(args.k) - 1
    genome_region = args.genome_region
    query_record, start_end = genome_region.split(':')
    query_start, query_end = map(int, start_end.split('-'))
    res = load_pq(in_file, query_record, query_start, query_end, k)
    mer_res = merge_intervals(res, num_docs)
    save_to_file(mer_res, query_record, out_file)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
