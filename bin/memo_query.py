#!/usr/bin/env python3
#
# Run: Extract bed-intervals overlapping query region
#  ./memo_query.py \
#    -b bed.parquet \
#    -r chr:start-end \
#    -n num_docs_in_pangenome \
#    -k kmer_size \
#    -o out.bed
#
# optional: -m flag for membership query

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import argparse
import os


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

def save_to_file(np_arr, query_record, out_file):
    ''' Save np array to BED file. '''
    # np.savetxt(out_file, np_arr, delimiter="\t", fmt=query_record+'\t%1i\t%1i\t%1i')
    np.savetxt(out_file, np_arr, fmt=query_record+'\t%1i\t%1i\t%1i')

def shadow_cast(genome_omems_arr, k):
    ''' Filter parquet bed file into k-mer shadow casted intervals. '''
    genome_omems_arr_subset = genome_omems_arr[genome_omems_arr[:,1] >= k]
    diff = genome_omems_arr_subset[:,1] - k
    genome_omems_arr_subset[:,1] = genome_omems_arr_subset[:,0]
    genome_omems_arr_subset[:,0] = diff
    return genome_omems_arr_subset[genome_omems_arr_subset[:,0] < genome_omems_arr_subset[:,1], :]

def merge_intervals(result, num_docs):
    ''' Merge intervals. '''
    merged_intervals = []
    for doc_idx in range(1, num_docs):
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


def adjust_margins(true_start, true_end, mb_start, mb_end, rec):
    ''' Adjust the boundaries of the output to match query interval. '''
    # assertions
    if true_start < mb_start:
        raise Exception("Error: Query does not capture interval start.")
    if len(rec) < (true_start - true_end):
        raise Exception("Error: Query does not capture interval end.")
    start_adj = true_start - mb_start
    region_len = true_end - true_start
    # print(len(rec), region_len)
    return rec[start_adj : start_adj + region_len]


def conservation_query(k, true_start, true_end, out_file, num_docs):
    ''' Get per-position document counts. '''
    mem_bed = pd.read_csv(out_file, sep='\t', header=None,
                          names=['chrm', 'start', 'end', 'order'])
    mb_start, mb_end = min(mem_bed['start']), max(mem_bed['end'])

    if true_start < mb_start:
        raise Exception("Error: Query does not capture interval start.")

    start_adj = true_start - mb_start
    region_len = true_end - true_start

    # adjust start of bed region to 0
    mem_bed['start'] -= mb_start
    mem_bed['end'] -= mb_start

    # make rec array
    rec_len = max(mb_end - mb_start + 1, true_end - true_start + 1)
    rec = [num_docs] * (rec_len + start_adj)


    for idx, values in mem_bed.iterrows():                       # iterate rows of bed file
        chrm, start, end, order = values
        start, end, order = map(int, [start, end, order])
        for idx in range(start, end):
            if order < rec[idx]:
                rec[idx] = order

    rec = rec[start_adj : start_adj + region_len]      # NOTE: queries end to end of region, even if k-mers expands past?
    print(*rec, sep='\n')


def membership_query(k, true_start, true_end, out_file, num_docs):
    ''' Get per-position presence/absence, per document. '''
    mem_bed = pd.read_csv(out_file, sep='\t', header=None,
                          names=['chrm', 'start', 'end', 'order'])
    mb_start, mb_end = min(mem_bed['start']), max(mem_bed['end'])

    # adjust the positions of the mem_bed to start at 0
    mem_bed['start'] -= mb_start
    mem_bed['end'] -= mb_start

    rec_len = max(mb_end - mb_start + 1, true_end - true_end + 1)
    rec = np.ones([num_docs, rec_len])

    for idx, values in mem_bed.iterrows():
        chrm, start, end, order = values
        start, end, order = map(int, [start, end, order])
        rec[order, start:end] = 0

    rec = adjust_margins(true_start, true_end, mb_start, mb_end, rec)
    for _ in rec.T:
        _ = map(int, _)
        print(*_, sep=' ')



################################################################################

def parse_arguments():
    """ Parse and return the command-line arguments. """
    parser = argparse.ArgumentParser(description="Extract and query overlap MEMs for k-mer presence/absence.")
    parser.add_argument('-b', '--pq_bed_file', dest='in_file', help='parquet bed file', required=True)
    parser.add_argument('-o', '--out_file', dest='out_file', help='output bed file', required=True)
    parser.add_argument('-n', '--ndocs', dest='num_docs', help='total number of genomes in the pangenome', required=True)
    parser.add_argument('-k', '--kmer_size', dest='k', help='k-mer size', required=True)
    parser.add_argument('-r', '--genome_region', dest='genome_region', help='genome region, formatted as record:start-end', required=True)
    parser.add_argument('-m', '--membership_query', dest='membership_query', action='store_true', default=False,
                        help='Perform membership query instead of conservation query')
    args = parser.parse_args()
    return args

def main(args):
    in_file = args.in_file
    out_file = args.out_file
    num_docs = int(args.num_docs)
    k = int(args.k)
    k_adj = k - 1
    genome_region = args.genome_region
    query_record, start_end = genome_region.split(':')
    query_start, query_end = map(int, start_end.split('-'))

    # filter pq and shadow cast
    genome_mems_arr = filter_pq(in_file, query_record, query_start, query_end)
    # save_to_file(genome_mems_arr, query_record, out_file + '.uncasted')
    genome_mems_arr_casted = shadow_cast(genome_mems_arr, k_adj)

    # save shadow casted bed file
    save_to_file(genome_mems_arr_casted, query_record, out_file)

    # shadow cast to per-position counts
    if args.membership_query:
        membership_query(k, query_start, query_end, out_file, num_docs)
    else:
        conservation_query(k, query_start, query_end, out_file, num_docs)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
