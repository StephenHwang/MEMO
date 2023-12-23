#!/usr/bin/env python3
#
# Name: dap_to_bed.py
# Description: This script converts full document profile to a bed file of
#              matching statistics for given genome.
# Date: Mar 12, 2023
#
# Run:
#   ./dap_to_ms_bed.py --mems --fai input.fna.fai --dap full_dap.txt > file.out
#
# Example data:
#   /home/stephen/Documents/projects/langmead_lab/analysis/order_mems/bacteria_5/e_coli_pivot/input.fna.fai
#   /home/stephen/Documents/projects/langmead_lab/analysis/order_mems/bacteria_5/e_coli_pivot/full_dap.txt
#   NZ_CP015023.1 NZ_CP015022.1

import warnings
import argparse
import os

def read_file(path):
    ''' Read file from path line-by-line. '''
    with open(path, 'r') as inFile:
        for line in inFile:
            yield line.strip()

def parse_fai(fai_path):
    ''' Return fai file as a list of intervals. '''
    fai_stream = read_file(fai_path)
    intervals = []
    csum = 0
    for row in fai_stream:
        header, length, *_ = row.split() # get header and length, discard rest of fai
        intervals.append((header, csum, csum := csum+int(length)))
    return intervals

def pos_to_record(pos, record_intervals):
    ''' Return document header and document start position for given position
    relative to the concatenated length of the fasta file. '''
    for interval in record_intervals:
        record, start, end = interval
        if start <= pos < end:
            return record, start
    if pos >= end:
        # warnings.warn('Position beyond all intervals; ensure your fai file is from fasta of initial query.')
        raise Exception('Position beyond all intervals; ensure your fai file is from fasta of initial query.')

def get_new_record(dap_row, record_intervals, sort_lcps):
    ''' Return parsed dap record. '''
    pos, *lcp_array = map(int, dap_row.split(' '))
    record, offset = pos_to_record(pos, record_intervals)
    lcp_array = lcp_array[1:]                 # exclude self-match
    if sort_lcps:
        lcp_array.sort(reverse=True)
    return record, pos - offset, lcp_array    # return LCPs

def overlaps(a, b):
    ''' Return the overlap interval beteween a and b, else None. '''
    interval_start = max(a[0], b[0])
    interval_end = min(a[1], b[1])
    if interval_end > interval_start:
        return interval_start, interval_end

def print_dap_as_ms_bed(dap_stream, record_intervals, sort_lcps):
    '''
    Parse document array profile and print the matching statistcs to stdout
    as bed-records (0-index, half-open).

    Example DAP row and output:
        pos, document_array (self-match, rest of docs)
        [3,  65535, 13, 13, 13, 12]
        chr1 3 16 1
        chr1 3 16 2
        chr1 3 16 3
        chr1 3 15 4
    '''
    for dap_row in dap_stream:
        header, pos, doc_array = get_new_record(dap_row, record_intervals, sort_lcps)
        for annot_idx, lcp in enumerate(doc_array):
            print('\t'.join(map(str, [header, pos, pos+lcp, annot_idx+1])))

def print_dap_as_mem_bed(dap_stream, record_intervals, print_overlaps, sort_lcps):
    '''
    Parse document array profile and print the MEMs (or MEM overlap intervals)
    to stdout as bed-records.

    Example DAP row:
        pos, document_array
        [3,  65535, 13, 13, 13, 12]

    Rule constituting a mem (~or~):
        1. prev <= curr
    '''
    # get first and second record
    header_prev, pos_prev, doc_array_prev = get_new_record(next(dap_stream), record_intervals, sort_lcps)
    header_curr, pos_curr, doc_array_curr = get_new_record(next(dap_stream), record_intervals, sort_lcps)

    # track prev interval per annot_idx if printing interval overlaps
    # TODO: check this
    if print_overlaps:
        prev_ms_list = [[pos_prev, pos_prev+doc_array_prev[annot_idx]] for annot_idx in range(len(doc_array_prev))]

    # iterate through rest of records
    for dap_row_next in dap_stream:
        header_next, pos_next, doc_array_next = get_new_record(dap_row_next, record_intervals, sort_lcps)

        # At junction in between records in the bed, reset prev and curr intervals.
        # TODO: do we need to print at reset?? probably?
        if not (header_prev == header_curr == header_next):
            header_prev, pos_prev, doc_array_prev = header_curr, pos_curr, doc_array_curr
            header_curr, pos_curr, doc_array_curr = header_next, pos_next, doc_array_next
            # marks transition to new doc (instance where now the next header is in query)
            if print_overlaps:
                prev_ms_list = [[pos_next, pos_next+doc_array_next[annot_idx]] for annot_idx in range(len(doc_array_prev))]
            continue

        # print MEMs with annot
        for annot_idx, doc_arrays in enumerate(zip(doc_array_prev, doc_array_curr, doc_array_next)):
            lcp_prev, lcp_curr, lcp_next = doc_arrays
            if (lcp_prev <= lcp_curr): # check if interval is a MEM
                start_curr, end_curr = pos_curr, pos_curr + lcp_curr
                if print_overlaps:
                    start_prev, end_prev = prev_ms_list[annot_idx]   # TODO: is this suppose to be annot_idx?, verify is of sam eorder
                    overlap_interval = overlaps((start_prev, end_prev), (start_curr, end_curr))
                    if overlap_interval:     # prev interval does overlap with current interval
                        print('\t'.join(map(str, [header_curr, overlap_interval[0], overlap_interval[1], annot_idx+1])))
                    prev_ms_list[annot_idx] = [start_curr, end_curr]   # update prev MEM
                else:
                    print('\t'.join(map(str, [header_curr, start_curr, end_curr, annot_idx+1])))

        # re-assign and move onto next record
        header_prev, pos_prev, doc_array_prev = header_curr, pos_curr, doc_array_curr
        header_curr, pos_curr, doc_array_curr = header_next, pos_next, doc_array_next

################################################################################

def parse_arguments():
    """ Parse and return the command-line arguments. """
    parser = argparse.ArgumentParser(description="Takes in .fai and full document array profile and converts to bed-style MEM intervals to stdout.")
    parser.add_argument('--fai', dest='fai_path', help='path to fai file', required=True)
    parser.add_argument('--dap', dest='dap_path', help='path to full document profile', required=True)
    parser.add_argument("--ms", action="store_true", default=False, dest="print_ms", help="Extract matching statistics and print to stdout (it can either MSs or MEMs, not both).")
    parser.add_argument("--mem", action="store_true", default=False, dest="print_mems", help="Extract MEMs and print to stdout (it can either MSs or MEMs, not both).")
    parser.add_argument("--overlap", action="store_true", default=False, dest="print_overlaps", help="extract overlap MEMs (can only be used when extracting MEMs).")
    parser.add_argument("--sort_lcps", action="store_true", default=False, dest="sort_lcps", help="sort LCP row to extract order MS/MEMs.")
    args = parser.parse_args()
    return args

def check_args(args):
    """ Check that the command-line arguments are valid. """
    # Verify the files exist
    if not os.path.isfile(args.fai_path):
        raise Exception("The fai file does not exist.")
        exit(1)
    if not os.path.isfile(args.dap_path):
        raise Exception("The dap file does not exist.")
        exit(1)
    # Verify the files are the right type
    if not args.fai_path.endswith(".fai"):
        raise Exception("The fai file has the incorrect file extension.")
        exit(1)
     # Verify only one of the options are chosen: MEMs or MS
    if (args.print_ms + args.print_mems) != 1:
        raise Exception("Error: Either print MSs or MEMs, not both.")
        exit(1)
    # Verify if select overlaps, MEMs is also selected
    if args.print_overlaps and args.print_ms:
        raise Exception("Error: Can only print overlaps if printing MEMs.")
        exit(1)

def main(args):
    fai_path = args.fai_path
    dap_path = args.dap_path
    record_intervals = parse_fai(fai_path)
    dap_stream  = read_file(dap_path)
    if args.print_ms:                # print MSs
        print_dap_as_ms_bed(dap_stream, record_intervals, args.sort_lcps)
    elif args.print_mems:  # print MEMs
        print_dap_as_mem_bed(dap_stream, record_intervals, args.print_overlaps, args.sort_lcps)


if __name__ == "__main__":
    args = parse_arguments()
    check_args(args)
    main(args)
