#!/usr/bin/env python3
#
# Name: dap_to_bed.py
# Description: This script converts full document profile to a bed file of
#              matching statistics for given query genomes.
# Date: Mar 12, 2023
#
# Run:
#   ./dap_to_ms_bed.py --mems --fai input.fna.fai --dap full_dap.txt --query NZ_CP015023.1 NZ_CP015022.1 > file.out
#
# Example data:
#   /home/stephen/Documents/projects/langmead_lab/analysis/order_mems/bacteria_5/e_coli_pivot/input.fna.fai
#   /home/stephen/Documents/projects/langmead_lab/analysis/order_mems/bacteria_5/e_coli_pivot/full_dap.txt
#   NZ_CP015023.1 NZ_CP015022.1

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
        chrm, length, *_ = row.split()
        intervals.append((chrm, csum, csum := csum+int(length)))
    return intervals

def pos_to_record(pos, record_intervals):
    ''' Return document header and document start position for given position. '''
    for interval in record_intervals:
        record, start, end = interval
        if start <= pos < end:
            return record, start
        else:
            # raise Exception('Error: position greater than all record intervals')
            continue

def get_new_record(dap_row, record_intervals, sort_lcps):
    ''' Return parsed dap record. '''
    pos, *doc_array = map(int, dap_row.split(' '))
    header_record_start = pos_to_record(pos, record_intervals)
    if header_record_start is None:
        return [None] * 3
    else:
        header, record_start = header_record_start
        if sort_lcps:
            doc_array.sort(reverse=True)
        pos -= record_start
        return header, pos, doc_array[1:]

def overlaps(a, b):
    """
    Return the overlap interval beteween a and b.
    """
    interval_start = max(a[0], b[0])
    interval_end = min(a[1], b[1])
    if interval_end > interval_start:
        return interval_start, interval_end

def print_dap_as_ms_bed(dap_stream, query, record_intervals, sort_lcps):
    '''
    Parse document array profile and print the matching statistcs to stdout
    as bed-records.

    Example DAP row:
        pos, document_array
        [3,  65535, 13, 13, 13, 12]
    '''
    for dap_row in dap_stream:
        header, pos, doc_array = get_new_record(next(dap_row), record_intervals, sort_lcps)
        if header not in query:   # skip rows not in query
            continue
        for order_ms, lcp in enumerate(doc_array, start=1):
            print('\t'.join(map(str, [header, pos, pos+lcp, order_ms])))

def print_dap_as_mem_bed(dap_stream, query, record_intervals, print_overlaps, sort_lcps):
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

    # track prev interval per order_ms if printing interval overlaps
    if print_overlaps:
        prev_ms_list = [[pos_prev, pos_prev+doc_array_prev[order_ms]] for order_ms in range(len(doc_array_prev))]

    # iterate through rest of records
    for dap_row_next in dap_stream:
        header_next, pos_next, doc_array_next = get_new_record(dap_row_next, record_intervals, sort_lcps)

        # skip if documents not in query or if comparing different documents
        if (header_next not in query) or not (header_prev == header_curr == header_next):
            # re-assign and move onto next record
            header_prev, pos_prev, doc_array_prev = header_curr, pos_curr, doc_array_curr
            header_curr, pos_curr, doc_array_curr = header_next, pos_next, doc_array_next
            # marks transition to new doc (instance where now the next header is in query)
            if print_overlaps and (header_next in query):
                prev_ms_list = [[pos_next, pos_next+doc_array_next[order_ms]] for order_ms in range(len(doc_array_prev))]
            continue

        # print order MEMs
        for order_ms, doc_arrays in enumerate(zip(doc_array_prev,
                                                  doc_array_curr,
                                                  doc_array_next),
                                              start=1):
            lcp_prev, lcp_curr, lcp_next = doc_arrays

            # check if interval is a MEM
            if (lcp_prev <= lcp_curr):
                start_curr, end_curr = pos_curr, pos_curr + lcp_curr
                if not print_overlaps:    # print normal MEM intervals
                    print('\t'.join(map(str, [header_curr, start_curr, end_curr, order_ms])))
                else:                     # print overlap MEM intervals
                    start_prev, end_prev = prev_ms_list[order_ms-1]
                    intervals = overlaps((start_prev, end_prev),
                                         (start_curr, end_curr))
                    if intervals is not None:
                        print('\t'.join(map(str, [header_curr, intervals[0], intervals[1], order_ms])))
                    prev_ms_list[order_ms-1] = [start_curr, end_curr]   # update prev MEM

        # re-assign and move onto next record
        header_prev, pos_prev, doc_array_prev = header_curr, pos_curr, doc_array_curr
        header_curr, pos_curr, doc_array_curr = header_next, pos_next, doc_array_next

################################################################################

def parse_arguments():
    """ Parse and return the command-line arguments. """
    parser = argparse.ArgumentParser(description="Takes in .fai and full document array profile and converts to bed-style MEM intervals to stdout.")
    parser.add_argument('--fai', dest='fai_path', help='path to fai file', required=True)
    parser.add_argument('--dap', dest='dap_path', help='path to full document profile', required=True)
    parser.add_argument('--query', dest="query", nargs='+', help='query regions', required=True)
    parser.add_argument("--ms", action="store_true", default=False, dest="print_ms", help="Extract matching statistics and print to stdout (it can either MSs or MEMs, not both).")
    parser.add_argument("--mems", action="store_true", default=False, dest="print_mems", help="Extract MEMs and print to stdout (it can either MSs or MEMs, not both).")
    parser.add_argument("--overlap", action="store_true", default=False, dest="print_overlaps", help="extract overlap MEMs (can only be used when extracting MEMs).")
    parser.add_argument("--sort_lcps", action="store_true", default=False, dest="sort_lcps", help="sort LCP row to extract order MS/MEMs.")
    args = parser.parse_args()
    return args

def check_args(args):
    """ Check that the command-line arguments are valid. """
    # Verify the files exist
    if not os.path.isfile(args.fai_path):
        print("Error: the fai file does not exist.")
        exit(1)
    if not os.path.isfile(args.dap_path):
        print("Error: the dap file does not exist.")
        exit(1)
    # Verify the files are the right type
    if not args.fai_path.endswith(".fai"):
        print("Error: the lengths file has the incorrect file extension.")
        exit(1)
    if not args.dap_path.endswith(".txt"):
        print("Error: the dap file has the incorrect file extension.")
        exit(1)
     # Verify only one of the options are chosen: MEMs or MS
    if (args.print_ms + args.print_mems) != 1:
        print("Error: exactly one type needs to be chosen (--ms or --mems).")
        exit(1)
    # Verify if select overlaps, MEMs is also selected
    if args.print_overlaps and args.print_ms:
        print("Error: can only print overlaps if printing MEMs.")
        exit(1)


def main(args):
    # input paths
    fai_path = args.fai_path
    dap_path = args.dap_path
    query = set(args.query)
    record_intervals = parse_fai(fai_path)
    dap_stream  = read_file(dap_path)
    if args.print_ms:                # print MSs
        print_dap_as_ms_bed(dap_stream, query, record_intervals, args.sort_lcps)
    elif args.print_mems:  # print MEMs
        print_dap_as_mem_bed(dap_stream, query, record_intervals, args.print_overlaps, args.sort_lcps)


if __name__ == "__main__":
    args = parse_arguments()
    check_args(args)
    main(args)
