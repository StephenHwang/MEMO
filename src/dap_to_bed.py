#!/usr/bin/env python3
#
# Name: dap_to_bed.py
# Description: This script converts full document profile to a bed file of
#              matching statistics for given genome.
# Date: Jun 11, 2024
#
# Run:
#   ./dap_to_bed.py --mem --order --fai input.fa.fai --dap dap.txt > file.out

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

def make_fai_dict(record_intervals):
    ''' Make dictionary of query headers and lengths. '''
    fai_dict = dict()
    for header, start, end in record_intervals:
        fai_dict[header] = end - start
    return fai_dict

def print_dap_as_ms_bed(dap_stream, record_intervals, sort_lcps):
    '''
    Parse document array profile and print the matching statistcs to stdout
    as bed-records (0-index, half-open).

    Example DAP row and output:
        pos, document_array
        [3,  13, 13, 13, 12]
        chr1 3 16 1
        chr1 3 16 2
        chr1 3 16 3
        chr1 3 15 4
    '''
    for dap_row in dap_stream:
        header, pos, doc_array = get_new_record(dap_row, record_intervals, sort_lcps)
        for annot_idx, lcp in enumerate(doc_array, start=1):
            print('\t'.join(map(str, [header, pos, pos+lcp, annot_idx])))

class print_dap_as_mem_bed:
    '''
    Parse document array profile and print the MEMs (or MEM overlap intervals) to stdout as bed-records.

    Rule constituting a MEM:
        - prev_lcp <= curr_lcp

    Example DAP row:
        pos, document_array
        [3,  13, 13, 13, 12]
    '''

    def __init__(self, dap_stream, record_intervals, print_overlaps, sort_lcps):
        ''' Initialize DAP and MEM attributes. '''
        self.dap_stream = dap_stream
        self.record_intervals = record_intervals
        self.fai_dict = make_fai_dict(record_intervals)
        self.print_overlaps = print_overlaps
        self.sort_lcps = sort_lcps
        self.prev_mem_intervals_by_order = {}

    def pos_to_record(self, pos):
        ''' Return document header and document start position for given position
        relative to the concatenated length of the fasta file. '''
        for record, start, end in self.record_intervals:    # case switch statements?
            if start <= pos < end:
                return record, start
        if pos >= end:
            raise Exception('Position beyond all intervals; ensure your fai file is from fasta of initial query.')

    def get_new_record(self, dap_row):
        ''' Return parsed dap record row. '''
        pos, *lcp_array = map(int, dap_row.split(' '))
        record, offset = self.pos_to_record(pos)
        if self.sort_lcps:
            lcp_array.sort(reverse=True)
        return record, pos - offset, lcp_array

    def overlaps(self, a, b):
        ''' Return the overlap interval beteween a and b, else None. '''
        interval_start = max(a[0], b[0])
        interval_end = min(a[1], b[1])
        if interval_end >= interval_start:        # handles bookend overlaps
            return interval_start, interval_end

    def print_interval(self, header, start, end, annot):
        ''' Print MEM or MEM overlap as a BED interval. '''
        if self.print_overlaps:  # prints the overlaps between consecutive MEMs
            # prev interval exists and consecutive MEMs overlap
            prev_interval = self.prev_mem_intervals_by_order.get(annot)
            if prev_interval and (overlap := self.overlaps(prev_interval, (start, end))):
                print('\t'.join(map(str, [header, *overlap, annot])))
            self.prev_mem_intervals_by_order[annot] = (start, end)
        else:
            print('\t'.join(map(str, [header, start, end, annot])))

    def print_current_dap_row(self, header, pos, doc_array):
        ''' Print whole dap record row as BED intervals. '''
        for annot_idx, lcp in enumerate(doc_array, start=1):
            self.print_interval(header, pos, pos+lcp, annot_idx)

    def dap_to_mem(self):
        ''' Convert DAP rows to MEMs. '''
        header_prev = None                        # Initial row MEM printing, as if transitioning records
        for dap_row_next in self.dap_stream:      # Iterate through record rows of DAP
            header_curr, pos_curr, doc_array_curr = self.get_new_record(dap_row_next)
            if header_prev == header_curr:        # same record
                for annot_idx, (lcp_prev, lcp_curr) in enumerate(zip(doc_array_prev, doc_array_curr), start=1):
                    if (lcp_prev <= lcp_curr):    # print if interval is a MEM
                        self.print_interval(header_curr, pos_curr, pos_curr+lcp_curr, annot_idx)
            else:                                 # new record; record chr end, reset and print initial new-record MEMs
                if header_prev:                   # print chr end
                    prev_chr_len = self.fai_dict[header_prev]
                    self.print_current_dap_row(header_prev, prev_chr_len, [prev_chr_len] * len(doc_array_curr))  # print end of chr BED point
                self.prev_mem_intervals_by_order = {}     # reset MEMs
                self.print_current_dap_row(header_curr, pos_curr, doc_array_curr)    # print initial MEMs
            header_prev, pos_prev, doc_array_prev = header_curr, pos_curr, doc_array_curr
        # print last chr, end of chr BED point
        prev_chr_len = self.fai_dict[header_prev]
        self.print_current_dap_row(header_prev, prev_chr_len, [prev_chr_len] * len(doc_array_curr) )


################################################################################

def parse_arguments():
    """ Parse and return the command-line arguments. """
    parser = argparse.ArgumentParser(description="Takes in .fai and full document array profile and converts to bed-style MEM intervals to stdout.")
    parser.add_argument('--fai', dest='fai_path', help='path to fai file', required=True)
    parser.add_argument('--dap', dest='dap_path', help='path to full document profile', required=True)
    parser.add_argument("--ms", action="store_true", default=False, dest="print_ms", help="Extract matching statistics and print to stdout (it can either MSs or MEMs, not both).")
    parser.add_argument("--mem", action="store_true", default=False, dest="print_mems", help="Extract MEMs and print to stdout (it can either MSs or MEMs, not both).")
    parser.add_argument("--overlap", action="store_true", default=False, dest="print_overlaps", help="extract overlap MEMs (can only be used when extracting MEMs).")
    parser.add_argument("--order", action="store_true", default=False, dest="sort_lcps", help="sort LCP row to extract order MS/MEMs.")
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
        print_dap_as_mem_bed(dap_stream, record_intervals, args.print_overlaps, args.sort_lcps).dap_to_mem()

if __name__ == "__main__":
    args = parse_arguments()
    check_args(args)
    main(args)
