#!/usr/bin/env python3
#
# Name: parquet_stats.py
# Description: This script outputs block size metadata. Note: this isn't the full compressed or uncompressed, does not count some overhead
# Date: Dec 17, 2023
#
# Run:
#   ./parquet_stats.py -f input.bed

import argparse
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.csv
import pandas as pd

def parquet_stats(in_path, print_row_block_stats=False):
    """ Get parquet stats. """
    parquet_file = pq.ParquetFile(in_path)

    num_row_groups = parquet_file.metadata.num_row_groups
    num_columns = parquet_file.metadata.num_columns

    print('Number row groups:', num_row_groups)
    # print('Number columns', num_columns
    # print()

    total_compressed_size = 0
    total_uncompressed_size = 0

    for row_group_idx in range(0, num_row_groups):
        row_group = parquet_file.metadata.row_group(row_group_idx)
        row_block_compressed_size = 0
        row_block_uncompressed_size = 0

        if print_row_block_stats: print('row_block:', row_group_idx+1, '/', num_row_groups)
        for col_idx in range(num_columns):
            row_group_col = row_group.column(col_idx)
            row_group_col_compressed_size = row_group_col.total_compressed_size
            row_group_col_uncompressed_size = row_group_col.total_uncompressed_size
            if print_row_block_stats: print('      row_block', row_group_idx+1, 'col', col_idx+1, 'r-c compressed:', row_group_col_compressed_size)
            if print_row_block_stats: print('      row_block', row_group_idx+1, 'col', col_idx+1, 'r-c uncompressed:', row_group_col_uncompressed_size)
            row_block_compressed_size += row_group_col_compressed_size
            row_block_uncompressed_size += row_group_col_uncompressed_size
        if print_row_block_stats: print('   row_group compressed_size:', row_block_compressed_size)
        if print_row_block_stats: print('   row_group uncompressed_size:', row_block_uncompressed_size)
        total_compressed_size += row_block_compressed_size
        total_uncompressed_size += row_block_uncompressed_size

    print('\nTotal:')
    print('total compresed size:', total_compressed_size)
    print('total uncompresed size:', total_uncompressed_size)
    print()


################################################################################

def parse_arguments():
    """ Parse and return the command-line arguments. """
    parser = argparse.ArgumentParser(description="Converts input bed file to Parquet file with optional arguments of block size and codec.")
    parser.add_argument('-f', '--file', dest='file', required=True, help='parquet bed file')
    parser.add_argument('-p', '--print_row_block_stats', dest='print_row_block_stats',
                        action='store_true', required=False, default=False, help='parquet bed file')
    args = parser.parse_args()
    return args

def main(args):
    ''' Compress input bed file. '''
    in_path = args.file
    print_row_block_stats = args.print_row_block_stats
    print('Input bed:', in_path)
    parquet_stats(in_path, print_row_block_stats)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
