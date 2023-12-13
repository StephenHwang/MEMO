#!/usr/bin/env python3
#
# Name: parquet_compress_bed.py
# Description: This script compresses a bed file to a Parquet file specified with block_size and codec.
# Date: Oct 19, 2023
#
# Run:
#   ./parquet_compress_bed.py -f input.bed

import argparse
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.csv


def compress_bed(in_path, out_path, block_size=500_000_000, codec='ZSTD'):
    """ Compress bed file. """
    parse_options = pyarrow.csv.ParseOptions(delimiter='\t')
    read_options = pyarrow.csv.ReadOptions(autogenerate_column_names=True, block_size=block_size)
    convert_options = pyarrow.csv.ConvertOptions()
    convert_options.column_types = {
        'f0': pa.utf8(),
        'f1': pa.int64(),
        'f2': pa.int64(),
        'f3': pa.int64(),
      }
    writer = None
    with pyarrow.csv.open_csv(in_path,
                              read_options=read_options,
                              parse_options=parse_options,
                              convert_options=convert_options) as reader:
        for next_chunk in reader:
            if next_chunk is None:
                break
            if writer is None:
                writer = pq.ParquetWriter(out_path, next_chunk.schema, compression=codec)
            next_table = pa.Table.from_batches([next_chunk])
            writer.write_table(next_table)
    writer.close()


################################################################################

def parse_arguments():
    """ Parse and return the command-line arguments. """
    parser = argparse.ArgumentParser(description="Converts input bed file to Parquet file with optional arguments of block size and codec.")
    parser.add_argument('-f', '--file', dest='file', required=True, help='bed file')
    parser.add_argument('-b', '--block_size', dest='block_size', default=500_000_000, help='Block size of rows to take in bytes [50_000_000]')
    parser.add_argument('-c', '--codec', dest='codec', default='ZSTD', help='Compression codec [ZSTD].')
    args = parser.parse_args()
    return args

def main(args):
    ''' Compress input bed file. '''
    in_path = args.file
    out_path = in_path.rstrip(".bed") + '.parquet'
    print('Input bed:', in_path)
    print('Output parquet:', out_path)
    print('Block size (bytes):', args.block_size)
    print('Code:', args.codec)
    compress_bed(in_path, out_path, block_size=int(args.block_size), codec=args.codec)
    print('DONE!')


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
