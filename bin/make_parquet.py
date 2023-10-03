#!/usr/bin/env python3

import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.csv

in_path = "/scratch4/blangme2/sjhwang/data/bacteria_5/mem_bed/tmp/e_coli_mem_tmp_no_header.bed"
# out_path = "/scratch4/blangme2/sjhwang/data/parquet/iterative_500m.parquet"
out_path = "tmp.parquet"

parse_options = pyarrow.csv.ParseOptions(delimiter='\t')
read_options = pyarrow.csv.ReadOptions(autogenerate_column_names=True, block_size=500_000_000)
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
            writer = pq.ParquetWriter(out_path, next_chunk.schema, compression='ZSTD')
        next_table = pa.Table.from_batches([next_chunk])
        writer.write_table(next_table)
writer.close()

