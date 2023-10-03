#!/usr/bin/env python3

import pyarrow.dataset as ds

in_file = "e_coli_mem_bed_pq.parquet"
out_file = "parquet_extracted.bed"

query_record = "NZ_CP015023.1"
query_start = 0
query_end = 5506800

# load parquet file into stream (question, speed bonus from sort, but if no sort maybe good)
e_coli_mem_bed_pq_ds = ds.dataset(in_file, format="parquet")

# extract and write to file
e_coli_mem_bed_pq_ds.to_table(filter=
    (ds.field('record') == query_record ) &
    (ds.field('start') >= query_start) &
    (ds.field('end') <= query_end)
  ).to_pandas().to_csv(out_file, sep="\t", index=False, header=False)

