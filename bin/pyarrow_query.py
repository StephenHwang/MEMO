#!/usr/bin/env python3

import numpy as np
import pyarrow.dataset as ds

in_file = "iterative_500m.parquet"    # "omem_bacteria.parquet"
out_file = "shadow_casted_2.bed"

query_record = "NZ_CP015023.1"
query_start = 0
query_end = 100   #5506800
k = 12
num_docs = 4


# load parquet file into stream
e_coli_mem_bed_pq_ds = ds.dataset(in_file, format="parquet")

# extract and write to file
x = np.array(
  e_coli_mem_bed_pq_ds.to_table(
    filter=(ds.field('f0') == query_record) &
           (ds.field('f1') >= query_start) &
           (ds.field('f2') <= query_end),
      columns=['f1', 'f2', 'f3']
    ).to_pandas(),
    np.uint
  )

# filter and shadow cast
k -= 1 # adjust k buffer
def cast_shadows(row):
    start, end, doc_idx = row
    diff = end - k
    if diff < start:
        end, start = start, diff * (diff>0)
        return start, end, doc_idx
    else:
        return [np.nan] * 3
result = np.apply_along_axis(cast_shadows, axis=1, arr=x)
result = result[~np.isnan(result).any(axis=1)].astype('uint')

# merge intervals (of the same order)
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

# sort and print
fully_merged_sorted_array = np.sort(np.concatenate(merged_intervals), axis=0)

print('Saving shadow casted bed to:', out_file)
np.savetxt(out_file, fully_merged_sorted_array, delimiter="\t", fmt=query_record+'\t%1i\t%1i\t%1i')

