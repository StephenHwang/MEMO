# MEMO: MEM-based pangenome indexing for k-mer queries
Maximal Exact Match Ordered (MEMO) is a pangenome indexing method based on maximal exact matches (MEMs) between sequences. A single MEMO index can handle arbitrary-length-k k-mer queries over pangenomic windows. MEMO performs membership queries for per-genome k-mer presence/absence and conservation queries for the number of genomes containing the k-mers in a window. MEMO achieves smaller index size and faster queries compared to k-mer-based approaches like KMC3 and PanKmer.

MEMO relies on <a href="https://github.com/maxrossi91/moni">MONI</a> by Massimiliano Rossi for finding pairwise matching statistics between the user-selected pivot genome and all other genomes in a pangenome.


## Installation
TBD

## Usage

### Index Creation
TBD

### Querying k-mer membership and conservation
Once you have created your indexes, specify your length-k, genomic region, and th etotal number of genomes in your genome (inclusive of pivot) to k-mer conservation or  query k-mer membership `-m`.

MEMO membership query:
```sh
src/memo_query.py \
  -m \
  -b membership.parquet \
  -r chr:start-end \
  -k k \
  -n num_genomes \
  -o memo_membership.txt
```

MEMO conservation query:
```sh
src/memo_query.py \
  -b conservation.parquet \
  -r chr:start-end \
  -k k \
  -n num_genomes \
  -o memo_conservation.txt
```

## Visualizing sequence conservation
<figure>
<img src="img/memo_hla_sequence_conservation.png" alt="hprc_hla_seq_conservation"/>
<figcaption> <p align="center">31-mer sequence conservation of the Human Leucocyte Antigen locus in the HPRC pangenome.</p></figcaption>
</figure>

From the MEMO conservation query, MEMO can visualize sequence conservation:
```sh
analysis/plot_conservation.py \
  -i memo_conservation.txt \
  -o out.png \
  -n num_genomes \
  -b num_bins
```

