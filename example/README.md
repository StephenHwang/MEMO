# MEMO: MEM-based pangenome indexing for k-mer queries
Here is an example dataset to test MEMO.

## Index Creation
Create a MEMO conservation index.
```sh
./memo index \
  -g genome_list.txt \
  -o work \
  -p test
```

## Querying k-mer membership and conservation
Once you have created your indexes, specify your length-k `k`, genomic region `-r`, and the total number of genomes in your genome (inclusive of pivot) `-n`. Then run `memo query` for the conservation query. To run the membership query, include the `-m` flag.
```sh
./memo query \
  -b test.parquet \
  -k 3 \
  -n 5 \
  -r ref_1:0-13 \
  -o memo_c_out.txt
```

## Visualizing sequence conservation
<figure>
<img src="img/memo_hla_sequence_conservation.png" alt="hprc_hla_seq_conservation"/>
<figcaption> <p align="center">31-mer sequence conservation of the Human Leucocyte Antigen locus in the HPRC pangenome.</p></figcaption>
</figure>

From the MEMO conservation query, MEMO can visualize sequence conservation:
```sh
./memo view \
  -i memo_c_out.txt \
  -o test.png \
  -n 5 \
  -b 4
```

