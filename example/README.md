# MEMO: MEM-based pangenome indexing for k-mer queries
Here is an example dataset to test MEMO.

Set up files:
```sh
realpath *.fa > genome_list.txt
mkdir work
cd ../src
```

Create a MEMO conservation index:
```sh
./memo index \
  -g ../example/genome_list.txt \
  -o ../example/work \
  -p test
```

Query 3-mer conservation (`ref_1:0-20`):
```sh
./memo query \
  -b ../example/work/test.parquet \
  -k 3 \
  -n 5 \
  -r ref_1:0-20 \
  -o ../example/work/memo_c_out.txt
```

Visualize sequence conservation with 3-mers:
```sh
./memo view \
  -i ../example/work/memo_c_out.txt \
  -o test.png \
  -n 5 \
  -b 4

