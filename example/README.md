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
```

## With Docker
Here is an example of the equivalent commands using Docker.

```sh
docker run --platform linux/amd64 \
  -v $(pwd):/data \
  memo:latest \
  memo index \
    -g /data/genome_list.txt \
    -o /data \
    -p test

docker run --platform linux/amd64 \
  -v $(pwd):/data \
  memo:latest \
  memo query \
    -b /data/test.parquet \
    -k 3 \
    -n 5 \
    -r ref_1:0-20 \
    -o /data/memo_c_out.txt

docker run --platform linux/amd64 \
  -v $(pwd):/data \
  memo:latest \
  memo view \
    -i /data/memo_c_out.txt \
    -o /data/test.png \
    -n 5 \
    -b 4
```

