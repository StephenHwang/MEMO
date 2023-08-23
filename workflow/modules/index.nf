/*
 * Module for creating order-MEM overlap index
 */


/*
 * Run doc_pfp to extract document array profile.
 */
process MONI_EXTRACT_DAP {
  label 'moni_extract_dap'
  debug true

  input:
  path moni
  path preprocess_moni_fasta
  path document_listing

  output:
  path "*.lenghts", emit: moni_length_files
  // path "${out_prefix}_dap.txt", emit: dap
  // path "${out_prefix}.fna", emit: fna

  script:
  """
  // pre-process
  document_listing_file = file($document_listing)
  allDocs = document_listing_file.readLines()
  for( doc : allDocs ) {
    // get the header
    doc_header = doc.getSimpleName()

    ./$preprocess_moni_fasta -r \
      < ${HEADER}.fa \
      > ${HEADER}_with_rc.fa

    ./$preprocess_moni_fasta \
      < ${HEADER}.fa \
      > ${HEADER}_processed.fa

    // build index (on fasta with rc)
    ./$moni build \
      -r ${HEADER}_with_rc.fa -f \
      -o index/${HEADER}

  // held out query
  cat *.processed.fa > all_docs_processed.fa

  // then for every one, run moni ms
  allDocs = document_listing_file.readLines()
  for( doc : allDocs ) {
    doc_header = doc.getSimpleName()
    ./moni ms \
      -i $doc_header \
      -p all_docs_processed.fa \
      -o moni_ms/$doc_header

  }

  // then convert to DAP
  """
}


























/*
 * Run doc_pfp to extract document array profile.
 */
process EXTRACT_DAP {
  label 'extract_dap'
  debug true

  input:
  path doc_pfp
  path document_listing
  val out_prefix
  val doc_pivot_id

  output:
  path "${out_prefix}_dap.txt", emit: dap
  path "${out_prefix}.fna", emit: fna

  script:
  """
  ./$doc_pfp build -r \
    -f $document_listing \
    -o ./$out_prefix \
    -e $doc_pivot_id \
    --no-heuristic

  sort -k1,1n ${out_prefix}.fna.extract_dap.txt > ${out_prefix}_dap.txt
  """
}

/*
 * Run doc_pfp to extract document array profile.
 */
process EXTRACT_DAP_TMP {
  label 'extract_dap'
  debug true

  input:
  path doc_pfp
  path example_full_dap
  path example_fna
  val out_prefix

  output:
  path "full_dap.txt", emit: dap
  path "${out_prefix}.fna", emit: fna

  script:
  """
  cp $example_full_dap full_dap_cp.txt
  cp $example_fna ${out_prefix}.fna
  """
}


process INDEX_FNA {
  label 'index_fna'
  debug true

  input:
  path fna_file

  output:
  path "${fna_file}.fai"

  script:
  """
  samtools faidx $fna_file
  """
}


/*
 * Define the `index` process
 *
 * Note: odd bug where for index.sh input -b must be the first parameter
 */
process INDEX {
  label 'index_records'
  debug true

  input:
  path omem_index
  path dap_to_ms_py
  path fna_fai
  path dap_txt
  val records
  val out_prefix
  val num_docs

  output:
  path "${out_prefix}.bed.gz", emit: bed_gz
  path "${out_prefix}.bed.gz.tbi", emit: bed_gz_tbi

  script:
  """
  ./$omem_index \
    -b ${out_prefix} \
    -f $fna_fai \
    -d $dap_txt \
    -r $records \
    -t $num_docs
  """
}

