/*
 * pipeline input parameters
 */



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
  path "full_dap.txt", emit: dap
  path "${out_prefix}.fna", emit: fna

  script:
  """
  ./$doc_pfp build -r \
    -f $document_listing \
    -o $out_prefix \
    -e $doc_pivot_id
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
  val num_docs
  val out_prefix

  output:
  path "${out_prefix}.bed.gz", emit: bed_gz
  path "${out_prefix}.bed.gz.tbi", emit: bed_gz_tbi

  script:
  """
  ./$omem_index \
    -b ${out_prefix}.bed \
    -f $fna_fai \
    -d $dap_txt \
    -r $records \
    -t $num_docs \
    -p
  """
}

