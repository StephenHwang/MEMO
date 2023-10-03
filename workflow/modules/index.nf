/*
 * Module for creating order-MEM overlap index
 */


/*
 * Pre-processed fastas
 */
process PROCESS_FASTA {
  label 'process_fastas'
  debug true
  executor 'local'
  conda '/home/shwang45/miniconda3/envs/python'

  input:
  tuple path(fasta), val(index)
  path preprocess_moni_fasta
  path document_listing
  val  doc_pivot_id

  output:
  path "*_clean_rc.fa", emit: processed_rc_fasta
  path "*_clean.fa", emit: processed_fasta, optional: true

  script:
  """
  header="\$(basename ${fasta} .fa)"
  pivot_file=\$(basename \$(sed '${doc_pivot_id}q;d' ${document_listing} | cut -f1 -d' ') .fa)

  # preprocess to base fa if is correct pivot
  if [ "\$header" = "\$pivot_file" ]
  then
    ./${preprocess_moni_fasta} \
      < ${fasta} \
      > "\${header}_clean.fa"
  fi

  # preprocess to clean and rc fasta
  ./${preprocess_moni_fasta} -r \
    < ${fasta} \
    > "\${header}_clean_rc.fa"
  """
}


/*
 * Concat records
 */
process DAP_PREPARE {
  label 'dap_prepare'
  executor 'local'
  debug true

  input:
  path "*"
  path document_listing
  val  doc_pivot_id

  output:
  path "query.fa", emit: query_fasta
  path "ref_records.lst", emit: ref_records

  script:
  """
  pivot_headers=\$(cat \$(cat $document_listing | grep "$doc_pivot_id\$" | cut -f1 -d' ') | grep "^>" | cut -f1 -d' ' | cut -f2 -d'>')
  echo "\$pivot_headers" | tr ' ' \\n > ref_records.lst
  cat *_clean.fa > query.fa
  """
}


/*
 * Find matching statistics with MONI
 */
process MONI_MS {
  label 'moni_ms'
  debug true
  executor 'local'
  // executor 'slurm'
  // queue 'defq'
  // memory '100 GB'
  // cpus 1
  // time '24h'

  input:
  path index_fasta
  path query_fasta
  path ref_records
  path moni
  path moni_deps
  path moni_src

  output:
  path "*_subset.lengths", emit: moni_length_files

  script:
  """
  mkdir index
  header="\$(basename $index_fasta .fa)"
  echo "MONI MS on: \$header"

  # build moni index
  # -t 8
  echo "  Building index"
  ./${moni} build \
    -r "\${header}.fa" -f \
    -o "index/\${header}"

  # find MS
  # -t 8
  echo "  Finding matching statistics"
  ./${moni} ms \
    -i "index/\${header}" \
    -p ${query_fasta} \
    -o "\${header}"

  seqtk subseq "\${header}.lengths" $ref_records | grep -v '^>' | tr ' ' '\n' | grep . > "\${header}_subset.lengths"
  """
}


/*
 * Run doc_pfp to extract document array profile.
 */
process MS_TO_DAP {
  label 'ms_to_dap'
  executor 'local'
  debug true

  input:
  path "*"
  path document_listing

  output:
  path "full_dap.txt", emit: full_dap

  script:
  """
  # make in order of doc_listing file
  # cut first column of document_listing, replace / delimeters with space, select names and then replace to align file suffix
  paste -d ' ' \$(cat ${document_listing} | cut -d' ' -f1 | tr '/' ' ' | awk '{print \$NF}' | sed 's/.fa/_clean_rc_subset.lengths/') | nl -v0 -w1 -s' ' > full_dap.txt
  """
}




////////////////////////////////////////////////////////////////////////////////

/*
 * Run doc_pfp to extract document array profile.
 */
process DOC_PFP_EXTRACT_DAP {
  label 'doc_pfp_extract_dap'
  executor 'local'
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


process INDEX_FNA {
  label 'index_fna'
  executor 'local'
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
  executor 'local'
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

