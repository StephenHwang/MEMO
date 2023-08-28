/*
 * Module for creating order-MEM overlap index
 */


/*
 * Pre-processed fastas
 */
process GATHER_FASTAS {
  label 'gather_fastas'
  debug true

  input:
  path document_listing

  output:
  path "*.fa", emit: raw_fastas

  script:
  """
  for line in \$(cut -f1 -d' ' $document_listing)
  do
    cp \$line .
  done
  """
}


/*
 * Pre-processed fastas 2
 */
process GATHER_FASTAS_2 {
  label 'gather_fastas_2'
  debug true

  input:
  path document_listing

  output:
  path "*.fa", emit: raw_fastas

  script:
  """
  for line in \$(cut -f1 -d' ' $document_listing)
  do
    cp \$line .
  done
  """
}






/*
 * Pre-processed fastas
 */
process PROCESS_FASTA {
  label 'process_fastas'
  debug true
  conda '/home/shwang45/miniconda3/envs/python'

  input:
  path fasta
  path preprocess_moni_fasta

  output:
  path "*_clean.fa", emit: processed_fasta

  script:
  """
  header="\$(basename ${fasta} .fa)"
  echo "\$header"

  # preprocess to base fa
  ./${preprocess_moni_fasta} \
    < ${fasta} \
    > "\${header}_clean.fa"
  """
}


/*
 * Pre-processed fastas
 */
process PROCESS_FASTA_RC {
  label 'process_rc_fastas'
  debug true
  conda '/home/shwang45/miniconda3/envs/python'

  input:
  path fasta
  path preprocess_moni_fasta

  output:
  path "*_clean_rc.fa", emit: processed_rc_fasta

  script:
  """
  header="\$(basename ${fasta} .fa)"
  echo "\$header"

  # preprocess to base fa
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

  # make query fasta
  # TODO: make in order of document_listing!!!
  cat *_clean.fa > query.fa
  """
}


/*
 * Find matching statistics with MONI
 */
process MONI_MS {
  label 'moni_ms'
  debug true

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
  echo "\$header"

  #  -r "\${header}_with_rc.fa" -f \
  # build moni index
  ./${moni} build \
    -r "\${header}.fa" -f \
    -o "index/\${header}"

  ls 
  echo "  ^ post-index" 

  # find MS
  ./${moni} ms \
    -i "index/\${header}" \
    -p ${query_fasta} \
    -o "\${header}"

  echo "Done \$header; start seqtk"
  ls

  seqtk subseq "\${header}.lengths" $ref_records | grep -v '^>' | tr ' ' '\n' | grep . > "\${header}_subset.lengths"
  """
}


/*
 * Run doc_pfp to extract document array profile.
 */
process MS_TO_DAP {
  label 'ms_to_dap'
  debug true

  input:
  path "*"
  path document_listing

  output:
  path "full_dap.txt", emit: full_dap

  script:
  """
  # in order of document_listings
  paste -d ' ' \$(ls *_subset.lengths) | nl -v0 -w1 -s' ' > full_dap.txt
  """
}

////////////////////////////////////////////////////////////////////////////////

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

