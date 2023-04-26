/*
 * pipeline input parameters
 */

// paths to software
// params.omem =  "/home/stephen/Documents/projects/langmead_lab/omem/src/omem"
params.dap_to_ms_py =  "/home/stephen/Documents/projects/langmead_lab/omem/src/dap_to_ms_bed.py"
params.omem_index =  "/home/stephen/Documents/projects/langmead_lab/omem/src/index.sh"
params.omem_query =  "/home/stephen/Documents/projects/langmead_lab/omem/src/query.sh"
params.doc_pfp = "/home/stephen/Documents/projects/langmead_lab/docprofiles/build_printout/pfp_doc64"

// example because currently skipping
params.example_full_dap = "/home/stephen/Documents/projects/langmead_lab/analysis/order_mems/bacteria_5/e_coli_pivot/full_dap.txt"
params.example_fna =      "/home/stephen/Documents/projects/langmead_lab/analysis/order_mems/bacteria_5/e_coli_pivot/input.fna"


// asdf
params.document_listing = "/home/stephen/Documents/projects/langmead_lab/omem/data/bacteria_pangeome_fasta/bacteria_pangenome_paths_with_annot.txt"
params.out_dir = "results"
params.out_prefix = "e_coli"

// parameters
params.pivot_idx = 1
params.num_docs = 4
params.K = 12
params.records = "NZ_CP015023.1 NZ_CP015022.1"
params.region = "NZ_CP015023.1:0-5506800"

/*
 * pipeline log
 */
log.info """\
  OMEM NF pipeline
  ===================================
  Output directory        : ${params.out_dir}
  Output prefix           : ${params.out_prefix}
  Document collection     : ${params.document_listing}
  Num non-pivot docs      : ${params.num_docs}
  Pivot doc idx           : ${params.pivot_idx}
  Index records           : ${params.records}
  Query region            : ${params.region}
  K                       : ${params.K}
  """
  .stripIndent()


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


/*
 * Define the `extract_region` process
 */
process EXTRACT_REGION {
  label 'extract_region'
  debug true

  input:
  path input_bed_gz
  path input_bed_gz_tbi
  val region

  output:
  path "region.bed"

  shell:
  '''
  QUERY_CHR="$(echo !{region} | cut -d':' -f1)"
  QUERY_START_END="$(echo !{region} | cut -d':' -f2)"
  QUERY_START="$(echo $QUERY_START_END | cut -d'-' -f1)"
  QUERY_END="$(echo $QUERY_START_END | cut -d'-' -f2)"

  echo -e "$QUERY_CHR\t$QUERY_START\t$QUERY_END" > query.bed
  tabix !{input_bed_gz} !{region} | \
    bedtools intersect -sorted -wa -f 1 -a 'stdin' -b query.bed \
    > region.bed
  '''
}


/*
 * Define the `index` process
 */
process QUERY {
  label 'query_region'
  debug true
  publishDir params.out_dir, mode:'copy'

  input:
  path omem_query
  val K
  val num_docs
  path region_bed

  output:
  path "omem_${K}mer.bed"

  script:
  """
  ./$omem_query \
    -k $K \
    -n $num_docs \
    -r $region_bed  \
    -p
  """
}

/*
 * Basic query workflow for casting `k-mers` in `chr:start-end`.
 */
workflow {

  // dap_ch = EXTRACT_DAP(params.doc_pfp,
                       // params.document_listing,
                       // params.out_prefix,
                       // params.pivot_idx)
  dap_ch = EXTRACT_DAP_TMP(params.doc_pfp,
                       params.example_full_dap,
                       params.example_fna,
                       params.out_prefix)

  index_fna_ch = INDEX_FNA(dap_ch.fna)

  index_ch = INDEX(params.omem_index,
                   params.dap_to_ms_py,
                   index_fna_ch,
                   dap_ch.dap,
                   params.records,
                   params.num_docs,
                   params.out_prefix)

  extract_region_ch = EXTRACT_REGION(index_ch.bed_gz,
                                     index_ch.bed_gz_tbi,
                                     params.region)

  query_ch = QUERY(params.omem_query,
                   params.K,
                   params.num_docs,
                   extract_region_ch)

}


workflow.onComplete {
    log.info ( workflow.success ? "\nWORKFLOW SUCCEEDED!" : "Oops .. something went wrong" )
}
