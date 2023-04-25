/*
 * pipeline input parameters
 */

// paths to software
params.omem = "$projectDir/../src/omem"
params.doc_pfp = "/home/stephen/Documents/projects/langmead_lab/docprofiles/build_printout/pfp_doc64"

// asdf
params.out_dir = "results"
params.out_prefix = "e_coli"

params.document_listing = "/scratch4/mschatz1/sjhwang/data/bacteria_5/paths/bacteria_pangenome_paths_with_annot.txt"

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
  path omem
  path fna_fai
  path dap_txt
  val records
  val num_docs
  path out_prefix

  output:
  path "${out_prefix}.gz"

  script:
  """
  $omem index \
    -f $fna_fai \
    -d $dap_txt \
    -r $records \
    -b $out_bed_name \
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
  val region

  output:
  path "region.bed"

  script:
  parse

  """
  QUERY_CHR=$(echo $region | cut -d':' -f1)
  QUERY_START_END=$(echo $region | cut -d':' -f2)
  QUERY_START=$(echo $QUERY_START_END | cut -d'-' -f1)
  QUERY_END=$(echo $QUERY_START_END | cut -d'-' -f2)

  echo -e "$QUERY_CHR\t$QUERY_START\t$QUERY_END" > query.bed
  tabix $input_bed_gz $region | \
    bedtools intersect -sorted -wa -f 1 -a 'stdin' -b query.bed \
    > region.bed
  """
}


/*
 * Define the `index` process
 */
process QUERY {
  label 'query_region'
  debug true
  publishDir params.out_dir, mode:'copy'

  input:
  path omem
  val K
  val num_docs
  path region_bed

  output:
  path "omem_${K}mer.bed"

  script:
  """
  $omem query \
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
  dap_ch = EXTRACT_DAP(params.doc_pfp,
                       params.document_listing,
                       params.out_prefix,
                       params.pivot_idx)

  index_fna_ch = INDEX_FNA(dap_ch.out.fna)

  index_ch = INDEX(params.omem,
                   index_fna_ch,
                   dap_ch.out.dap,
                   params.records,
                   params.num_docs,
                   params.out_prefix)

  extract_region_ch = EXTRACT_REGION(index_ch,
                                     params.region)

  query_ch = QUERY(params.omem,
                   params.K,
                   params.num_docs,
                   extract_region_ch)

}


workflow.onComplete {
    log.info ( workflow.success ? "\nDONE!" : "Oops .. something went wrong" )
}
