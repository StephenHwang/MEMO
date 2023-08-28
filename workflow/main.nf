#!/usr/bin/env nextflow

// Import modules
include {PROCESS_FASTA; DAP_PREPARE; MONI_MS; MS_TO_DAP}        from './modules/index.nf'
include {DOC_PFP_EXTRACT_DAP}                                           from './modules/index.nf'
include {INDEX_FNA; INDEX}                                      from './modules/index.nf'
include {EXTRACT_REGION; QUERY}                                 from './modules/query.nf'
include {SUM_XS_PER_K; FIND_K_STAR}                             from './modules/analyze.nf'

// Pipeline log
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
  K low                   : ${params.K_low}
  K high                  : ${params.K_high}
  """
  .stripIndent()


//////////////////////         Pipeline workflows           ///////////////////

workflow moni_ms_to_dap {
  main:
  // Channel of files and doc index
  file_list_ch = Channel.fromPath(params.document_listing).splitCsv(sep:" ")

  // Clean and rc fastas
  processed_fa_ch = PROCESS_FASTA(file_list_ch,
                                  params.preprocess_moni_fasta)

  // Find fasta headers of pivot and concatenate processed query fastas
  dap_prep_ch = DAP_PREPARE(processed_fa_ch.processed_fasta.collect(),
                            params.document_listing,
                            params.pivot_idx)

  // Find MS lengths for each fasta in the collection
  ms_ch = MONI_MS(processed_fa_ch.processed_rc_fasta.flatten(),
                  dap_prep_ch.query_fasta,
                  dap_prep_ch.ref_records,
                  params.moni,
                  params.moni_deps,
                  params.moni_src)

  // Make the full DAP
  dap_ch = MS_TO_DAP(ms_ch.collect(),
                     params.document_listing)
}


/*
 * Workflow for creating and indexing order MEMs
 */
workflow index {
  main:
    dap_ch = DOC_PFP_EXTRACT_DAP(params.doc_pfp,
                                 params.document_listing,
                                 params.out_prefix,
                                 params.pivot_idx)
    index_fna_ch = INDEX_FNA(dap_ch.fna)
    index_ch = INDEX(params.omem_index,
                     params.dap_to_ms_py,
                     index_fna_ch,
                     dap_ch.dap,
                     params.records,
                     params.out_prefix,
                     params.num_docs)
  emit:
    index_ch.bed_gz
    index_ch.bed_gz_tbi
}


/*
 * Workflow for extracting region from overlap order MEM index
 */
workflow extract {
  take:
    bed_gz
    bed_gz_tbi
  main:
    extract_region_ch = EXTRACT_REGION(params.omem_extract,
                                       bed_gz,
                                       bed_gz_tbi,
                                       params.region)
  emit:
    extract_region_ch
}


/*
 * Workflow for casting shadows over overlap order MEM bed file
 */
workflow query {
  take:
    extract_region_ch
  main:
    query_ch = QUERY(params.omem_query,
                     params.K,
                     params.num_docs,
                     extract_region_ch)
  emit:
    query_ch
}

/*
 * Workflow for indexing and finding k-shadows in region `chr:start-end`
 */
workflow query_region_with_k {
  index_ch = index()
  extract_ch = extract(index_ch)
  query(extract_ch)
}


/*
 * Workflow for varying casted K in region `chr:start-end`
 */
workflow vary_k {
  main:
    index_ch = index()
    extract_ch = extract(index_ch)

    // vary K
    k_channel = Channel.of(params.K_low..params.K_high)
    query_ch = QUERY(params.omem_query,
                     k_channel,
                     params.num_docs,
                     extract_ch)

    // sum Xs
    sum_Xs_ch = SUM_XS_PER_K(query_ch.collect(),
                           params.K_low,
                           params.K_high,
                           params.num_docs)

  emit:
    sum_Xs_ch
}



workflow find_k_star {
  vary_k_ch = vary_k()
  find_k_val = FIND_K_STAR(params.find_k_star_py,
                           vary_k_ch,
                           params.num_docs)
  find_k_val.view()
}




workflow panagram_plot {
}


/*
 * run:
 *   nextflow run main.nf -process.echo -entry <workflow> -resume
 *    ie. query_region_with_k
 */
workflow {
  moni_ms_to_dap()
  query_region_with_k()
  vary_k()
  find_k_star()
  panagram_plot()
}


workflow.onComplete {
    log.info ( workflow.success ? "\nWORKFLOW SUCCEEDED!" : "Oops .. something went wrong" )
}
