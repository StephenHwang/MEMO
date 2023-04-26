#!/usr/bin/env nextflow

// Import modules
include {EXTRACT_DAP_TMP}               from './modules/index.nf'
include {EXTRACT_DAP; INDEX_FNA; INDEX} from './modules/index.nf'
include {EXTRACT_REGION; QUERY}         from './modules/query.nf'

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
  """
  .stripIndent()


// ///////////////////         Pipeline workflows           ///////////////////

/*
 * Workflow for creating and indexing order MEMs
 */
workflow index {
  main:
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
  emit:
    index_ch.bed_gz
    index_ch.bed_gz_tbi
}

/*
 * Workflow for extracting and casting shadows for specified region of
 *   an order MEM index file.
 */
workflow query {
  take:
    bed_gz
    bed_gz_tbi
  main:
    extract_region_ch = EXTRACT_REGION(bed_gz,
                                       bed_gz_tbi,
                                       params.region)
    query_ch = QUERY(params.omem_query,
                     params.K,
                     params.num_docs,
                     extract_region_ch)
  emit:
    query_ch
}

/*
 * Workflow for indexing and finding k-shadows in region `chr:start-end`.
 */
workflow query_region_with_k {
  index_ch = index()
  query(index_ch)
}

/*
 * Workflow for finding k* for a region
 */
workflow find_k_star {
  // running multiple cast w K and then normalize arg max 
  // might need a python script
  // output is k* and omem-delta plot
}


workflow panagram_plot {
}



/*
 * run:
 *   nextflow run main.nf -process.echo -entry query_region_with_k
 */
workflow {
  query_region_with_k()
  find_k_star()
}

workflow.onComplete {
    log.info ( workflow.success ? "\nWORKFLOW SUCCEEDED!" : "Oops .. something went wrong" )
}
