/*
 * pipeline input parameters
 */
params.documents = "/scratch4/mschatz1/sjhwang/data/bacteria_5/paths/bacteria_pangenome_paths_with_annot.txt"
params.doc_idx = 1
params.doc_pfp = "/scratch4/mschatz1/sjhwang/src/docprofiles/build_printout/pfp_doc64"
params.outdir = "results"
params.src = "../src/"

/*
 * pipeline log
 */
log.info """\
  OMEM NF pipeline
  ===================================
  documents     : ${params.documents}
  outdir        : ${params.outdir}
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
  path documents
  val doc_pivot_id

  script:
  """
  ./$doc_pfp build \
    -f $documents \
    -o './dap' \
    -e $doc_pivot_id
  """
}

workflow {
  dap_ch = EXTRACT_DAP(params.doc_pfp,
                       params.documents,
                       params.doc_idx)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDONE!" : "Oops .. something went wrong" )
}
