/*
 * pipeline input parameters
 */
params.documents = "/scratch4/mschatz1/sjhwang/data/bacteria_5/paths/bacteria_pangenome_paths_with_annot.txt"
params.doc_idx = 1
params.doc_pfp = "/scratch4/mschatz1/sjhwang/src/docprofiles/build_printout/pfp_doc64"
params.outdir = "results"
params.src = "../src/"

params.fna_fai = "../data/example_dap/input.fna.fai"
params.dap_txt = "../data/example_dap/full_dap.txt"
params.records = "NZ_CP015023.1 NZ_CP015022.1" 
params.out_bed_name = "example_dap"


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


/*
 * Define the `index` process 
 */
process INDEX {
  input:
  path fna_fai
  path dap_txt
  path records
  path out_dir 
  path out_bed_name
  path src

  script:
  """
  $src/omem index \
    -f $fna_fai \
    -d $dap_txt \
    -r $records \
    -o $out_dir \
    -b $out_bed_name \
    -p
  """
}

















workflow {
  dap_ch = EXTRACT_DAP(params.doc_pfp,
                       params.documents,
                       params.doc_idx)

  index_ch = INDEX(params.doc_pfp,
                   params.documents,
                   params.doc_idx,
                   params.doc_idx,
                   params.doc_idx,
                   params.doc_idx,

  path fna_fai
  path dap_txt
  path records
  path out_dir 
  path out_bed_name
  path src



  // flow the channels

}

workflow.onComplete {
    log.info ( workflow.success ? "\nDONE!" : "Oops .. something went wrong" )
}
