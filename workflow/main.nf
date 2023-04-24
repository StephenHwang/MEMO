/*
 * pipeline input parameters
 */
params.documents = "/scratch4/mschatz1/sjhwang/data/bacteria_5/paths/bacteria_pangenome_paths_with_annot.txt"
params.doc_idx = 1
params.num_docs = 4
params.doc_pfp = "/scratch4/mschatz1/sjhwang/src/docprofiles/build_printout/pfp_doc64"
params.outdir = "results"
params.src = "../src/"

params.fna_fai = "../data/example_dap/input.fna.fai"
params.dap_txt = "../data/example_dap/full_dap.txt"
params.records = "NZ_CP015023.1 NZ_CP015022.1"
params.out_bed_name = "example_dap"

params.k = 12
params.region = "NZ_CP015023.1:0-5506800"


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
  label 'index_region'
  debug true

  input:
  path fna_fai
  path dap_txt
  val records
  path outdir
  path out_bed_name
  path src

  script:
  """
  $src/omem index \
    -f $fna_fai \
    -d $dap_txt \
    -r $records \
    -o $outdir \
    -b $out_bed_name \
    -p
  """
}

/*
 * Define the `index` process
 */
process QUERY {
  label 'query_region'
  debug true

  input:
  val k
  val num_docs
  path input_bed_gz
  path outdir
  val region   // NZ_CP015023.1:0-5506800
  path src

  script:
  """
  $src/omem query \
    -k $k \
    -n $num_docs \
    -o $outdir \
    -b $input_bed_gz \
    -r $region  \
    -p
  """
}

/*
 * Define the `index` process
 */
workflow {
  dap_ch = EXTRACT_DAP(params.doc_pfp,
                       params.documents,
                       params.doc_idx)

  index_ch = INDEX(params.doc_pfp,
                   params.documents,
                   params.records,
                   params.outdir,
                   params.out_bed_name,
                   params.src)

  query_ch = QUERY(prams.k,
                   prams.num_docs,
                   prams.out_bed_name, //+ gz
                   prams.outdir,
                   prams.region,       // NZ_CP015023.1:0-5506800
                   prams.src)

}


workflow.onComplete {
    log.info ( workflow.success ? "\nDONE!" : "Oops .. something went wrong" )
}
