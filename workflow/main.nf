/*
 * pipeline input parameters
 */
params.documents = "/home/stephen/Documents/projects/langmead_lab/omem/data/bacteria_pangeome_fasta/bacteria_pangenome_paths.txt"
params.doc_pfp = "/home/stephen/Documents/projects/langmead_lab/docprofiles/build/pfp_doc64"
params.outdir = "results"

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
 * Run doc_pfp to extract doc_pfp
 *   TODO: has to use the printout branch of docprofiles
 */
process EXTRACT_DAP {
  label 'extract_dap'
  debug true

  input:
  path documents
  path doc_pivot_id

  output:
  path 'dap'

  script:
  """
  params.doc_pfp \
    -f $documents \
    -o 'dap' \
    -e $doc_pivot_id
  """
}


/*
 * define the `index` process 
 * 
 */
process INDEX {
  input:
  path tmp

  output:
  path tmp

  script:
  """
  ./omem index \
    -f ../data/example_dap/input.fna.fai \
    -d ../data/example_dap/full_dap.txt \
    -r "NZ_CP015023.1 NZ_CP015022.1" \
    -o ../data/example_dap \
    -b omem.bed \
    -p
  """
}


/*
 * define the `query` process
 * 
 */
process QUERY {
  input:
  path tmp

  output:
  path tmp

  script:
  """
  ./omem query \
    -k 12 \
    -n 4 \
    -o ../data/example_dap \
    -b omem.bed.gz \
    -r NZ_CP015023.1:0-5506800 \
    -i \
    -p
  """
}


workflow {
  TODO
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDONE!" : "Oops .. something went wrong" )
}
