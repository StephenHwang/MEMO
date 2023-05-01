/*
 * pipeline input parameters
 */


/*
 * Define the `extract_region` process
 */
process EXTRACT_REGION {
  label 'extract_region'
  debug true

  input:
  path omem_extract
  path input_bed_gz
  path input_bed_gz_tbi
  val region

  output:
  path "omem_olaps_*.bed"

  script:
  """
  ./$omem_extract \
    -b $input_bed_gz \
    -r $region
  """
}


/*
 * Define the `query` process
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
    -r $region_bed
  """
}


