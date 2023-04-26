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


