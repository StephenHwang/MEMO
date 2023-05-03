/*
 * pipeline input parameters
 */


/*
 * Run doc_pfp to extract document array profile.
 */
process SUM_XS_PER_K {
  label 'sum_Xs_per_K'
  debug true

  input:
  path list_of_K_output_files
  val K_low
  val K_high
  val num_docs

  output:
  path "summary.txt"

  shell:
  '''
  FILE=summary.txt.tmp
  > $FILE
  for K in {!{K_low}..!{K_high}}
  do
    echo $K >> $FILE
    for ORDER_IDX in {1..!{num_docs}}
      do
        cat omem_${K}mer.bed | grep '${ORDER_IDX}\$' | awk '{print \$3 - \$2}' | paste -sd+ | bc >> $FILE
      done
  done

  # join every N+1 lines for each K
  cat $FILE | paste -d'\t' \$(printf -- '- %.0s' \$(seq 0 $NUM_ORDER_MEMS)) > summary.txt
  '''
}

