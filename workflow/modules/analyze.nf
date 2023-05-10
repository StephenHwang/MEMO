/*
 * Module for analyzing results
 */


/*
 * Cast range of K-shadows and summarize the number of Xs drawn per order
 */
process SUM_XS_PER_K {
  label 'sum_Xs_per_K'
  debug true
  publishDir params.out_dir, mode:'copy'

  input:
  path "*"
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
    echo -n $K >> $FILE
    for ORDER_IDX in {1..!{num_docs}}
      do
        echo -ne "\t" >> $FILE
        cat omem_${K}mer.bed | grep -w ${ORDER_IDX}\$ | awk '{print \$3 - \$2}' | paste -sd+ | bc | tr -d "\n" >> $FILE
      done
    echo "" >> $FILE
  done

  awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if(\$i ~ /^ *\$/) \$i = 0 }; 1' $FILE > summary.txt
  '''
}


/*
 * Cast range of K-shadows and summarize the number of Xs drawn per order
 */
process FIND_K_STAR() {
  label 'find_k_star'
  debug true
  publishDir params.out_dir, mode:'copy'

  input:
  path find_k_star_py 
  path "summary.txt"
  val num_docs

  output:
  path "k_star.txt"

  script:
  """
  ./$find_k_star_py \
    -f summary.txt \
    -d $num_docs \
    > k_star.txt
  """
}

