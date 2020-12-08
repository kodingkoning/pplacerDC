#!/bin/bash

module load bwpy
export PATH=$PATH:/u/sciteam/ekoning2/scratch/apples:/u/sciteam/ekoning2/scratch/pplacer-Linux-v1.1.alpha19:/u/sciteam/ekoning2/scratch/cs581-project/approach2/src:/u/sciteam/ekoning2/scratch/cs581-project/approach1/src:/u/sciteam/ekoning2/scratch/cs581-project/common:/u/sciteam/ekoning2/scratch/newick-utils-1.6/src:/u/sciteam/ekoning2/scratch/raxml-ng/bin

maxTreeSize=2500

cd ~/scratch/RNASim-VS/variable-size/data

## for each dataset size: 500, 1000, 5000, 10 000, 50 000, 100 000, 200 000 (start with 500)

T5=`mktemp -t time_XXXXXX1.txt`
T6=`mktemp -t time_XXXXXX2.txt`

for dataset in 50000; do
  cd $dataset
  #rm delta_error*

  touch time_approach1.txt
  touch time_approach2.txt
  touch time_apples.txt
  touch time_pplacer.txt
  touch delta_error_approach1-2.txt
  touch delta_error_approach2-2.txt
  touch delta_error_apples-2.txt
  touch delta_error_pplacer.txt

  echo "On datset $dataset" >> run_output-score.log
  for replicate in 0; do
    echo "On replicate $replicate" >> run_output-score.log
    #tail -n 150 $replicate/queries.txt > tmp4.txt
    #head -n 20  tmp3.txt > tmp3.txt
    while read query; do
      echo $query >> run_output-score.log
      python3.6 $(which treecompare.py) $replicate/true_topo.tree $replicate/$query/approach1.tre $replicate/RAxML_result.REF >> delta_error_approach1-6.txt
      python3.6 $(which treecompare.py) $replicate/true_topo.tree $replicate/$query/approach2.tre $replicate/RAxML_result.REF >> delta_error_approach2-6.txt
      python3.6 $(which treecompare.py) $replicate/true_topo.tree $replicate/$query/approach2.tre.apples $replicate/RAxML_result.REF >> delta_error_apples-6.txt
      # TODO: get the pplacer times and results in cases where pplacer may not be able to run
    done < tmp6.txt
  done
  
  #grep -oP 'apples took \K.*(?= s)' run_output.log >> time_apples.txt
  if [ $maxTreeSize -gt $dataset ]; then
    #grep -oP 'Running pplacer... took \K.*(?= s)' run_output.log >> time_pplacer.txt
    cp delta_error_approach1.txt delta_error_pplacer.txt
  fi

  cd ../
done

