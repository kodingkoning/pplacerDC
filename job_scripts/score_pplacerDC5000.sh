#!/bin/bash

module load bwpy
export PATH=$PATH:/u/sciteam/ekoning2/scratch/apples:/u/sciteam/ekoning2/scratch/pplacer-Linux-v1.1.alpha19:/u/sciteam/ekoning2/scratch/cs581-project/approach2/src:/u/sciteam/ekoning2/scratch/cs581-project/approach1/src:/u/sciteam/ekoning2/scratch/cs581-project/common:/u/sciteam/ekoning2/scratch/newick-utils-1.6/src:/u/sciteam/ekoning2/scratch/raxml-ng/bin

maxTreeSize=2500

cd ~/scratch/RNASim-VS/variable-size/data

## for each dataset size: 500, 1000, 5000, 10 000, 50 000, 100 000, 200 000 (start with 500)

T5=`mktemp -t time_XXXXXX1.txt`
T6=`mktemp -t time_XXXXXX2.txt`

for dataset in 5000; do
  cd $dataset
  rm delta_error_approach1.txt

  echo "On datset $dataset" 
  for replicate in 0 1 2 3 4; do
    echo "On replicate $replicate" 
    head -n 200 $replicate/queries.txt > tmp-approach1-$replicate.txt
    while read query; do
      echo $query 
      if [ ! -f $replicate/$query/backbone_pp.tre ]; then
        nw_prune $replicate/RAxML_result.REF $query &> $replicate/$query/backbone_pp.tre
      fi
      if [ ! -f $replicate/$query/backbone_true.tre ]; then
        nw_prune $replicate/true_topo.tree $query &> $replicate/$query/backbone_true.tre
      fi
      python3.6 $(which treecompare.py) $replicate/true_topo.tree $replicate/$query/approach1.tre $replicate/$query/backbone_pp.tre $replicate/$query/backbone_true.tre >> delta_error_approach1.txt
    done < tmp-approach1-$replicate.txt
  done
  cd ../
done

