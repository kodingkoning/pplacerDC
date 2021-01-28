#!/bin/bash

# Required in path: python3, pplacer, newick-utils-1.6/src, raxml-ng, global/src/shell

# evaluation scripts from: https://github.com/smirarab/global.git
# export WS_HOME=/location/of/global

maxTreeSize=2500

# cd to data directory
cd ../data

T5=`mktemp -t time_XXXXXX1.txt`

for dataset in 1000; do
  cd $dataset
  echo "On datset $dataset" 
  for replicate in 0 1 2 3 4; do
    echo "On replicate $replicate"
    while read query; do
      echo $query 
      mkdir -p $replicate/$query
      if [ ! -f $replicate/$query/backbone_pp.tre ]; then
        nw_prune $replicate/RAxML_result.REF $query &> $replicate/$query/backbone_pp.tre
      fi
      /usr/bin/time -o $T5 -f "%e\t%M" pplacer -m GTR -s $replicate/RAxML_info.REF -t $replicate/$query/backbone_pp.tre -o $replicate/$query/pplacer.jplace $replicate/aln_dna.fa -j 16
      cat $T5 >> time_pplacer-$replicate.txt
      cat $T5 > $replicate/$query/pplacer_time.txt

      guppy tog -o $replicate/$query/pplacer.tre $replicate/$query/pplacer.jplace
      n1=`compareTrees.missingBranch $replicate/true_topo.tree <(nw_topology $replicate/${query}/pplacer.tre) -simplify | awk '{printf $2}'`
      n2=`compareTrees.missingBranch $replicate/true_topo.tree <(nw_topology $replicate/${query}/backbone_pp.tre) -simplify | awk '{printf $2}'`

      python -c "print ($n1-$n2)" > $replicate/$query/pplacer_error.txt
      python -c "print ($n1-$n2)"

      # recommended when running many queries
      #rm $replicate/$query/pplacer.tre
      #rm $replicate/$query/pplacer.jplace
      #rm $replicate/$query/backbone_pp.tre 

      echo $query $(cat $T5)
    done < $replicate/queries.txt 
  done
  cd ../
done

