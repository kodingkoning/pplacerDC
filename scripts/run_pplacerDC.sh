#!/bin/bash
#
# Required in path: python3, pplacer, pplacerDC/approach1/src, pplacerDC/common, newick-utils-1.6/src, raxml-ng, global/src/shell

# evaluation scripts from: https://github.com/smirarab/global.git
# export WS_HOME=/location/of/global

maxTreeSize=2500

# cd to data directory
cd ../data

T5=`mktemp -t time_XXXXXX1.txt`

for dataset in 100 5000 10000 50000 100000; do
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
      /usr/bin/time -o $T5 -f "%e\t%M" pplacerDC.py -j 16 -m $maxTreeSize -s $replicate/RAxML_info.REF -t $replicate/$query/backbone_pp.tre -q $query -r $replicate/aln_dna.fa -o $replicate/$query/pplacerDC.tre 
      cat $T5 >> time_pplacerDC-$replicate.txt
      cat $T5 >> $replicate/$query/pplacerDC_time.txt
      echo $query $(cat $T5) 

      n1=`compareTrees.missingBranch $replicate/true_topo.tree <(nw_topology $replicate/${query}/pplacerDC.tre) -simplify | awk '{printf $2}'`
      n2=`compareTrees.missingBranch $replicate/true_topo.tree <(nw_topology $replicate/${query}/backbone_pp.tre) -simplify | awk '{printf $2}'`

      python -c "print ($n1-$n2)" > $replicate/$query/pplacerDC_error.txt
      python -c "print ($n1-$n2)"

      # recommended when running many queries
      #rm $replicate/$query/pplacerDC.tre
      #rm $replicate/$query/backbone_pp.tre 
    done < $replicate/queries.txt 
  done
  cd ../
done

