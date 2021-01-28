#!/bin/bash
#
# Required in path: python3, pplacer, epa-ng, newick-utils-1.6/src, global/src/shell

# evaluation scripts from: https://github.com/smirarab/global.git
# export WS_HOME=/location/of/global

maxTreeSize=2500

# cd to data directory
cd ../data

T5=`mktemp -t time_XXXXXX1.txt`

for dataset in 1000 5000 10000; do
  cd $dataset
  echo "On datset $dataset" 
  for replicate in 0 1 2 3 4; do
    echo "On replicate $replicate"
    while read query; do
      mkdir -p $replicate/$query
      echo $query 
      if [ ! -f $replicate/$query/backbone_epa.tre ]; then
        nw_prune $replicate/RAxML_result.REF8 $query &> $replicate/$query/backbone_epa.tre
      fi
      if [ ! -f $replicate/$query/query.fa ]; then
        faSomeRecords.py --fasta $replicate/aln_dna.fa --records $query --outfile $replicate/$query/query.fa
      fi
      /usr/bin/time -o $T5 -f "%e\t%M" epa-ng --ref-msa $replicate/aln_dna.fa --tree $replicate/$query/backbone_epa.tre --query $replicate/$query/query.fa --outdir $replicate/$query --model $replicate/RAxML_info.REF8 --redo -T 16
      echo $query $(cat $T5)
      cat $T5 >> time_epa-$replicate.txt
      cat $T5 > $replicate/$query/time_epa.txt

      guppy tog -o $replicate/$query/epa.tre $replicate/$query/epa_result.jplace
      n1=`compareTrees.missingBranch $replicate/true_topo.tree <(nw_topology $replicate/$query/epa.tre) -simplify | awk '{printf $2}'`
      n2=`compareTrees.missingBranch $replicate/true_topo.tree <(nw_topology $replicate/$query/backbone_epa.tre) -simplify | awk '{printf $2}'`
      python -c "print ($n1-$n2)" > $replicate/$query/epa_error.txt
      python -c "print ($n1-$n2)"

      # recommended when running many queries
      #rm $replicate/$query/epa_result.jplace
      #rm $replicate/$query/backbone_epa.tre
      #rm $replicate/$query/epa.tre
    done < $replicate/queries.txt 
  done
  cd ../
done

