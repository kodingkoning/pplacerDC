#!/bin/bash
#
# Required in path: python3, pplacer, epa-ng, newick-utils-1.6/src, global/src/shell

# evaluation scripts from: https://github.com/smirarab/global.git
# export WS_HOME=/location/of/global

maxTreeSize=2500

# cd to data directory
cd ../data

T5=`mktemp -t time_XXXXXX1.txt`

for dataset in 1000 5000 10000 50000 100000; do
  cd $dataset
  echo "On datset $dataset" 
  for replicate in 0 1 2 3 4; do
    echo "On replicate $replicate"
    while read query; do
      mkdir -p $replicate/$query
      echo $query 
      if [ ! -f $replicate/$query/backbone_app.tre ]; then
        nw_prune $replicate/true_me.fasttree $query &> $replicate/$query/backbone_app.tre
      fi
      if [ ! -f $replicate/$query/query.fa ]; then
        faSomeRecords.py --fasta $replicate/aln_dna.fa --records $query --outfile $replicate/$query/query.fa
      fi
      if [ ! -f $replicate/$query/aln_dna.fa ]; then
        faSomeRecords.py --exclude --fasta $replicate/aln_dna.fa --records $query --outfile $replicate/$query/aln_dna.fa
      fi
      /usr/bin/time -o $T5 -f "%e\t%M" run_apples.py -T 16 -t $replicate/$query/backbone_app.tre -q $replicate/$query/query.fa -s $replicate/$query/aln_dna.fa -o $replicate/$query/apples.jplace
      cat $T5 >> time_apples-$replicate.txt
      cat $T5 > $replicate/$query/apples_time.txt
      echo $query $(cat $T5)

      guppy tog -o $replicate/$query/apples.tre $replicate/$query/apples.jplace
      n1=`compareTrees.missingBranch $replicate/true_topo.tree <(nw_topology $replicate/${query}/apples.tre) -simplify | awk '{printf $2}'`
      n2=`compareTrees.missingBranch $replicate/true_topo.tree <(nw_topology $replicate/${query}/backbone_app.tre) -simplify | awk '{printf $2}'`

      python -c "print ($n1-$n2)" > $replicate/$query/apples_error.txt
      python -c "print ($n1-$n2)"

      # recommended when running many queries
      #rm $replicate/$query/backbone_app.tre
      #rm $replicate/$query/apples.tre
      #rm $replicate/$query/apples.jplace

    done < $replicate/queries.txt 
  done
  cd ../
done

