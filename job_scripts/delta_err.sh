#!/bin/bash

T5=`mktemp -t time_XXXXXX1.txt`
T6=`mktemp -t time_XXXXXX2.txt`
#pathToScalePPlacer=$(which scalepplacer.py)

cd ../../variable-size/data/500 # TODO: remember to set this to the right path to the data
# If you use different data than the RNASim-VS, make sure the change the name of the raxml result/info, the alignment file, the input tree, and select the appropriate names of replicates

rm time_approach1.txt
rm delta_error_approach1.txt
rm run_output.log

touch time_approach1.txt
touch delta_error_approach1.txt
touch run_output.log
#for replicate in 0 1 2 3 4; do
#for replicate in 0 1 2; do
for replicate in 0; do
  echo "On replicate $replicate"
  cd $replicate
  #head -n 20 queries.txt > tmp.txt
  head -n 3 queries.txt > tmp.txt
  cat tmp.txt
  while read query; do
      nw_prune RAxML_result.REF $query &> input.tre
      #time pplacer -m GTR -s RAxML_info.REF -t input.tre -o pp.jplace aln_dna.fa -j 1
      /usr/bin/time -o $T5 -f "%e" scalepplacer.py -j 16 -m 500 -s RAxML_info.REF -t input.tre -q $query -r aln_dna.fa -o approach1.tre >> ../run_output.log
      /usr/bin/time -o $T6 -f "%e" pplacerAPPLES.py -n 1 -j 16 -m 500 -s RAxML_info.REF -t input.tre -q $query -r aln_dna.fa -o approach2.tre >> ../run_output.log
      cat $T5 >> ../time_approach1.txt
      cat $T6 >> ../time_approach2.txt
      python3 $(which treecompare.py) true_topo.tree approach1.tre RAxML_result.REF >> ../delta_error_approach1.txt
      python3 $(which treecompare.py) true_topo.tree approach2.tre RAxML_result.REF >> ../delta_error_approach2.txt
  #done < queries.txt
  done < tmp.txt
  cd ../
done

