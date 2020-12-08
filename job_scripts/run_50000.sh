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
  #rm -f time_approach1.txt
  #rm -f time_approach2.txt
  #rm -f delta_error_approach1.txt
  #rm -f delta_error_approach2.txt
  #rm -f time_apples.txt
  #rm -f time_pplacer.txt
  #rm -f delta_error_apples.txt
  #rm -f delta_error_pplacer.txt
  #rm -f run_output.log

  touch time_approach1.txt
  touch time_approach2.txt
  touch time_apples.txt
  touch time_pplacer.txt
  touch delta_error_approach1.txt
  touch delta_error_approach2.txt
  touch delta_error_apples.txt
  touch delta_error_pplacer.txt
  touch run_output.log

  echo "On datset $dataset" >> run_output.log
  for replicate in 0; do
    echo "On replicate $replicate" >> run_output.log
    #head -n 200 $replicate/queries.txt > tmp.txt
    while read query; do
      echo $query >> run_output.log
      nw_prune $replicate/RAxML_result.REF $query &> input.tre
      mkdir -p $replicate/$query
      /usr/bin/time -o $T5 -f "%e" scalepplacer.py -j 16 -m $maxTreeSize -s $replicate/RAxML_info.REF -t input.tre -q $query -r $replicate/aln_dna.fa -o $replicate/$query/approach1.tre >> run_output.log
      /usr/bin/time -o $T6 -f "%e" pplacerAPPLES.py -j 16 -m $maxTreeSize -n 1 -s $replicate/RAxML_info.REF -t input.tre -q $query -r $replicate/aln_dna.fa -o $replicate/$query/approach2.tre >> run_output.log
      cat $T5 >> time_approach1.txt
      cat $T6 >> time_approach2.txt

      #python3.6 $(which treecompare.py) $replicate/true_topo.tree $replicate/$query/approach1.tre $replicate/RAxML_result.REF >> delta_error_approach1.txt
      #python3.6 $(which treecompare.py) $replicate/true_topo.tree $replicate/$query/approach2.tre $replicate/RAxML_result.REF >> delta_error_approach2.txt
      #python3.6 $(which treecompare.py) $replicate/true_topo.tree $replicate/$query/approach2.tre.apples $replicate/RAxML_result.REF >> delta_error_apples.txt
      # TODO: get the pplacer times and results in cases where pplacer may not be able to run
    done < tmp.txt
  done
  
  grep -oP 'apples took \K.*(?= s)' run_output.log >> time_apples.txt
  if [ $maxTreeSize -gt $dataset ]; then
    grep -oP 'Running pplacer... took \K.*(?= s)' run_output.log >> time_pplacer.txt
    #cp delta_error_approach1.txt delta_error_pplacer.txt
  fi

  cd ../
done

