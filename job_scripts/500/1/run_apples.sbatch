#!/bin/bash
#
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=secondary
#SBATCH --output=apples500.o%j
#SBATCH --error=apples500.e%j
#

module unload python/2
module load python/3

export PATH=/usr/local/python/3.7.0/bin:$PATH:/home/ekoning2/scratch/apples:/home/ekoning2/scratch/pplacer-Linux-v1.1.alpha19:/home/ekoning2/scratch/pplacerDC/approach2/src:/home/ekoning2/scratch/pplacerDC/approach1/src:/home/ekoning2/scratch/pplacerDC/common:/home/ekoning2/scratch/newick-utils-1.6/src:/home/ekoning2/scratch/raxml-ng/bin:/home/ekoning2/scratch/epa-ng/bin

maxTreeSize=2500

cd ~/scratch/RNASim-VS/variable-size/data

T5=`mktemp -t time_XXXXXX1.txt`
T6=`mktemp -t time_XXXXXX2.txt`

for dataset in 500; do
  cd $dataset

  echo "On datset $dataset" 
  for replicate in 1; do
    echo "On replicate $replicate"
    rm -f time_apples-$replicate.txt 
    head -n 200 $replicate/queries.txt > tmp-apples-$replicate.txt
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
      echo $query $(cat $T5)
    done < tmp-apples-$replicate.txt
    rm tmp-apples-$replicate.txt
  done
  cd ../
done

