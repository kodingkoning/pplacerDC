#!/bin/bash
#
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=secondary
#SBATCH --output=one/pplacerDC100000.o%j
#SBATCH --error=one/pplacerDC100000.e%j
#

module unload python/2
module load python/3

export PATH=/usr/local/python/3.7.0/bin:$PATH:/home/ekoning2/scratch/apples:/home/ekoning2/scratch/pplacer-Linux-v1.1.alpha19:/home/ekoning2/scratch/pplacerDC/approach2/src:/home/ekoning2/scratch/pplacerDC/approach1/src:/home/ekoning2/scratch/pplacerDC/common:/home/ekoning2/scratch/newick-utils-1.6/src:/home/ekoning2/scratch/raxml-ng/bin:/home/ekoning2/scratch/epa-ng/bin

maxTreeSize=2500

cd ~/scratch/RNASim-VS/variable-size/data

T5=`mktemp -t time_XXXXXX1.txt`

dataset=200000
replicate=0
query=L103

echo $dataset $replicate $query

cd $dataset

mkdir -p $replicate/$query

if [ ! -f $replicate/$query/backbone_pp.tre ]; then
  nw_prune $replicate/RAxML_result.REF $query &> $replicate/$query/backbone_pp.tre
fi
/usr/bin/time -o $T5 -f "%e\t%M" scalepplacer.py -j 16 -m $maxTreeSize -s $replicate/RAxML_info.REF -t $replicate/$query/backbone_pp.tre -q $query -r $replicate/aln_dna.fa -o $replicate/$query/approach1.tre 
cat $T5 >> time_pplacerDC-$replicate.txt
echo $query $(cat $T5) 

