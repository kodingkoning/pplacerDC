#!/bin/bash
#
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=secondary
#SBATCH --output=raxml200K0.o%j
#SBATCH --error=raxml200K0.e%j
#

module load python/3

export PATH=$PATH:/home/ekoning2/scratch/apples/apples-1.2.0:/home/ekoning2/scratch/pplacer-Linux-v1.1.alpha19:/home/ekoning2/scratch/pplacerDC/approach2/src:/home/ekoning2/scratch/pplacerDC/approach1/src:/home/ekoning2/scratch/pplacerDC/common:/home/ekoning2/scratch/newick-utils-1.6/src:/home/ekoning2/scratch/raxml-ng/bin:/home/ekoning2/scratch/epa-ng/bin:/home/ekoning2/scratch/raxml7/RAxML-7.2.6

maxTreeSize=2500

cd ~/scratch/RNASim-VS/variable-size/data

for dataset in 200000; do
  cd $dataset

  echo "On datset $dataset" 
  for replicate in 0; do
    echo "On replicate $replicate"
    cd $replicate
    raxmlHPC-PTHREADS -f e -t ./true_unweighted.fasttree -m GTRGAMMA -s ./aln_dna.phy -n REF7 -p 1984 -T 16
    cd ..
  done
  cd ../
done

