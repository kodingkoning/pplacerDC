#!/bin/bash

cd ~/Documents/CS581/finalProject/cs581-project/data/train/1000M1/1000M1/

for i in {0..19}; # for R0 through 19
do
    cd R${i}
    ~/Documents/CS581/homework/03/FastTree -nosupport -gtr -gamma -nt -log tree.log < rose.aln.true.fasta > tree.nwk
    cd ..
done