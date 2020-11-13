#!/bin/bash

cd ~/Documents/CS581/finalProject/cs581-project/data/train/1000M1/1000M1/

for i in {0..0}; # for R0 through 19
do
    cd R${i}
    rm seqs.txt
    while read p; do
        if [[ $p == ">"* ]] ;
        then
            echo "${p:1}" >> seqs.txt
        fi
    done <rose.aln.true.fasta
    echo "" >> seqs.txt
    shuf -n 200 seqs.txt > queries.txt
    cd ..
done
