#!/bin/bash

dir=$1

while read query; do
#    then
    mkdir -p ${dir}/${query}
    /home/erk24/Documents/CS581/finalProject/cs581-project/approach2/src/newick-utils-1.6/src/nw_prune ${dir}/gtr-gamma-raxml.REF $query > ${dir}/${query}/backbone.tree
    faSomeRecords -exclude ${dir}/rose.aln.true.fasta <(echo "${query}") ${dir}/${query}/ref.fa
    faSomeRecords ${dir}/rose.aln.true.fasta <(echo "${query}") ${dir}/${query}/query.fa
#    fi

done < ${dir}/queries.txt
