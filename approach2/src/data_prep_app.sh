#!/bin/bash

dir=$1
newick_utils=/home/ekoning2/scratch/newick-utils-1.6/src

conda activate ppa
while read query; do
#    then
    mkdir -p ${dir}/${query}
    ${newick_utils}/nw_prune ${dir}/gtr-gamma-raxml.REF $query > ${dir}/${query}/backbone.tree
    faSomeRecords -exclude ${dir}/rose.aln.true.fasta <(echo "${query}") ${dir}/${query}/ref.fa
    faSomeRecords ${dir}/rose.aln.true.fasta <(echo "${query}") ${dir}/${query}/query.fa
#    fi

done < ${dir}/queries.txt
