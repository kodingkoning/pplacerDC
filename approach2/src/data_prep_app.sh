#!/bin/bash

dir=$1
newick_utils=/home/ekoning2/scratch/newick-utils-1.6/src


while read query; do
#    then
    mkdir -p ${dir}/${query}
    ${newick_utils}/nw_prune ${dir}/RAxML_result.REF7 $query > ${dir}/${query}/backbone.tree
    #${newick_utils}/nw_prune ${dir}/gtr-gamma-raxml.REF $query > ${dir}/${query}/backbone.tree
    # query.fa and ref.fa are required for running APPLES, so generate before running APPLES
    #faSomeRecords -exclude ${dir}/rose.aln.true.fasta <(echo "${query}") ${dir}/${query}/ref.fa
    #faSomeRecords ${dir}/rose.aln.true.fasta <(echo "${query}") ${dir}/${query}/query.fa
#    fi

done < ${dir}/queries.txt
#done < ${dir}/1query.txt
