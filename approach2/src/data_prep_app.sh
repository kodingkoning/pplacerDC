#!/bin/bash

dir=$1

while read query; do
#    then
    faSomeRecords -exclude ${dir}/aln_shuf.fa <(echo "${query}") ${dir}/${query}/ref.fa
    faSomeRecords ${dir}/aln_shuf.fa <(echo "${query}") ${dir}/${query}/query.fa
#    fi

done < ${dir}/queries.txt
