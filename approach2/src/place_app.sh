#!/bin/bash
#$1 = dir
#$2 = number of threads
#$3 = fastme flag

dir=$1

# apples
while read query; do
    run_apples.py -t ${dir}/${query}/backbone.tree -s ${dir}/${query}/ref.fa -q ${dir}/${query}/query.fa -T $2 -o ${dir}/${query}/apples.jplace
done < ${dir}/queries.txt
