#!/bin/bash
#$1 = dir
#$2 = number of threads
#$3 = fastme flag

dir=$1

# apples
while read query; do
    echo "for query ${query}:"
    time run_apples.py -t ${dir}/${query}/backbone.tree -s ${dir}/${query}/ref.fa -q ${dir}/${query}/query.fa -T $2 -o ${dir}/${query}/apples.jplace
    echo
done < ${dir}/queries.txt

# for each of the placed trees in apples.jplace, compare to the true tree (rose.mt) and the backbone tree it was removed from (tree.nwk)
while read query; do
    # compare to backbone tree (tree.nwk)
    sed -E 's/[\s\S]*\"tree\": \"\s*([^\"]*)\"[\s\S]*/\1/gm;t;d' ${dir}/${query}/apples.jplace | sed -E 's/(\{[0-9]*\}|,$)//gm;t;d' &> ${dir}/${query}/apples.tree
    python3 treecompare.py ${dir}/tree.nwk ${dir}/${query}/apples.tree
done < ${dir}/queries.txt
