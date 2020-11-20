#!/bin/bash
#$1 = dir
#$2 = number of threads
#$3 = fastme flag

dir=$1

# apples
# while read query; do
#     run_apples.py -t ${dir}/${query}/backbone.tree -s ${dir}/${query}/ref.fa -q ${dir}/${query}/query.fa -T $2 -o ${dir}/${query}/apples.jplace
# done < ${dir}/queries.txt

# for each of the placed trees (pplacerAPPLES.tree), compare to the true tree (rose.mt) and the backbone tree it was removed from (tree.nwk)
while read query; do
    # compare to backbone tree (tree.nwk)
    # sed -E 's/[\s\S]*\"tree\": \"\s*([^\"]*)\"[\s\S]*/\1/gm;t;d' ${dir}/${query}/apples.jplace | sed -E 's/\{[0-9]*\}//gm;t;d' &> ${dir}/${query}/apples.tree
    sed -E 's/[\s\S]*\"tree\": \"\s*([^\"]*)\"[\s\S]*/\1/gm;t;d' ${dir}/${query}/apples.jplace | sed -E 's/(\{[0-9]*\}|,$)//gm;t;d' &> ${dir}/${query}/apples.tree
    # sed -E 's/[\s\S]*\"tree\": \"\s*([^\"]*)\"[\s\S]*/\1/gm;t;d' -E 's/\{[0-9]*\}//gm;t;d' ${dir}/${query}/apples.jplace &> ${dir}/${query}/apples.tree
    # sed -E 's/\{[0-9]*\}//gm;t;d'<(sed -E 's/[\s\S]*\"tree\": \"\s*([^\"]*)\"[\s\S]*/\1/gm;t;d' ${dir}/${query}/apples.jplace) &> ${dir}/${query}/apples.tree
    # sed -E 's/[\s\S]*\"tree\": \"\s*([^\"]*)\"[\s\S]*/\1/gm;t;d' ${dir}/${query}/apples.jplace &> ${dir}/${query}/apples.tree.temp
    # sed -E 's/\{[0-9]*\}//gm;t;d' ${dir}/${query}/apples.tree.temp &> ${dir}/${query}/apples.tree
    # rm ${dir}/${query}/apples.tree.temp
    python3 treecompare.py ${dir}/tree.nwk ${dir}/${query}/apples.tree
    # compare trees ${dir}/${query}/pplacerAPPLES.tree vs ${dir}/tree.nwk and ${dir}/rose.mt
done < ${dir}/queries.txt
