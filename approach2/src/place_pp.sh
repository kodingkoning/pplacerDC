#!/bin/bash
# raxml info file
#$1 directory to find all input files

# get info file from: raxml7

export LC_ALL=C

dir=$1
raxml_info=${dir}/RAxML_info.REF7

# while read query; do
#         # pplacer on aln_qry
#        pplacer -m GTR -s ${raxml_info} -t ${dir}/${query}/backbone.tree -o ${dir}/${query}/pp.jplace ${dir}/rose.aln.true.fasta -j 1
# #     ./pplacerAPPLES.py -t input.tre -q ${query} -s ${dir}/RAxML_info.REF7  -r ${dir}/rose.aln.true.fasta -o ${dir}/${query}/pplacerAPPLES.tree -j ${threads} -m 500 -n 1

# done < ${dir}/queries.txt

# for each of the placed trees (pplacerAPPLES.tree), compare to the true tree (rose.mt) and the backbone tree it was removed from (tree.nwk)
while read query; do
    # compare to backbone tree (tree.nwk)
#     sed -E 's/\{\"tree\":\s*\"\s*([^\"]*)\"[\s\S]*/\1/gm;t;d' ${dir}/${query}/pp.jplace | sed -E 's/(\{[0-9]*\}|,$)//gm;t;d' &> ${dir}/${query}/pp.tree
#     sed -E 's/.*\{\"tree\":\s*\"([^\"]*)\"[\s\S]*/\1/gm;t;d' ${dir}/${query}/pp.jplace &> ${dir}/${query}/pp.tree
    sed -E 's/(\{[0-9]*\}|,$)//gm;t;d' ${dir}/${query}/pp.jplace | sed -E 's/\s*\"([^\"]*;)\"[\s\S]*/\1/gm;t;d' &> ${dir}/${query}/pp.tree
    python3 treecompare.py ${dir}/tree.nwk ${dir}/${query}/pp.tree
done < ${dir}/queries.txt
