#!/bin/bash
#$1 = dir
#$2 = number of threads

dir=$1
threads=$2

#raxml=raxmlHPC-PTHREADS
prune=nw_prune
#raxml=/home/erk24/workspace/RAxML7/RAxML-7.2.7-ALPHA/raxmlHPC-PTHREADS #from https://cme.h-its.org/exelixis/web/software/raxml/index.html 
#prune=/home/ekoning2/scratch/newick-utils-1.6/src/nw_prune #from http://cegg.unige.ch/newick_utils
#
## APPLES-pplacer
#if test ! -f ${dir}/RAxML_info.REF7; then
#	# rm ${dir}/RAxML_info.REF7
#    ${raxml} -f e -t ${dir}/tree.nwk -m GTRGAMMA -s ${dir}/rose.aln.true.fasta -n REF7 -p 1984 -T 8 -w ${dir}
#fi

while read query; do
    echo "For query ${query}:"
    ${prune} ${dir}/RAxML_result.REF7 ${query} &> input.tre
    echo "scalepplacer.py -t input.tre -q ${query} -s ${dir}/RAxML_info.REF7  -r ${dir}/rose.aln.true.fasta -o ${dir}/${query}/scalepplacer.tree -j ${threads} -m 500"
    time scalepplacer.py -t input.tre -q ${query} -s ${dir}/RAxML_info.REF7  -r ${dir}/rose.aln.true.fasta -o ${dir}/${query}/scalepplacer.tree -j ${threads} -m 500
    echo
#done < ${dir}/queries.txt
done < ${dir}/1query.txt

# for each of the placed trees (scalepplacer.tree), compare to the true tree (rose.mt) and the backbone tree it was removed from (tree.nwk)
#while read query; do
#    # compare to backbone tree (tree.nwk)
#    python3 treecompare.py ${dir}/tree.nwk ${dir}/${query}/scalepplacer.tree
#    # compare trees ${dir}/${query}/scalepplacer.tree vs ${dir}/tree.nwk and ${dir}/rose.mt
#done < ${dir}/queries.txt

