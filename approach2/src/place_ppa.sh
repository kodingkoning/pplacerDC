#!/bin/bash
#$1 = dir
#$2 = number of threads

dir=$1

# APPLES-pplacer
raxml8 -f e -t ${dir}/tree.nwk -m GTRGAMMA -s ${dir}/rose.aln.true.fasta -n REF8 -p 1984 -T 8 -w ${dir}
while read query; do
    # python3 pplacerAPPLES.py -t ${dir}/${query}/backbone.tree -q ${query} -s ${dir}/RAxML_info.REF8  -r ${dir}/rose.aln.true.fasta -o ${dir}/${query}/pplacerAPPLES.tree -j 4 -m 500 -n 1
done < ${dir}/queries.txt
