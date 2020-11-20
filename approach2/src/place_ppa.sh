#$1 = dir
#$2 = number of threads
#$3 = raxml executable

# access RAxML v7.2.7 from: https://cme.h-its.org/exelixis/web/software/raxml/index.html
# extract files and then make prefered version (used Makefile.PHREADS.gcc here)

dir=$1
threads=$2
raxml=$3

# APPLES-pplacer
# if test ! -f ${dir}/RAxML_info.REF7; then
# 	# rm ${dir}/RAxML_info.REF7
#     python3 fasta2phylip.py -r -i ${dir}/rose.aln.true.fasta -o ${dir}/truealignment.phylip
#     ${raxml} -f e -t ${dir}/tree.nwk -m GTRGAMMA -s ${dir}/truealignment.phylip -n REF7 -p 1984 -T 8 -w ${dir}
# fi
while read query; do
    echo
    echo "For query ${query}:"
    /home/erk24/Documents/CS581/finalProject/cs581-project/approach2/src/newick-utils-1.6/src/nw_prune ${dir}/RAxML_result.REF7 ${query} &> input.tre
    ./pplacerAPPLES.py -t input.tre -q ${query} -s ${dir}/RAxML_info.REF7  -r ${dir}/rose.aln.true.fasta -o ${dir}/${query}/pplacerAPPLES.tree -j ${threads} -m 500 -n 1
    # ./pplacerAPPLES.py -t ${dir}/${query}/backbone.tree -q ${query} -s ${dir}/RAxML_info.REF8  -r ${dir}/rose.aln.true.fasta -o ${dir}/${query}/pplacerAPPLES.tree -j ${threads} -m 500 -n 1
done < ${dir}/queries.txt

# for each of the placed trees (pplacerAPPLES.tree), compare to the true tree (rose.mt) and the backbone tree it was removed from (tree.nwk)
while read query; do
    # compare to backbone tree (tree.nwk)
    python3 treecompare.py ${dir}/tree.nwk ${dir}/${query}/pplacerAPPLES.tree
    # compare trees ${dir}/${query}/pplacerAPPLES.tree vs ${dir}/tree.nwk and ${dir}/rose.mt
done < ${dir}/queries.txt

