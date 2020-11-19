#$1 = dir
#$2 = number of threads
#$3 = raxml executable

# access RAxML v7.2.7 from: https://cme.h-its.org/exelixis/web/software/raxml/index.html
# extract files and then make prefered version (used Makefile.PHREADS.gcc here)

dir=$1
threads=$2
raxml=$3

# TODO: make script to do this for each of the 20 replicates

# APPLES-pplacer
if test ! -f ${dir}/RAxML_info.REF7; then
	# rm ${dir}/RAxML_info.REF7
    ${raxml} -f e -t ${dir}/tree.nwk -m GTRGAMMA -s ${dir}/rose.aln.true.fasta -n REF7 -p 1984 -T 8 -w ${dir}
fi
while read query; do
    /home/erk24/Documents/CS581/finalProject/cs581-project/approach2/src/newick-utils-1.6/src/nw_prune ${dir}/RAxML_result.REF7 ${query} &> ${query}-input.tre
    ./pplacerAPPLES.py -t ${query}-input.tre -q ${query} -s ${dir}/RAxML_info.REF8  -r ${dir}/rose.aln.true.fasta -o ${dir}/${query}/pplacerAPPLES.tree -j ${threads} -m 500 -n 1
    # ./pplacerAPPLES.py -t ${dir}/${query}/backbone.tree -q ${query} -s ${dir}/RAxML_info.REF8  -r ${dir}/rose.aln.true.fasta -o ${dir}/${query}/pplacerAPPLES.tree -j ${threads} -m 500 -n 1
done < ${dir}/queries.txt
