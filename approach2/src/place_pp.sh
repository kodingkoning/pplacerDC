#!/bin/bash
# raxml info file
#$1 directory to find all input files

# get info file from: raxml-ng-mpi --msa rose.aln.true.fasta --prefix name --threads 10 --seed 12345 --model GTR+G --tree pars{1}

export LC_ALL=C

dir=$1
raxml_info=${dir}/gtr-gamma-raxml.info

while read query; do
        # pplacer on aln_qry
       pplacer -m GTR -s ${raxml_info} -t ${dir}/${query}/backbone.tree -o ${query}/pp.jplace ${dir}/rose.aln.true.fasta -j 1

done < ${dir}/queries.txt
