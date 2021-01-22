#!/bin/bash
#$1 dir
#$2 query size

# Modified from measure_app.sh in RNASim-VS Supplementary Materials of https://doi.org/10.1093/sysbio/syz063

export dir=$1
export size=$2

export R1=${dir}/results2_me_app.csv
export R2=${dir}/results2_pp.csv
export R3=${dir}/results2_epa.csv
export R4=${dir}/results2_ppdc.csv

echo -n "" > ${dir}/results2_me_app.csv
echo -n "" > ${dir}/results2_pp.csv
echo -n "" > ${dir}/results2_epa.csv
echo -n "" > ${dir}/results2_ppdc.csv

f() {
    query=$1
	# apples

	guppy tog -o $dir/$query/FM_MLSE.tree $dir/$query/apples.jplace

        n1=`compareTrees.missingBranch ${dir}/true_topo.tree <(nw_topology ${dir}/${query}/FM_MLSE.tree) -simplify | awk '{printf $2}'`
        n2=`compareTrees.missingBranch ${dir}/true_topo.tree <(nw_topology ${dir}/${query}/backbone_app.tre) -simplify | awk '{printf $2}'`
     
        python -c "print ($n1-$n2)" >> $R1

        # pplacer
	if (( $size < 5000 )); then
          guppy tog -o $dir/$query/pplacer.tree $dir/$query/pplacer.jplace

          n1=`compareTrees.missingBranch ${dir}/true_topo.tree <(nw_topology ${dir}/${query}/pplacer.tree) -simplify | awk '{printf $2}'`
          n2=`compareTrees.missingBranch ${dir}/true_topo.tree <(nw_topology ${dir}/${query}/backbone_pp.tre) -simplify | awk '{printf $2}'`
      
          python -c "print ($n1-$n2)" >> $R2
        fi

	# epa
	if (( $size < 50000 )); then
          guppy tog -o $dir/$query/tog_epa.tree $dir/$query/epa_result.jplace

          n1=`compareTrees.missingBranch ${dir}/true_topo.tree <(nw_topology ${dir}/${query}/tog_epa.tree) -simplify | awk '{printf $2}'`
          n2=`compareTrees.missingBranch ${dir}/true_topo.tree <(nw_topology ${dir}/${query}/backbone_epa.tre) -simplify | awk '{printf $2}'`
      
          python -c "print ($n1-$n2)" >> $R3
        fi

	# pplacerDC

        n1=`compareTrees.missingBranch ${dir}/true_topo.tree <(nw_topology ${dir}/${query}/approach1.tre) -simplify | awk '{printf $2}'`
        n2=`compareTrees.missingBranch ${dir}/true_topo.tree <(nw_topology ${dir}/${query}/backbone_pp.tre) -simplify | awk '{printf $2}'`
     
	python -c "print ($n1-$n2)" >> $R4
   }
  export -f f; xargs -P 16 --process-slot-var=index -I@ bash -c 'f @' < ${dir}/queries.txt

