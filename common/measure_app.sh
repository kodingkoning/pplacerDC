#!/bin/bash
#$1 dir
#$2 query size

export dir=$1

# gamma variable
#gamma="$(grep -F "alpha[0]" ${1} | awk '{print $2}')"


export R1=${dir}/results_me_app.csv # default APPLES
export R2=${dir}/results_pp.csv # default pplacer
export R3=${dir}/results_ppa.csv # APPLES with pplacer
export R4=${dir}/results_pdc.csv # pplacer divide and conquer

echo -n "" > ${dir}/results_me_app.csv
echo -n "" > ${dir}/results_pp.csv
echo -n "" > ${dir}/results_ppa.csv
echo -n "" > ${dir}/results_pdc.csv

f() {
    query=$1
        # apples
        #(head -n 1 ${dir}/dist.mat; grep -P "^${query}\s+" ${dir}/dist.mat) > $dir/$query/vec.pd
        #grep -P "${query}\s+" ${dir}/dist.mat > $dir/$query/vec.pd
        #python ~/apples/run_apples.py -t ${dir}/${query}/backbone_app.tree -s ${dir}/${query}/ref.fa -q ${dir}/${query}/query.fa -T 8  -o ${dir}/${query}/apples.jplace

        guppy tog -o $dir/$query/apples.tree $dir/$query/apples.jplace

        n1=`compareTrees.missingBranch ${dir}/tree.nwk <(nw_topology ${dir}/${query}/apples.tree) -simplify | awk '{printf $2}'`
        n2=`compareTrees.missingBranch ${dir}/tree.nwk <(nw_topology ${dir}/${query}/backbone_app.tree) -simplify | awk '{printf $2}'`
      
       python -c "print ($n1-$n2)" >> $R1

        guppy tog -o $dir/$query/pp.tree $dir/$query/pp.jplace

        n1=`compareTrees.missingBranch ${dir}/tree.nwk <(nw_topology ${dir}/${query}/pp.tree) -simplify | awk '{printf $2}'`
        n2=`compareTrees.missingBranch ${dir}/tree.nwk <(nw_topology ${dir}/${query}/backbone_pp.tree) -simplify | awk '{printf $2}'`
      
       python -c "print ($n1-$n2)" >> $R2

        guppy tog -o $dir/$query/ppa.tree $dir/$query/ppa_result.jplace

        n1=`compareTrees.missingBranch ${dir}/tree.nwk <(nw_topology ${dir}/${query}/ppa.tree) -simplify | awk '{printf $2}'`
        n2=`compareTrees.missingBranch ${dir}/tree.nwk <(nw_topology ${dir}/${query}/backbone_ppa.tree) -simplify | awk '{printf $2}'`
      
       python -c "print ($n1-$n2)" >> $R3
       
       
       guppy tog -o $dir/$query/pdc.tree $dir/$query/pdc_result.jplace
        n1=`compareTrees.missingBranch ${dir}/tree.nwk <(nw_topology ${dir}/${query}/pdc.tree) -simplify | awk '{printf $2}'`
        n2=`compareTrees.missingBranch ${dir}/tree.nwk <(nw_topology ${dir}/${query}/backbone_pdc.tree) -simplify | awk '{printf $2}'`
       python -c "print ($n1-$n2)" >> $R4

   }
  export -f f; xargs -P 28 -I@ bash -c 'f @' < ${dir}/queries.txt
