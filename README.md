# pplacerDC: a New Scalable Maximum Likelihood Phylogenetic Placement Method

## Elizabeth Koning, Malachi Phillips, and Tandy Warnow

## How to Use

pplacerDC does not need to be installed, but should be run from the pplacerDC.py script found in approach1/

Example scripts can be found in scripts/

The command to run pplacerDC is:

`pplacerDC.py -m $maxTreeSize -s $RAxML_info.REF7 -t $backboneTree -q $query -r alignment -o outputFile -j $threads`

All options can also be found by running `pplacerDC.py -h`

Options:

- maxTreeSize: the maximum size of a tree to pass to pplacer
- backboneTree: the tree in Newick format
- RAxML_info.REF7: the RAxML info file from RAxML 7.2.6 based on all the sequences in the tree and the query sequence.
- query: the name of the query, as it appears in the alignment
- alignment: the MSA of all the sequences in the tree and the query sequence
- threads: the maximum number of threads to use

## Summary

pplacerDC is a method for phylogenetic placement of genetic sequences into an existing reference tree. It is based on pplacer, a highly accurate placement method, but that cannot reliably place sequences into trees beyond 5,000 taxa. pplacerDC improves on pplacer by using a divide and conquer approach to use the accuracy of pplacer in trees that would otherwise be too large for pplacer, and is tested on trees up to 100,000 taxa.

## Requirements

1. Python >= 3.7
2. pplacer v1.1.alpha19 (https://github.com/matsen/pplacer/releases/tag/v1.1.alpha19)
3. raxml-ng (https://github.com/amkozlov/raxml-ng.git) 

## Recommended Software for Pre-processing

1. Newick Utilities 1.6 (http://cegg.unige.ch/newick\_utils)
2. RAxML 7.2.6 (https://cme.h-its.org/exelixis/web/software/raxml)

## Citing

When using pplacerDC, please cite [this paper](https://doi.org/10.1145/3459930.3469516):

Elizabeth Koning, Malachi Phillips, and Tandy Warnow. 2021. PplacerDC: a new scalable phylogenetic placement method. In _Proceedings of the 12th ACM Conference on Bioinformatics, Computational Biology, and Health Informatics (BCB '21)_. Association for Computing Machinery, New York, NY, USA, Article 3, 1–9. [https://doi.org/10.1145/3459930.3469516](https://doi.org/10.1145/3459930.3469516)


