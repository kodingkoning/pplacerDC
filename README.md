# pplacerDC: a New Scalable Maximum Likelihood Phylogenetic Placement Method

## Elizabeth Koning, Malachi Phillips, and Tandy Warnow

## Summary

pplacerDC is a method for phylogenetic placement of genetic sequences into an existing reference tree. It is based on pplacer, a highly accurate placement method, but that cannot reliably place sequences into trees beyond 5,000 taxa. pplacerDC improves on pplacer by using a divide and conquer approach to use the accuracy of pplacer in trees that would otherwise be too large for pplacer, and is tested on trees up to 100,000 taxa.

## Requirements

1. Python >= 3.7
2. pplacer v1.1.alpha19 (https://github.com/matsen/pplacer/releases/tag/v1.1.alpha19)
3. raxml-ng (https://github.com/amkozlov/raxml-ng.git) 

## Recommended Software for Pre-processing

1. Newick Utilities 1.6 (http://cegg.unige.ch/newick\_utils)
2. RAxML 7.2.6 (https://cme.h-its.org/exelixis/web/software/raxml)

## How to Use

pplacerDC.py can be found in approach1/

Example scripts can be found in scripts/

All options can be found by running `pplacerDC.py -h`

