---
title: CS581 Project Proposal
author: Elizabeth Koning, Malachi Phillips
date: October 22, 2020
geometry: margin=2cm
header-includes: |
    \usepackage[ruled,vlined]{algorithm2e}
nocite: |
  @*
---
# CS 581 Project: Imrpoving Scalability and Accuracy in Phylogenetic Placement

## Elizabeth Koning and Malachi Phillips

Taxonomic identification and phylogenetic profiling software TIPP @nguyen_tipp_2014 
requires solving a phylogenetic placement problem in order
to insert a new alignment into a tree.
Currently, software such as pplacer [@matsen_pplacer_2010] is an accurate,
maximum likelihood phylogenetic placement method.
However, pplacer does not scale to problems larger than 1,000 sequences.
Therefore, alternative phylogenetic placement methods, such as APPLES [@balaban_apples_2020],
were developed for solving phylogenetic placement.
While APPLES has been tested to run on 200,000 sequences, the accuracy of the
method for smaller problems is typically worse than pplacer.
This course project aims to improve phylogenetic placement software already
available by using pplacer in a divide-and-conquer approach.
