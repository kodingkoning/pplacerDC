# Project Approach 1: Divide and Conquer pplacer

### Step 1: Divide the backbone tree

The software will divide the backbone into disjoint subtrees with fewer than 1000 sequences each. This will be done similarly to the first step of PASTA described in [@PASTA], using centroid decomposition. We will start with the maximum size of the subtrees at 1000 sequences. We will also consider if it may improve the accuracy to have smaller trees, such as 500 sequences, though based on the placement accuracy by reference tree size in the comparison to APPLES in [@balaban_apples_2020], we expect that accuracy will be higher for larger reference trees.

### Step 2: Run pplacer

pplacer will be run independently on each of the subtrees from step 1. If parallelism is later added to the project, then this step would be embarassingly parallel.

### Step 3: Compare pplacer placements

The placements from pplacer runs will be compared to find the tree where the query sequence fits best. The uncertainties or scores from pplacer will not be able to be directly compared, so we will need to either compare the subtrees to each other or find a system to compare the scores so that they reflect the global likelihood of the sequence placement rather than only relative to the other positions in the subtree.

If the comparison of the subtrees is able to be done as a comparison of the subtrees without using all of the pplacer outputs, then that comparison step may replace step 2, and pplacer will only need to be run on the region around the position from running APPLES.

As the query sequence will not fit well in most of the trees, we still need to assess how we will identify which tree is the appropriate one for the tree. If this is best done by comparing the trees to each other rather than the individual placements, then we would use Divide-and-Conquer Approach #2.

### Step 4: Return result in backbone tree

Based on the location of the best score on the best subtree, it will return the corresponding location on the original backbone tree.
