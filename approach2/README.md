
# Approach 2: APPLES with pplacer

### Step 1: Run APPLES on backbone tree

APPLES would be run on the query sequence and the entire backbone tree.

### Step 2: Define the region of APPLES' placement

Then, the software will find the "region" around the placement from APPLES. The region would be less than 1000 sequences and would be a subtree of the backbone tree. We still need to specify the algorithm to obtain this region. The region could be defined as the bipartion where side including the placed sequence has fewer than 1000 sequences but is as close to the maximum as possible. This relies on the assumption that the majority of APPLES' average error is from having a small amount of error on each placement, rather than a large amount of error on a few runs where it places the sequence in the wrong region.

### Step 3: Run pplacer on the subtree

Next, the software will run pplacer on the subtree from Step 2 and the query sequence.

### Step 4: Return the result in the backbone tree

As in the other approaches, the final result is the placement of the query sequence in the original backbone tree based on its subtree placement.

# Data

## Testing Data

We still need to choose the data we will use to test and debug.

## Benchmarking Data

For simulated data, we will use the same datasets as APPLES used from Guo et al. (2009), and use the RNASim-VS sets. These sets are of size 500, 1000, 5000, 10,000, 50,000, 100,000 and 200,000.

For biological data, we will need advice on selecting appropriate datasets.
