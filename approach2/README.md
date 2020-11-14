
# Approach 2: APPLES with pplacer

### Step 1: Run APPLES on backbone tree

APPLES would be run on the query sequence and the entire backbone tree.

### Step 2: Define the region of APPLES' placement

Then, the software will find the "region" around the placement from APPLES. The region would be less than 1000 sequences and would be a subtree of the backbone tree. We still need to specify the algorithm to obtain this region. The region could be defined as the bipartion where side including the placed sequence has fewer than 1000 sequences but is as close to the maximum as possible. This relies on the assumption that the majority of APPLES' average error is from having a small amount of error on each placement, rather than a large amount of error on a few runs where it places the sequence in the wrong region.

### Step 3: Run pplacer on the subtree

Next, the software will run pplacer on the subtree from Step 2 and the query sequence.

### Step 4: Return the result in the backbone tree

As in the other approaches, the final result is the placement of the query sequence in the original backbone tree based on its subtree placement.


### Running the software

In order to run APPLES on the training dataset, the 1000M1 dataset has already been prepared and is ready to run place_app.sh. The scripts get_seqs.sh and data_prep_app.sh both prepared the data for APPLES and the APPLES+pplacer approach. The script to run the combined approach is still in progress, but will be at place_ppa.sh. (Final name yet to be decided.)

### Notes about commands run

The true alignment is used for both the backbone tree estimatio and the placement. For 1000M1, the true alignment is called rose.aln.true.fasta.

For the backbone tree estimation, as in the supplementary materials for APPLES's RNAsim-VS data, we used FastTree (they used FastTreeMP, which is parallelized using OpenMP and is non-deterministic, while the sequential version is deterministic) for estimating the backbone topology, using the following command and FastTree version 2.1.11 SSE3:

```FastTree -nosupport -gtr -gamma -nt -log tree.log < aln_dna.fa > tree.nwk```

We used rose.aln.true.fasta for aln_dna.fa, and the output tree is tree.nwk.
