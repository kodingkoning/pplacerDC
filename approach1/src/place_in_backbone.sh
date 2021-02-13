#!/bin/bash

module load python/3
module unload python/2
PATH=$PATH:/home/ekoning2/scratch/pplacer-Linux-v1.1.alpha19

inputTree=$1 	# backbone tree without query leaf
outputTree=$2	# output file
jplaceIn=$3	# jplace file for subtree with placement
jplaceOut=$4	# output for jplace with full tree and placement

# TODO: adapt sed and grep to work in Python for the script, or use subprocess and change later when improving the code more

# 1. Add index numbers to the backbone tree
echo "Step 1"
# NOTE: this only works if the branches have lengths, which is not required
# TODO: make sure that there aren't issues with spaces, etc, in any of the regex matches
sed -E 's/([a-zA-Z0-9]+:[0-9]+\.?[0-9]*)/\1\{\n\}/g' $inputTree | awk '{print $0,NR}' | tr -d '\n ' | sed 's/;.*$/;/' > $outputTree
# TODO: sketch in Python -- insert index # after each one of the matching locations -- split into a list and then insert elements and then patch back together?

# 2. Find the index of the query in the original placement
echo "Step 2"
originalIndex=`./get_placement.py $jplaceIn`
siblingLeaf=`grep -oP "[a-zA-Z0-9:\.]*\{$originalIndex\}" $jplaceIn` # For python, should be able to use regex for this and newIndex 
leafName=`python -c "print(\"$siblingLeaf\".split(':')[0].split('{')[0])"`
newIndex=`grep -oP "$leafName[a-zA-Z0-9:\.]*{\K[0-9]*" $outputTree` 

# 3. combine the jplace info (minus the tree) from $jplaceIn with the tree from $outputTree, and replace the index number in the new jplace file with $newIndex
echo "Step 3"
./modify_jplace.py $jplaceIn $jplaceOut $newIndex $outputTree

# run guppy 
echo "Step 4"
# NOTE: does not work for trees without branch lengths
guppy tog $jplaceOut -o $outputTree

