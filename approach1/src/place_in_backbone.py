#!/usr/bin/env python3
import sys
import re
import json

inputTreeFile=sys.argv[1] 	# backbone tree without query leaf
outputTreeFile=sys.argv[2]	# output file
jplaceInFile=sys.argv[3]	# jplace file for subtree with placement
jplaceOutFile=sys.argv[4]	# output for jplace with full tree and placement

# 1. Add index numbers to the backbone tree
# NOTE: this only works if the branches have lengths, which is not required
# TODO: make sure that there aren't issues with spaces, etc, in any of the regex matches
#sed -E 's/([a-zA-Z0-9]+:[0-9]+\.?[0-9]*)/\1\{\n\}/g' $inputTree | awk '{print $0,NR}' | tr -d '\n ' | sed 's/;.*$/;/' > $outputTree

with open(inputTreeFile,"r") as f:
    inputTree = f.read().strip()

inputTree = re.sub(r'([a-zA-Z0-9]+:[0-9]+\.?[0-9]*)', r'\1{\n}', inputTree)
inputTree = inputTree.split('\n')
#inputTree = re.split(r'([^a-zA-Z0-9]+[a-zA-Z0-9]+:[0-9]+\.?[0-9]*)', inputTree) # NOTE: this version has double the splits because it gets empty strings
#inputTree = re.split(r'([a-zA-Z0-9]+:[0-9]+\.?[0-9]*)', inputTree) 
#print(len(inputTree))
treeIndices = list(range(len(inputTree)))
#inputTree = list(zip(inputTree, treeIndices))
inputTree = [str(val) for pair in zip(inputTree, treeIndices) for val in pair]
inputTree.pop() # last item isn't a leaf
#print(inputTree)
inputTree = ''.join(inputTree)

#with open(outputTreeFile,'w') as f:
#    f.write(str(inputTree))

# 2. Find the index of the query in the original placement
filename = sys.argv[1]
with open(jplaceInFile,'r') as f:
    data = json.load(f)
originalIndex = data['placements'][0]['p'][0][1]
# originalIndex=`./get_placement.py $jplaceIn`
siblingLeaf = re.search(r'[a-zA-Z0-9:\.]*\{{{}}}'.format(originalIndex),json.dumps(data,indent=1)).group(0)
#siblingLeaf=`grep -oP "[a-zA-Z0-9:\.]*\{$originalIndex\}" $jplaceIn` # For python, should be able to use regex for this and newIndex 
#leafName=`python -c "print(\"$siblingLeaf\".split(':')[0].split('{')[0])"`
leafName= siblingLeaf.split(':')[0].split('{')[0]
#print(siblingLeaf)
#print(leafName)
#newIndex=`grep -oP "$leafName[a-zA-Z0-9:\.]*{\K[0-9]*" $outputTree`
#indexRegex = r'{0}[a-zA-Z0-9:\.]*{{[0-9]*'.format(leafName)
newIndex= re.search(r'{0}[a-zA-Z0-9:\.]*{{([0-9]*)}}'.format(leafName), inputTree).group(1)
#print(newIndex)

# 3. combine the jplace info (minus the tree) from $jplaceIn with the tree from $outputTree, and replace the index number in the new jplace file with $newIndex
#./modify_jplace.py $jplaceIn $jplaceOut $newIndex $outputTree
placement_p = data['placements'][0]['p'][0]
data['placements'][0]['p'] = [placement_p]
data['placements'][0]['p'][0][1] = int(newIndex)
data['tree'] = inputTree
#with open(outputTreeFile,"r") as f:
#    data['tree'] = f.readlines()[0]
with open(jplaceOutFile,"w") as f:
    json.dump(data,f,indent=1)

# run guppy 
# NOTE: does not work for trees without branch lengths
#guppy tog $jplaceOut -o $outputTree

