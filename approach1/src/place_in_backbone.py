#!/usr/bin/env python3
import sys
import re
import json
import script_executor as se # for guppy

def place_in_backbone(inputTreeFile, outputTreeFile, jplaceInFile, jplaceOutFile):
    with open(inputTreeFile,"r") as f:
        inputTree = f.read().strip()

    # 1. add index numbers to the backbone tree
    inputTree = re.sub(r'([a-zA-Z0-9]+:[0-9]+\.?[0-9]*)', r'\1{\n}', inputTree)
    inputTree = inputTree.split('\n')
    treeIndices = list(range(len(inputTree)))
    inputTree = [str(val) for pair in zip(inputTree, treeIndices) for val in pair]
    inputTree.pop() # last item isn't a leaf
    inputTree = ''.join(inputTree)

    # 2. find the index of the query in the original placement
    filename = sys.argv[1]
    with open(jplaceInFile,'r') as f:
        data = json.load(f)
    originalIndex = data['placements'][0]['p'][0][1]
    siblingLeaf = re.search(r'[a-zA-Z0-9:\.]*\{{{}}}'.format(originalIndex),json.dumps(data,indent=1)).group(0)
    leafName= siblingLeaf.split(':')[0].split('{')[0]
    newIndex= re.search(r'{0}[a-zA-Z0-9:\.]*{{([0-9]*)}}'.format(leafName), inputTree).group(1)

    # 3. combine the jplace info (minus the tree) from $jplaceIn with the tree from $outputTree, and replace the index number in the new jplace file with $newIndex
    placement_p = data['placements'][0]['p'][0]
    data['placements'][0]['p'] = [placement_p]
    data['placements'][0]['p'][0][1] = int(newIndex)
    data['tree'] = inputTree
    with open(jplaceOutFile,"w") as f:
        json.dump(data,f,indent=1)

    # 4. run guppy 
    # NOTE: does not work for trees without branch lengths
    se.run_guppy_tog(jplaceOutFile, outputTreeFile)

if __name__ == "__main__":
  inputTreeFile=sys.argv[1]     # backbone tree without query leaf
  outputTreeFile=sys.argv[2]    # output file
  jplaceInFile=sys.argv[3]      # jplace file for subtree with placement
  jplaceOutFile=sys.argv[4]     # output for jplace with full tree and placement
  place_in_backbone(inputTreeFile, outputTreeFile, jplaceInFile, jplaceOutFile)

