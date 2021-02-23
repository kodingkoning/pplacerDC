#!/usr/bin/env python3
import sys
import re
import json
import script_executor as se # for guppy

DEBUG = False

def place_in_backbone(inputTreeFile, outputTreeFile, jplaceInFile, jplaceOutFile):
    with open(inputTreeFile,"r") as f:
        inputTree = f.read().strip()

    # 1. add index numbers to the backbone tree's branches
    inputTree = re.sub(r'([a-zA-Z0-9]*:[0-9]+\.?[0-9]*)', r'\1{\n}', inputTree)
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
    if DEBUG: print(f"originalIndex = {originalIndex}")

    # find the node with the number in the subtree
    index = data['tree'].find(f'{{{originalIndex}}}')
    if DEBUG: print(f"index in subtree = {index}")
    colonIndex = data['tree'][0:index].rfind(':')
    if DEBUG: print(f"colon index in subtree = {colonIndex}")
    if DEBUG: print(data['tree'][0:colonIndex][-20:])
    closeParens=0
    startIndex = colonIndex
    # Formats: right before the colon is either L102 or (L102:0.0002,L103:0.003)
    for c in reversed(data['tree'][0:startIndex]):
        if c == '(':
            closeParens -= 1
        elif c == ')':
            closeParens += 1
        if closeParens == 0:
            break
        startIndex -= 1
    siblingNames = []
    newIndex = -1
    if data['tree'][colonIndex-1] != ')':
        if DEBUG: print("sibling is leaf")
        if DEBUG: print(data['tree'][:startIndex][-30:])
        siblingName = re.search(r'[a-zA-Z0-9]+$', data['tree'][:startIndex]).group(0)
        #the newIndex is the next {x} -- format will be siblingName:#.#*{newIndex}
        if DEBUG: print(rf'siblingName+regex = {siblingName}:[0-9]+\.?[0-9]*\{{([0-9]*)')
        newIndex = re.search(rf'{siblingName}:[0-9]+\.?[0-9]*\{{([0-9]*)',inputTree).group(1)
        #maxIndex = inputTree.find(siblingName)
        #newIndex = re.search(r'\{([0-9]*)\}',inputTree[maxIndex:]).group(1)
    else:
        if DEBUG: print("sibling is NOT leaf")
        # get everything between ( and ) and find all the leaf names in that. Then find the {[0-9]*} after the last one
        if DEBUG: print(data['tree'][startIndex:colonIndex])
        siblingNames = re.findall(r'[^:a-zA-Z0-9{}.]([a-zA-Z0-9]+)', data['tree'][startIndex:colonIndex])
        if DEBUG: print(siblingNames)

        maxIndex = 0
        minIndex = len(inputTree)
        for name in siblingNames:
            idx = inputTree.find(name)
            maxIndex = max(maxIndex, idx)
            minIndex = min(minIndex, idx)
        # after this, max index is going to be at the end of the NAME of the last leaf in the sibling
        # we need to get the the end of the last {} of the last leaf -- this means finding the ) that goes with the first ( 
        if DEBUG: print(f"maxIndex = {maxIndex}")
        maxIndex += inputTree[maxIndex:].find('}')
        if DEBUG: print(f"maxIndex = {maxIndex}")
        if DEBUG: print(f"minIndex = {minIndex}")
        openParens = data['tree'][startIndex:colonIndex].count('(') + 1
        openParens = data['tree'][startIndex:colonIndex].count(')')
        if DEBUG: print(f"openParens = {openParens}") 
        idx = minIndex
        while openParens > 0:
            if inputTree[idx] == ')':
                openParens -= 1
            idx += 1
            #print(f"openParens = {openParens}, idx = {idx}")
        newIndex = re.search(r'\{([0-9]*)\}',inputTree[idx:]).group(1)

        if DEBUG: print(inputTree[minIndex:idx])
        if DEBUG: print(f"before = {inputTree[minIndex-5:minIndex]}")
        if DEBUG: print(f"after = {inputTree[maxIndex:idx+5]}")

    if DEBUG: print(f"newIndex = {newIndex}")

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

