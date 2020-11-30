import dendropy
import subprocess
import script_executor as se
import string
import random

DEBUG = False

def read_tree(tree_file):
    handle = open(tree_file, "r")
    tree = dendropy.Tree.get_from_stream(handle,schema="newick",preserve_underscores=True)
    return tree

def read_list(fileName):
    with open(fileName) as fileHandle:
        content = fileHandle.readlines()
    content=[x.strip() for x in content]
    return content
def find_matching_node(subTree, tree, querySequence, siblingToQuery):
  def find_valid_leaf_node(x):
    if x.is_leaf():
      if x.taxon.label != querySequence:
        return x
    for child in x.child_nodes():
      return find_valid_leaf_node(child)
  validLeaf = find_valid_leaf_node(siblingToQuery)
  assert validLeaf
  mainTreeLeaf = None
  for leaf in tree.leaf_nodes():
    if leaf.taxon.label == validLeaf.taxon.label:
      mainTreeLeaf = leaf
      break
  assert mainTreeLeaf
  def advance_up(nodeSubTree, nodeTree, goal):
    if nodeSubTree == goal:
      if DEBUG: print("Found target")
      return nodeTree
    if nodeSubTree.parent_node and nodeTree.parent_node:
      if DEBUG: print("Did not find target, advancing up a level")
      return advance_up(nodeSubTree.parent_node, nodeTree.parent_node, goal)
  target = advance_up(validLeaf, mainTreeLeaf, siblingToQuery)
  return target

def modify_backbone_tree_with_placement(resultTree, backBoneTree, querySequence):
  queryNode = None
  for leaf in resultTree.leaf_nodes():
    if leaf.taxon.label == querySequence:
      queryNode = leaf
      break
  if queryNode == None:
    print("Expected non null queryNode")
 
  assert len(queryNode.sibling_nodes()) == 1, "Expected query node to only have a single sibling"
  siblingToQuery = queryNode.sibling_nodes()[0]
  siblingToQuery = find_matching_node(resultTree, backBoneTree, querySequence, siblingToQuery)
  parent = siblingToQuery._get_parent_node()

  """
  Currently have:

  x---------------------------x
  parent                      eventual siblingToQuery
  Transform into:
           q  query node
           |
           |
           |
  x--------x------------------x
  parent   unmarked node      siblingToQuery
  """

  "Create unmarked node, add it into the main tree"
  unmarkedNode = dendropy.Node(label="unmarkedNode")

  "Disconnect parent from siblingToQuery"
  parent.remove_child(siblingToQuery)
  parent.add_child(unmarkedNode)
  unmarkedNode.add_child(queryNode)
  unmarkedNode.add_child(siblingToQuery)
# Adapted from code provided by Vlad
def sampleCompact(leaf, size):
    taxons = [] # Don't include query node
    node = leaf.parent_node # Don't include parent of query
                            # this isn't a real node in the backbone
    while len(taxons) < size:
        shuffled_nodes = node.sibling_nodes().copy()
        random.shuffle(shuffled_nodes)
        for bro in shuffled_nodes:
            taxons.extend(collectSubtreeTaxa(bro, size - len(taxons)))
            if len(taxons) >= size:
                return taxons
        
        node = node.parent_node      
    return taxons

def collectSubtreeTaxa(node, numTaxa):
    if node.is_leaf():
        return [node.taxon.label]
    
    taxons = []
    shuffledChildren = node.child_nodes().copy()
    random.shuffle(shuffledChildren)
    for child in shuffledChildren:
        taxons.extend(collectSubtreeTaxa(child, numTaxa - len(taxons)))
        if len(taxons) >= numTaxa:
            return taxons
    return taxons
