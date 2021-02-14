from tree import PhylogeneticTree
import dendropy
import subprocess
import script_executor as se
import string
import random

DEBUG = False

# TODO: check if dendropy-related methods are used and either delete unneeded functions or replace with TreeSwift 

def read_tree(tree_file):
    handle = open(tree_file, "r")
    tree = PhylogeneticTree(
      dendropy.Tree.get_from_stream(handle,
                                    schema="newick",
                                    preserve_underscores=True))
    return tree
def read_list(fileName):
    with open(fileName) as fileHandle:
        content = fileHandle.readlines()
    content=[x.strip() for x in content]
    return content

def validate_result_tree(treeWithPlacement, tree, querySequence):
    taxaTreeWithPlacement = [i.taxon.label for i in treeWithPlacement.leaf_nodes()]
    treeTaxa = [i.taxon.label for i in tree.leaf_nodes()]
    for taxa in taxaTreeWithPlacement:
      expression = taxa in treeTaxa or taxa == querySequence
      assert expression, f"Expected {taxa} to be in the main tree, or match the query sequence {querySequence}"
      if DEBUG: print(f"{expression}")

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
