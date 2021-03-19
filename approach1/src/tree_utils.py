from tree import PhylogeneticTree
import dendropy #TODO: delete this when possible
import treeswift
import subprocess
import script_executor as se
import string
import random

DEBUG = False

def read_tree(tree_file):
    tree = treeswift.read_tree(tree_file, "Newick")
    return PhylogeneticTree(tree)
def read_list(fileName):
    with open(fileName) as fileHandle:
        content = fileHandle.readlines()
    content=[x.strip() for x in content]
    return content

def validate_result_tree(treeWithPlacement, tree, querySequence):
    taxaTreeWithPlacement = [label for i in treeWithPlacement.labels(leaves=True,internal=False)]
    treeTaxa = [label for i in tree.labels(leaves=True,internal=False)]
    for taxa in taxaTreeWithPlacement:
      expression = taxa in treeTaxa or taxa == querySequence
      assert expression, f"Expected {taxa} to be in the main tree, or match the query sequence {querySequence}"
      if DEBUG: print(f"{expression}")

def find_matching_node(subTree, tree, querySequence, siblingToQuery):
  def find_valid_leaf_node(x):
    if x.is_leaf():
      if x.get_label() != querySequence:
        return x
    for child in x.child_nodes():
      return find_valid_leaf_node(child)
  validLeaf = find_valid_leaf_node(siblingToQuery)
  assert validLeaf
  mainTreeLeaf = None
  for leaf in tree.traverse_inorder(leaves=True,internal=False):
    if leaf.get_label() == validLeaf.get_label():
      mainTreeLeaf = leaf
      break
  assert mainTreeLeaf
  def advance_up(nodeSubTree, nodeTree, goal):
    if nodeSubTree == goal: # TODO -- make sure that the equality operator works okay with the TreeSwift leaves now
      if DEBUG: print("Found target")
      return nodeTree
    if nodeSubTree.get_parent() and nodeTree.get_parent():
      if DEBUG: print("Did not find target, advancing up a level")
      return advance_up(nodeSubTree.get_parent(), nodeTree.get_parent(), goal)
  target = advance_up(validLeaf, mainTreeLeaf, siblingToQuery)
  return target

def sibling_nodes(node):
    return node.get_parent().child_nodes()

def modify_backbone_tree_with_placement(resultTree, backBoneTree, querySequence):
  queryNode = None
  for leaf in resultTree.traverse_inorder(leaves=True,internal=False):
    if leaf.get_label() == querySequence:
      queryNode = leaf
      break
  if queryNode == None:
    print("Expected non null queryNode")
 
  assert len(queryNode.sibling_nodes()) == 1, "Expected query node to only have a single sibling"
  siblingToQuery = queryNode.sibling_nodes()[0] 
  if siblingToQuery.get_label() == queryNode.get_label(): #TODO: do this better
      siblingToQuery = queryNode.sibling_nodes()[1]
  siblingToQuery = find_matching_node(resultTree, backBoneTree, querySequence, siblingToQuery)
  parent = siblingToQuery.get_parent()

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
  unmarkedNode = treeswift.Node(label="unmarkedNode")

  "Disconnect parent from siblingToQuery"
  parent.remove_child(siblingToQuery)
  parent.add_child(unmarkedNode)
  unmarkedNode.add_child(queryNode)
  unmarkedNode.add_child(siblingToQuery)

