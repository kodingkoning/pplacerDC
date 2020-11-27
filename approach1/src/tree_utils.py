from sepp_utils.tree import PhylogeneticTree
import dendropy
import subprocess
import script_executor as se
import string
import random

DEBUG = False

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

def label_internal_nodes_subtree_impl(subTreeNode, node, queryNode, subTreeNodeToNode):
    if subTreeNode in subTreeNodeToNode:
      return # nothing to do
    if DEBUG: print(f"Mapping {subTreeNode} to {node}")
    subTreeNodeToNode[subTreeNode] = node
    queryNodeParent = queryNode._get_parent_node()
    if subTreeNode._get_parent_node() and node._get_parent_node():
      subTreeParent = subTreeNode._get_parent_node()
      parent = node._get_parent_node()
      # skip parent of the query node, since that doesn't exist in the backBoneTree
      if subTreeParent == queryNodeParent:
        if DEBUG: print("Skipping parent of query node in sub tree")
        assert subTreeParent._get_parent_node()
        subTreeParent = subTreeParent._get_parent_node()
      label_internal_nodes_subtree_impl(subTreeParent, parent, queryNode, subTreeNodeToNode)

def label_internal_nodes_subtree(subTree, tree, querySequence, subTreeNodeToNode):
    if DEBUG: print("Labeling internal nodes...")
    if DEBUG: print(f"subTree has {len(subTree.leaf_nodes())} leaf nodes and {len(subTree.nodes())} nodes")
    # Grab queryNode
    queryNode = None
    for leaf in subTree.leaf_nodes():
      if leaf.taxon.label == querySequence:
        queryNode = leaf
        break
    assert queryNode, "Expected queryNode to be non-null"

    for leaf in tree.leaf_nodes():
      for subTreeLeaf in subTree.leaf_nodes():
        #if DEBUG: print(f"Comparing {leaf.taxon.label} and {subTreeLeaf.taxon.label}")
        if leaf.taxon.label == subTreeLeaf.taxon.label:
          #if DEBUG: print(f"Taxons matched for {leaf} and {subTreeLeaf}")
          subTreeParent = subTreeLeaf._get_parent_node()
          parent = leaf._get_parent_node()
          subTreeNodeToNode[subTreeLeaf] = leaf
          label_internal_nodes_subtree_impl(subTreeParent, parent, queryNode, subTreeNodeToNode)
    if DEBUG: print(f"subTreeNodeToNode has {len(subTreeNodeToNode.keys())}")
def validate_result_tree(treeWithPlacement, tree, querySequence):
    taxaTreeWithPlacement = [i.taxon.label for i in treeWithPlacement.leaf_nodes()]
    treeTaxa = [i.taxon.label for i in tree.leaf_nodes()]
    for taxa in taxaTreeWithPlacement:
      expression = taxa in treeTaxa or taxa == querySequence
      assert expression, f"Expected {taxa} to be in the main tree, or match the query sequence {querySequence}"
      if DEBUG: print(f"{expression}")

def modify_backbone_tree_with_placement(resultTree, backBoneTree, querySequence):
  subTreeNodeToNode={}
  label_internal_nodes_subtree(resultTree, backBoneTree, querySequence, subTreeNodeToNode)
  queryNode = None
  for leaf in resultTree.leaf_nodes():
    if leaf.taxon.label == querySequence:
      queryNode = leaf
      break
  if queryNode == None:
    print("Expected non null queryNode")
 
  assert len(queryNode.sibling_nodes()) == 1, "Expected query node to only have a single sibling"
  siblingToQuery = queryNode.sibling_nodes()[0]
  siblingToQuery = subTreeNodeToNode[siblingToQuery]
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
