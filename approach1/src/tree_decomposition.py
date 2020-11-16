from sepp_utils.alignment import MutableAlignment
from sepp_utils.tree import PhylogeneticTree
import dendropy
import subprocess
import os # TODO: may not be portable
import numpy as np
import script_executor as se
import string
import random

DEBUG = False

def read_tree(tree_file):
    tree = PhylogeneticTree(
      dendropy.Tree.get_from_stream(tree_file,
                                    schema="newick",
                                    preserve_underscores=True))
    return tree

def read_alignment_and_tree(alignment_file, tree_file):
    alignment = MutableAlignment()
    alignment.read_file_object(alignment_file)
    tree = PhylogeneticTree(
      dendropy.Tree.get_from_stream(tree_file,
                                    schema="newick",
                                    preserve_underscores=True))
    return (alignment, tree)

# TODO: system calls may not be reliable on all systems
# Useful for pruning a query sequence from the tree for leave-one-out
# or leave-many-out studies
def read_list(fileName):
    with open(fileName) as fileHandle:
        content = fileHandle.readlines()
    content=[x.strip() for x in content]
    return content

def label_internal_nodes_subtree_impl(subTreeNode, node, subTreeNodeToNode):
    if subTreeNode in subTreeNodeToNode:
      return
    if DEBUG: print(f"Mapping {subTreeNode} to {node}")
    subTreeNodeToNode[subTreeNode] = node
    if subTreeNode.adjacent_nodes():
      subTreeParent = subTreeNode.adjacent_nodes()[0]
      parent = node.adjacent_nodes()[0]
      if DEBUG: print(f"Mapping {subTreeNode} to {node}")
      label_internal_nodes_subtree_impl(subTreeParent, parent, subTreeNodeToNode)

def label_internal_nodes_subtree(subTree, tree, subTreeNodeToNode):
    visited=()
    visitedSubTree=()
    if DEBUG: print("Labeling internal nodes...")
    for leaf in tree.leaf_nodes():
      for subTreeLeaf in subTree.leaf_nodes():
        if DEBUG: print(f"Comparing {leaf.taxon.label} and {subTreeLeaf.taxon.label}")
        if leaf.taxon.label == subTreeLeaf.taxon.label:
          if DEBUG: print(f"Taxons matched for {leaf} and {subTreeLeaf}")
          subTreeParent = subTreeLeaf.adjacent_nodes()[0]
          parent = leaf.adjacent_nodes()[0]
          label_internal_nodes_subtree_impl(subTreeParent, parent, subTreeNodeToNode)
def validate_result_tree(treeWithPlacement, tree, querySequence):
    taxaTreeWithPlacement = [i.taxon.label for i in treeWithPlacement.leaf_nodes()]
    treeTaxa = [i.taxon.label for i in tree.leaf_nodes()]
    for taxa in taxaTreeWithPlacement:
      expression = taxa in treeTaxa or taxa == querySequence
      assert expression, "Expected {taxa} to be in the main tree, or match the query sequence {querySequence}"
      if DEBUG: print(f"{expression}")
