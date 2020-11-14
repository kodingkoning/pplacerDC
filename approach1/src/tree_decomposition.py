from sepp_utils.alignment import MutableAlignment
from sepp_utils.tree import PhylogeneticTree
import dendropy
import subprocess
import os # TODO: may not be portable
import numpy as np
import script_executor as se

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
