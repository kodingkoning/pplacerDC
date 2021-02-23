from tree import PhylogeneticTree
import dendropy
import subprocess
import script_executor as se
import string
import random

DEBUG = False

# TODO: replace use of dendropy
def read_tree(tree_file):
    handle = open(tree_file, "r")
    tree = PhylogeneticTree(
      dendropy.Tree.get_from_stream(handle,
                                    schema="newick",
                                    preserve_underscores=True))
    return tree

