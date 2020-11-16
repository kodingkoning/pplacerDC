from tree_utils import *
from script_executor import *
from util import *
import uuid
import numpy as np
import shutil


########### TODO ##############
# - Clean up this code

VALIDATE = False
DEBUG = False

timer = Timer()

timer.tic("Setup")

os.chdir("test")
base_dir = "/home/malachi/work/classes/variable-size/data/500/0"
tree_size = 500
raxml_info_file = f"{base_dir}/RAxML_result.REF"
querySequences = read_list("queries.txt")
assert len(querySequences) == 1, "Expected only a single query sequence!"
querySequence = querySequences[0]
# generate backbone tree by removing query sequence
prunedTree = "pruned-tree.tre"
se.prune_query(raxml_info_file, querySequence, prunedTree)

tree_file_handle = open(prunedTree, "r")
alignment_file = f"{base_dir}/aln_dna.fa"
alignment_file_handle = open(f"{base_dir}/aln_dna.fa", "r")
alignment, tree = read_alignment_and_tree(alignment_file_handle,
     tree_file_handle)

maxSize = 250 # e.g., this will be a tuneable parameter
decomposed_trees = tree.decompose_tree(maxSize, strategy="centroid", minSize=1)
timer.toc("Setup")

# Constants midrun
outputLocation = "output.jplace"
resultTree = "mytestoutput.tre"
raxml_info_file = f"{base_dir}/RAxML_info.REF"
query_alignment_file = "foobar.fa" # test <-----
temporaryResultTree = "mytestoutput.tre"
temporaryBackBoneTree = "the-result.tre" # backbone tree with a tentative placement
bestTree = "bestTree.tre"

maxScore = -np.inf

for i, tree_key in enumerate(decomposed_trees.keys()):
  if i > 0:
    break
  tree_object = decomposed_trees[tree_key]
  print(f"Tree ({tree_key}) has {tree_object.count_nodes()} nodes and {tree_object.count_leaves()} leafs!")
  outputTreeFile = f"decomposed-tree-{tree_key}.tre"
  tree_object.get_tree().write(file=open(outputTreeFile, 'w'),
      schema="newick")
  timer.tic("FASTA pruning")
  generate_fasta_file(tree_object, querySequence, alignment_file, query_alignment_file)
  timer.toc("FASTA pruning")
  timer.tic("pplacer")
  run_pplacer(raxml_info_file, outputTreeFile, alignment_file, query_alignment_file, outputLocation)
  timer.toc("pplacer")
  # overwrite the current file
  timer.tic("guppy")
  place_sequence_in_subtree(outputLocation, temporaryResultTree)
  timer.toc("guppy")

  timer.tic("tree mod")
  resultTree = read_tree(temporaryResultTree).get_tree()
  backBoneTree = read_tree(prunedTree).get_tree()
  backBoneTreeCopy = None
  if VALIDATE: backBoneTreeCopy = dendropy.Tree(backBoneTree)
  modify_backbone_tree_with_placement(resultTree, backBoneTree, querySequence)
  backBoneTree.write(file=open(temporaryBackBoneTree, "w"), schema="newick")
  timer.toc("tree mod")

  if VALIDATE:
      timer.tic("tree validation")
      validate_result_tree(theTree, theTreeCopy, querySequence)
      timer.toc("tree validation")

  timer.tic("raxml")
  score = score_raxml(temporaryBackBoneTree, alignment_file)
  if score > maxScore:
    maxScore = score
    shutil.move(temporaryBackBoneTree,bestTree)
    
  timer.toc("raxml")
  if DEBUG: print(f"ML score = {score}")




timer.dump()
