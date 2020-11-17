from tree_utils import *
from script_executor import *
from util import *
import uuid
import numpy as np
import shutil
import concurrent
import concurrent.futures
from multiprocessing.pool import ThreadPool as Pool


########### TODO ##############
# - Clean up this code

VALIDATE = False
DEBUG = False

timer = Timer()
timer.tic("Program execution")

timer.tic("Setup")
tmpdir = str(uuid.uuid4())
os.mkdir(tmpdir)
querySequences = read_list("queries.txt")
#shutil.copyfile("faSomeRecords.py", f"{tmpdir}/faSomeRecords.py")
oldDir = os.getcwd()
os.chdir(f"{tmpdir}")
base_dir = "/home/malachi/work/classes/variable-size/data/500/0"
tree_size = 150
raxml_info_file = f"{base_dir}/RAxML_result.REF"
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

maxSize = 125 # e.g., this will be a tuneable parameter
decomposed_trees = tree.decompose_tree(maxSize, strategy="centroid", minSize=1)
nTrees = len(decomposed_trees.keys())
timer.toc("Setup")


maxScore = -np.inf
bestTree = "bestTree.tre"
raxml_info_file = f"{base_dir}/RAxML_info.REF"


scores = [0 for i in decomposed_trees.keys()]

for i, tree_key in enumerate(decomposed_trees.keys()):
  #threadLocalTimer = t_timers[i]
  # Constants midrun
  outputLocation = f"output-{i}.jplace"
  resultTree = f"mytestoutput-{i}.tre"
  query_alignment_file = f"foobar-{i}.fa" # test <-----
  temporaryResultTree = f"mytestoutput-{i}.tre"
  temporaryBackBoneTree = f"the-result-{i}.tre" # backbone tree with a tentative placement
  tree_object = decomposed_trees[tree_key]
  print(f"Tree ({tree_key}) has {tree_object.count_nodes()} nodes and {tree_object.count_leaves()} leafs!")
  outputTreeFile = f"decomposed-tree-{i}.tre"
  tree_object.get_tree().write(file=open(outputTreeFile, 'w'),
      schema="newick")
  #threadLocalTimer.tic("FASTA pruning")
  generate_fasta_file(tree_object, querySequence, alignment_file, query_alignment_file)
  #threadLocalTimer.toc("FASTA pruning")
  #threadLocalTimer.tic("pplacer")
  run_pplacer(raxml_info_file, outputTreeFile, query_alignment_file, outputLocation)
  #threadLocalTimer.toc("pplacer")
  # overwrite the current file
  #threadLocalTimer.tic("guppy")
  place_sequence_in_subtree(outputLocation, temporaryResultTree)
  #threadLocalTimer.toc("guppy")

  #threadLocalTimer.tic("tree mod")
  resultTree = read_tree(temporaryResultTree).get_tree()
  backBoneTree = read_tree(prunedTree).get_tree()
  backBoneTreeCopy = None
  if VALIDATE: backBoneTreeCopy = dendropy.Tree(backBoneTree)
  modify_backbone_tree_with_placement(resultTree, backBoneTree, querySequence)
  backBoneTree.write(file=open(temporaryBackBoneTree, "w"), schema="newick")
  #threadLocalTimer.toc("tree mod")

  if VALIDATE:
      #threadLocalTimer.tic("tree validation")
      validate_result_tree(theTree, theTreeCopy, querySequence)
      #threadLocalTimer.toc("tree validation")

  #threadLocalTimer.tic("raxml")
  if nTrees == 1:
      scores[i] = 1.0 # only a single tree, don't bother
  else:
      scores[i] = score_raxml(temporaryBackBoneTree, alignment_file)
    
  #threadLocalTimer.toc("raxml")
  if DEBUG: print(f"ML score = {score}")

#timer.tic("Threaded region")
#nThreads = 1
##executor = concurrent.futures.ProcessPoolExecutor(nThreads)
##futures = [executor.submit(run_subtree, (i,tree_key,scores)) for i, tree_key in enumerate(decomposed_trees.keys())]
##concurrent.futures.wait(futures)
#for i, tree_key in enumerate(decomposed_trees.keys()):
#  run_subtree((i,tree_key, scores))
#
#timer.toc("Threaded region")

#for i, tree_key in enumerate(decomposed_trees.keys()):
#  run_subtree(i,tree_key)

# Do maxLoc reduction to find best tree
bestOne = None
maxScore = -np.inf
for i, score in enumerate(scores):
  if score > maxScore:
    maxScore = score
    bestOne = f"the-result-{i}.tre"
# move best on to old dir, remove old dir
shutil.move(bestOne, f"{oldDir}/{bestTree}")
os.chdir(oldDir)
shutil.rmtree(tmpdir)

timer.toc("Program execution")
timer.dump()
