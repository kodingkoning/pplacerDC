#!/usr/bin/python3
import numpy as np
import sys
import argparse
from tree_utils import *
from script_executor import *
from util import *
import uuid
import numpy as np
import shutil
import concurrent
import concurrent.futures
from multiprocessing.pool import ThreadPool as Pool

VALIDATE = False
DEBUG = False
def run_program(args):
    numThreads = int(args.numThreads)
    maxSubTreeSize = int(args.max)
    outputTree = args.output
    raxml_info_file = args.info
    inputTree = args.tree
    querySequence = args.query
    msaFile = args.ref_msa
    verbose = args.verbose
    if verbose:
      VALIDATE = True
      DEBUG = True
    timer = Timer()
    timer.tic("Program execution")
    
    timer.tic("Setup")
    tmpdir = str(uuid.uuid4())
    os.mkdir(tmpdir)
    oldDir = os.getcwd()
    os.chdir(f"{tmpdir}")
    se.prune_query(raxml_info_file, querySequence, prunedTree)
    
    tree_file_handle = open(inputTree, "r")
    tree = read_tree(tree_file_handle)
    
    decomposed_trees = tree.decompose_tree(maxSize, strategy="centroid", minSize=1)
    nTrees = len(decomposed_trees.keys())
    timer.toc("Setup")
    
    maxScore = -np.inf
    bestTree = "bestTree.tre"
    
    scores = [0 for i in decomposed_trees.keys()]
    
    def run_subtree(item):
      i, tree_key, scores = item
      outputLocation = f"output-{i}.jplace"
      resultTree = f"mytestoutput-{i}.tre"
      query_alignment_file = f"foobar-{i}.fa"
      temporaryResultTree = f"mytestoutput-{i}.tre"
      temporaryBackBoneTree = f"the-result-{i}.tre"
      tree_object = decomposed_trees[tree_key]
      print(f"Tree ({tree_key}) has {tree_object.count_nodes()} nodes and {tree_object.count_leaves()} leafs!")
      outputTreeFile = f"decomposed-tree-{i}.tre"
      tree_object.get_tree().write(file=open(outputTreeFile, 'w'),
          schema="newick")
      generate_fasta_file(tree_object, querySequence, alignment_file, query_alignment_file)
      run_pplacer(raxml_info_file, outputTreeFile, alignment_file, query_alignment_file, outputLocation)
      # overwrite the current file
      place_sequence_in_subtree(outputLocation, temporaryResultTree)
    
      resultTree = read_tree(temporaryResultTree).get_tree()
      backBoneTree = read_tree(prunedTree).get_tree()
      backBoneTreeCopy = None
      if VALIDATE: backBoneTreeCopy = dendropy.Tree(backBoneTree)
      modify_backbone_tree_with_placement(resultTree, backBoneTree, querySequence)
      backBoneTree.write(file=open(temporaryBackBoneTree, "w"), schema="newick")
    
      if VALIDATE:
          validate_result_tree(theTree, theTreeCopy, querySequence)
    
      if nTrees == 1:
          scores[i] = 1.0 # only a single tree, don't bother
      else:
          scores[i] = score_raxml(temporaryBackBoneTree, alignment_file)
        
      if DEBUG: print(f"ML score = {score}")
    
    timer.tic("Threaded region")
    executor = concurrent.futures.ProcessPoolExecutor(numThreads)
    futures = [executor.submit(run_subtree, (i,tree_key,scores)) for i, tree_key in enumerate(decomposed_trees.keys())]
    concurrent.futures.wait(futures)
    timer.toc("Threaded region")
    
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Place a query sequence into a backbone tree.')
    # optional args
    parser.add_argument('-j', '--numThreads', default=1, help='Number of threads (default: 1)')
    parser.add_argument('-m', '--max', default=500, help='Maximum size of subtree to hand to pplacer')
    parser.add_argument('-o', '--output', default='tree.tre', help='Resultant tree with placement')
    parser.add_argument('-v', '--verbose', default=False, help='Run in verbose mode')
    
    requiredNamed = parser.add_argument_group("required named arguments")
    
    # required args
    requiredNamed.add_argument('-s', '--info', help='RAxML v7 or v8 info file [NOT raxml-ng]', required=True)
    requiredNamed.add_argument('-t', '--tree', help='Input tree', required=True)
    requiredNamed.add_argument('-q', '--query', help='Query taxa to place into tree.', required=True)
    requiredNamed.add_argument('-r', '--ref-msa', help='Reference MSA', required=True)
    
    args = parser.parse_args()
    run_program(args)
