#!/usr/bin/python3
import sys
import math
import argparse
from tree_utils import *
from script_executor import *
from util import *
import uuid
import shutil
import concurrent
import concurrent.futures
from multiprocessing.pool import ThreadPool as Pool

def run_program(args):
    numThreads = int(args.numThreads)
    maxSubTreeSize = int(args.max)
    outputTree = args.output
    raxml_info_file = args.info
    inputTree = args.tree
    querySequence = args.query
    msaFile = args.ref_msa
    verbose = args.verbose
    VALIDATE = False
    DEBUG = False
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

    inputTree = inputTree if (inputTree.startswith("/") or inputTree.startswith("~")) else f"{oldDir}/{inputTree}"
    msaFile = msaFile if (msaFile.startswith("/") or msaFile.startswith("~")) else f"{oldDir}/{msaFile}"
    raxml_info_file = raxml_info_file if raxml_info_file.startswith("/") or raxml_info_file.startswith("~") else f"{oldDir}/{raxml_info_file}"
    outputTree = outputTree if outputTree.startswith("/") or outputTree.startswith("~") else f"{oldDir}/{outputTree}"
    
    if DEBUG: print("Reading tree...")
    tree = read_tree(inputTree)
    
    if DEBUG: print("Decomposing tree...")
    decomposed_trees = tree.decompose_tree(maxSubTreeSize, strategy="centroid", minSize=1)
    nTrees = len(decomposed_trees.keys())
    if DEBUG: print(f"Has {nTrees} sub trees...")
    timer.toc("Setup")
    
    maxScore = -math.inf
    
    scores = [-math.inf for i in decomposed_trees.keys()]
    
    def run_subtree(item):
      i, tree_key = item
      if DEBUG: print(f"On thread id = {i}")
      outputLocation = f"output-{i}.jplace"
      resultTree = f"mytestoutput-{i}.tre"
      query_alignment_file = f"foobar-{i}.fa"
      temporaryResultTree = f"mytestoutput-{i}.tre"
      temporaryBackBoneTree = f"the-result-{i}.tre"
      tree_object = decomposed_trees[tree_key]
      if DEBUG: print(f"Tree ({tree_key}) has {tree_object.count_nodes()} nodes and {tree_object.count_leaves()} leafs!")
      outputTreeFile = f"decomposed-tree-{i}.tre"
      tree_object.get_tree().write(file=open(outputTreeFile, 'w'),
          schema="newick")
      generate_fasta_file(tree_object, querySequence, msaFile, query_alignment_file)
      run_pplacer(raxml_info_file, outputTreeFile, query_alignment_file, outputLocation)
      # overwrite the current file
      place_sequence_in_subtree(outputLocation, temporaryResultTree)
    
      resultTree = read_tree(temporaryResultTree).get_tree()
      backBoneTree = read_tree(inputTree).get_tree()
      backBoneTreeCopy = None
      if VALIDATE: backBoneTreeCopy = dendropy.Tree(backBoneTree)
      if DEBUG: print(f"Modifying with query sequence {querySequence}")
      modify_backbone_tree_with_placement(resultTree, backBoneTree, querySequence)
      backBoneTree.write(file=open(temporaryBackBoneTree, "w"), schema="newick")
    
      if VALIDATE:
          validate_result_tree(backBoneTree, backBoneTreeCopy, querySequence)
    
      if nTrees == 1:
          scores[i] = 1.0 # only a single tree, don't bother
      else:
          scores[i] = score_raxml(temporaryBackBoneTree, msaFile)
        
      if DEBUG: print(f"ML score = {scores[i]}")

    timer.tic("Threaded region")
    with concurrent.futures.ThreadPoolExecutor(max_workers=numThreads) as executor:
        executor.map(run_subtree, [(i, tree_key) for i, tree_key in enumerate(decomposed_trees.keys())])
    timer.toc("Threaded region")

    print(scores)
    
    # Do maxLoc reduction to find best tree
    bestOne = None
    maxScore = -math.inf
    for i, score in enumerate(scores):
      if score > maxScore:
        maxScore = score
        bestOne = f"the-result-{i}.tre"
    # move best on to old dir, remove old dir
    print(f"Moving file {bestOne} to {outputTree}")
    shutil.move(bestOne, outputTree)
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
    requiredNamed.add_argument('-s', '--info', help='RAxML v7 or v8 info file [NOT raxml-ng]. Path must be relative to directory flag.', required=True)
    requiredNamed.add_argument('-t', '--tree', help='Input tree. Path must be relative to directory flag.', required=True)
    requiredNamed.add_argument('-q', '--query', help='Query taxa to place into tree. Path must be relative to directory flag.', required=True)
    requiredNamed.add_argument('-r', '--ref-msa', help='Reference MSA. Path must be relative to directory flag.', required=True)
    
    args = parser.parse_args()
    run_program(args)
