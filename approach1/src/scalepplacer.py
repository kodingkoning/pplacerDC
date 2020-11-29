#!/usr/bin/python3
import sys
import math
import argparse
from tree_utils import *
import script_executor as se
import uuid
from util import *
import shutil
import concurrent
import concurrent.futures
from multiprocessing.pool import ThreadPool as Pool
import os.path # for debugging
from os import path
DEBUG_THREAD = False
class ThreadLocalStorage(object):
  def __init__(self, tid, outputTree, raxml_info_file, inputTree, querySequence, msaFile):
    i = tid
    self.outputLocation = f"output-{i}.jplace"
    self.resultTree = f"mytestoutput-{i}.tre"
    self.query_alignment_file = f"foobar-{i}.fa"
    self.temporaryResultTree = f"mytestoutput-{i}.tre"
    self.temporaryBackBoneTree = f"the-result-{i}.tre"
    self.outputTreeFile = f"decomposed-tree-{i}.tre"

    self.outputTree = outputTree
    self.raxml_info_file = raxml_info_file
    self.inputTree = inputTree
    self.querySequence = querySequence
    self.msaFile = msaFile
    self.score = -math.inf

def make_fasta_files(item):
  tid, tree_object, threadLocalStorage = item
  if DEBUG_THREAD: print(f"Making fasta files on thread {tid}")
  tree_object.get_tree().write(file=open(threadLocalStorage.outputTreeFile, 'w'),
      schema="newick")
  se.generate_fasta_file(tree_object, threadLocalStorage.querySequence, threadLocalStorage.msaFile, threadLocalStorage.query_alignment_file)
  return
def execute_pplacer(item):
  tid, tree_object, threadLocalStorage = item
  if DEBUG_THREAD: print(f"Running pplacer on thread {tid}")
  se.run_pplacer(threadLocalStorage.raxml_info_file, threadLocalStorage.outputTreeFile, threadLocalStorage.query_alignment_file, threadLocalStorage.outputLocation)
  return
def place_query(item):
  tid, tree_object, threadLocalStorage = item
  if DEBUG_THREAD: print(f"Placing query sequence on thread {tid}")
  se.place_sequence_in_subtree(threadLocalStorage.outputLocation, threadLocalStorage.temporaryResultTree)
  return
def modify_trees(item):
  tid, tree_object, threadLocalStorage = item
  if DEBUG_THREAD: print(f"Modifying trees on thread {tid}")
  resultTree = read_tree(threadLocalStorage.temporaryResultTree).get_tree()
  backBoneTree = read_tree(threadLocalStorage.inputTree).get_tree()
  modify_backbone_tree_with_placement(resultTree, backBoneTree, threadLocalStorage.querySequence)
  backBoneTree.write(file=open(threadLocalStorage.temporaryBackBoneTree, "w"), schema="newick")
  return
def score_trees(item):
  tid, tree_object, threadLocalStorage = item
  if DEBUG_THREAD: print(f"Scoring trees on thread {tid}")
  se.score_raxml(threadLocalStorage.temporaryBackBoneTree, threadLocalStorage.msaFile, tid)
  return

def parse_score(tid):
    regex = "Final LogLikelihood: (.+)"
    raxmlLog = f"raxml-prefix-{tid}.raxml.log"
    score = se.field_by_regex(regex, raxmlLog)[0]
    return score

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
    if not os.path.exists("temp_files"):
      os.mkdir("temp_files")
    tmpdir = "temp_files/"+str(uuid.uuid4())
    print(tmpdir)
    os.mkdir(tmpdir)
    oldDir = os.getcwd()
    os.chdir(f"{tmpdir}")

    inputTree = inputTree if (inputTree.startswith("/") or inputTree.startswith("~")) else f"{oldDir}/{inputTree}"
    msaFile = msaFile if (msaFile.startswith("/") or msaFile.startswith("~")) else f"{oldDir}/{msaFile}"
    raxml_info_file = raxml_info_file if raxml_info_file.startswith("/") or raxml_info_file.startswith("~") else f"{oldDir}/{raxml_info_file}"
    outputTree = outputTree if outputTree.startswith("/") or outputTree.startswith("~") else f"{oldDir}/{outputTree}"
    
    if DEBUG:
      print(f"inputTree: "+ str(path.exists(inputTree)))
      print(f"tmpdir: "+ str(path.isdir(tmpdir)))
      print(f"msaFile: "+ str(path.exists(msaFile)))
      print(f"raxml_info_file: "+ str(path.exists(raxml_info_file)))
    if DEBUG: print("Reading tree...")
    tree = read_tree(inputTree)
    
    if DEBUG: print("Decomposing tree...")
    decomposed_trees = tree.decompose_tree(maxSubTreeSize, strategy="centroid", minSize=1)
    nTrees = len(decomposed_trees.keys())
    if DEBUG: print(f"Has {nTrees} sub trees...")
    timer.toc("Setup")
    
    maxScore = -math.inf
    
    scores = [-math.inf for i in decomposed_trees.keys()]

    # Thread local data
    threadData = [ThreadLocalStorage(tid,
      outputTree,
      raxml_info_file,
      inputTree,
      querySequence,
      msaFile
      ) for tid, _ in enumerate(decomposed_trees.keys())]

    # Thread local variants do *one* process at a time
    SERIAL = False # This provides better debugging output
    if not SERIAL:
      timer.tic("Generating fasta files for pplacer...")
      if DEBUG: print("Generating fasta files for pplacer...")
      with concurrent.futures.ProcessPoolExecutor(max_workers=numThreads) as executor:
          executor.map(make_fasta_files, [(i, decomposed_trees[tree_key], threadData[i]) for i, tree_key in enumerate(decomposed_trees.keys())])
      timer.toc("Generating fasta files for pplacer...")

      timer.tic("Running pplacer...")
      if DEBUG: print("Running pplacer...")
      with concurrent.futures.ProcessPoolExecutor(max_workers=numThreads) as executor:
          executor.map(execute_pplacer, [(i, decomposed_trees[tree_key], threadData[i]) for i, tree_key in enumerate(decomposed_trees.keys())])
      timer.toc("Running pplacer...")

      timer.tic("Placing query into subtree...")
      if DEBUG: print("Placing query into subtree...")
      with concurrent.futures.ProcessPoolExecutor(max_workers=numThreads) as executor:
          executor.map(place_query, [(i, decomposed_trees[tree_key], threadData[i]) for i, tree_key in enumerate(decomposed_trees.keys())])
      timer.toc("Placing query into subtree...")

      timer.tic("Modify main tree...")
      if DEBUG: print("Modifying main tree...")
      with concurrent.futures.ProcessPoolExecutor(max_workers=numThreads) as executor:
          executor.map(modify_trees, [(i, decomposed_trees[tree_key], threadData[i]) for i, tree_key in enumerate(decomposed_trees.keys())])
      timer.toc("Modify main tree...")

      timer.tic("Scoring trees...")
      if DEBUG: print("Scoring trees...")
      with concurrent.futures.ProcessPoolExecutor(max_workers=numThreads) as executor:
          executor.map(score_trees, [[i, decomposed_trees[tree_key], threadData[i]] for i, tree_key in enumerate(decomposed_trees.keys())])
      timer.toc("Scoring trees...")
    else:
      timer.tic("Generating fasta files for pplacer...")
      if DEBUG: print("Generating fasta files for pplacer...")
      for i, tree_key in enumerate(decomposed_trees.keys()):
        make_fasta_files((i, decomposed_trees[tree_key], threadData[i]))
      timer.toc("Generating fasta files for pplacer...")

      timer.tic("Running pplacer...")
      if DEBUG: print("Running pplacer...")
      for i, tree_key in enumerate(decomposed_trees.keys()):
        execute_pplacer((i, decomposed_trees[tree_key], threadData[i]))
      timer.toc("Running pplacer...")

      timer.tic("Placing query into subtree...")
      if DEBUG: print("Placing query into subtree...")
      for i, tree_key in enumerate(decomposed_trees.keys()):
        place_query((i, decomposed_trees[tree_key], threadData[i]))
      timer.toc("Placing query into subtree...")

      timer.tic("Modify main tree...")
      if DEBUG: print("Modifying main tree...")
      for i, tree_key in enumerate(decomposed_trees.keys()):
        modify_trees((i, decomposed_trees[tree_key], threadData[i]))
      timer.toc("Modify main tree...")

      timer.tic("Scoring trees...")
      if DEBUG: print("Scoring trees...")
      for i, tree_key in enumerate(decomposed_trees.keys()):
        score_trees((i, decomposed_trees[tree_key], threadData[i]))
      timer.toc("Scoring trees...")

    for i, threadStorage in enumerate(threadData):
      scores[i] = parse_score(i)

    if DEBUG: print(scores)
    
    # Do maxLoc reduction to find best tree
    bestOne = None
    maxScore = -math.inf
    for i, score in enumerate(scores):
      if score > maxScore:
        maxScore = score
        bestOne = f"the-result-{i}.tre"
    # move best on to old dir, remove old dir
    if DEBUG: print(f"Moving file {bestOne} to {outputTree}")
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
