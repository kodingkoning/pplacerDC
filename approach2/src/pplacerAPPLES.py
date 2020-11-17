#!/usr/bin/python3
import dendropy
from script_executor import *
from tree_utils import *
import uuid
import os
from util import *
import numpy
import argparse

def run_program(args):
    tree = args.tree
    alignment = args.ref_msa
    query = args.query
    raxml_info_file = args.info
    maxPplacer = int(args.maxPplacer)
    nClades = int(args.numClades)
    numThreads = int(args.numThreads)
    timer = Timer()
    timer.tic("execution")
    timer.tic("setup")
    tmpdir = str(uuid.uuid4())
    os.mkdir(tmpdir)
    oldDir = os.getcwd()
    os.chdir(tmpdir)

    tree = f"{oldDir}/{tree}"
    alignment = f"{oldDir}/{alignment}"
    raxml_info_file = f"{oldDir}/{raxml_info_file}"

    backboneTree = read_tree(tree)

    applesPlacement = "apples.jplace"
    Tapples = f"result-{nClades}.tre"
    queryAlignment = "query.fa"
    generate_fasta_file(backboneTree, query, alignment, queryAlignment)
    timer.toc("setup")
    
    timer.tic("apples")
    run_apples(alignment, tree, queryAlignment, applesPlacement, numThreads)
    timer.toc("apples")
    place_sequence_in_subtree(applesPlacement, Tapples)

    scores = [-np.inf for i in range(nClades+1)] # last position is apples

    scores[-1] = run_raxml(Tapples, alignment)

    def execute_with_random_clade(item):
        threadIdx = item
        cladeNodes = sampleCompact(queryNode, maxPplacer)
        outputSubTree = f"clade-subtree-{threadIdx}.tre"
        taxa = [node.taxon.label for node in cladeNodes]
        subTree = backboneTree.extract_subtree_with_taxa_labels(taxa)
        subTree.write(file=open(outputSubTree, "w"), schema="newick")
        generate_fasta_file(subTree, query, alignment, queryAlignment)
        placementOutput = f"pplacer-{threadIdx}.jplace"
        run_pplacer(raxml_info_file, outputSubTree, queryAlignment, placementOutput)
        subTreeWithPlacement = f"subtree-with-placement-{threadIdx}.tre"
        place_sequence_in_subtree(placementOutput, subTreeWithPlacement)
        resultTree = read_tree(subTreeWithPlacement)
        newbackboneTree = read_tree(tree)
        modify_backbone_tree_with_placement(resultTree, newbackboneTree, query)
        temporaryBackBoneTree = f"result-{threadIdx}.tre"
        newbackboneTree.write(file=open(temporaryBackBoneTree, "w"), schema="newick")
        scores[threadIdx] = score_raxml(temporaryBackBoneTree, alignment)
        return
    # TODO: turn on threading
    for thread in range(nClades):
        execute_with_random_clade(thread)

    # Do maxLoc reduction to find best tree
    bestOne = None
    maxScore = -np.inf
    for i, score in enumerate(scores):
      if score > maxScore:
        maxScore = score
        bestOne = f"result-{i}.tre"
    # move best on to old dir, remove old dir
    shutil.move(bestOne, f"{oldDir}/{bestTree}")
    os.chdir(oldDir)
    shutil.rmtree(tmpdir)
    
    timer.toc("execution")
    timer.dump()

    return

def runWithTrainingData():
    train_data = "../../data/train/1000M1" # testing data includes R0..R19 (20 reps)
    print("Running APPLES+pplacer with: " + train_data)

    # for each of the sequences in each of the sets in the training data, perform a take-one-out run to see the placement
    # take-one-out will only show the accuracy, not how well it works with many sequences at once... need to 

    print("Done, but you haven't implemented anything yet")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Place a query sequence into a backbone tree.')
    # optional args
    parser.add_argument('-j', '--numThreads', default=1, help='Number of threads (default: 1)')
    parser.add_argument('-m', '--maxPplacer', default=500, help='Maximum size of subtree to hand to pplacer')
    parser.add_argument('-n', '--numClades', default=1, help='Number of (random) clades to consider')
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
