from sepp_utils.alignment import MutableAlignment
from sepp_utils.tree import PhylogeneticTree
import dendropy
import subprocess
import os # TODO: may not be portable

def generate_fasta_file(subtree, querySequence):
    concatSequences = ""
    leaves = subtree.leaf_node_names()
    for leaf in leaves:
      concatSequences.join(" " + leaf)

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
def prune_query(referenceTree, querySequence, outputTree):
    fileHandle = open(outputTree, "w")
    p = subprocess.call(["nw_prune", referenceTree, querySequence], stdout=fileHandle)
def run_pplacer(raxml_info_file, backbone_tree, reference_aln, queries, output):
    # Ubuntu 18.0.4 workaround for a bad assertion in loadLocale.c:129
    # THIS IS A BODGE!!!
    # see: https://askubuntu.com/questions/1081901/what-is-the-correct-way-to-fix-an-assertion-in-loadlocale-c
    # for more details
    os.environ["LC_ALL"] = "C"
    subprocess.call(["pplacer",
                       "-m", "GTR", # model
                       "-s", raxml_info_file, # raxml info file location
                       "-t", backbone_tree, # backbone tree
                       "-o", f"{output}", # output location
                       #"-r", reference_aln, # reference alignment
                       "-j", "1", # Run on single thread
                       queries    # Name of file containing query sequences
                       ]
                       )

if __name__ == "__main__":
  os.chdir("test")
  base_dir = "/home/malachi/work/classes/variable-size/data/500/0/"
  tree_size = 500
  raxml_info_file = f"{base_dir}/RAxML_result.REF"
  querySequences = read_list("queries.txt")
  assert len(querySequences) == 1, "Expected only a single query sequence!"
  querySequence = querySequences[0]
  # generate backbone tree by removing query sequence
  prunedTree = "pruned-tree.tre"
  prune_query(raxml_info_file, querySequence, prunedTree)

  tree_file_handle = open(prunedTree, "r")
  alignment_file = f"{base_dir}/aln_dna.fa"
  alignment_file_handle = open(f"{base_dir}/aln_dna.fa", "r")
  alignment, tree = read_alignment_and_tree(alignment_file_handle,
       tree_file_handle)

  maxSize = 100 # e.g., this will be a tuneable parameter
  decomposed_trees = tree.decompose_tree(maxSize, strategy="centroid", minSize=1)
  assert len(decomposed_trees.keys())*maxSize >= tree_size, "Expected at least 5 sub trees in the decomposition!"
  for i, tree_key in enumerate(decomposed_trees.keys()): # For testing!
    if i > 0:
      continue # only run on first tree...
    tree_object = decomposed_trees[tree_key]
    numNodes = tree_object.count_nodes()
    numLeafs = tree_object.count_leaves()
    print(f"Tree ({tree_key}) has {numNodes} nodes and {numLeafs} leafs!")
    # output to a file through dendropy.Tree object
    outputTreeFile = f"decomposed-tree-{tree_key}.tre"
    tree_object.get_tree().write(file=open(outputTreeFile, 'w'),
        schema="newick")
    #subprocess.call(["cd", "test"])
    print(os.getcwd())
    outputLocation = "output.jplace"
    raxml_info_file = f"{base_dir}/RAxML_info.REF"
    query_alignment_file = "aln_dna.fa" # test <-----
    run_pplacer(raxml_info_file, outputTreeFile, alignment_file, query_alignment_file, outputLocation)
    #run_pplacer(raxml_info_file, outputTreeFile, alignment_file, alignment_file, outputLocation)


