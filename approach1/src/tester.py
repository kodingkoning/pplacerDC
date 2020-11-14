from tree_decomposition import *
from script_executor import *


########### TODO ##############
# - Clean up this code

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

maxSize = 5 # e.g., this will be a tuneable parameter
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
  query_alignment_file = "foobar.fa" # test <-----
  generate_fasta_file(tree_object, querySequence, alignment_file, query_alignment_file)
  run_pplacer(raxml_info_file, outputTreeFile, alignment_file, query_alignment_file, outputLocation)
  # overwrite the current file
  resultTree = "mytestoutput.tre"
  place_sequence_in_subtree(outputLocation, "mytestoutput.tre")

  # Find leaf node matching query sequence
  treeWithPlacement = read_tree(open(resultTree, "r")).get_tree()

  # Reload original tree, since it seems the initial
  # tree has been modified by the centroid decomposition
  theTree = read_tree(open(prunedTree, "r")).get_tree()

  subTreeNodeToNode={}
  label_internal_nodes_subtree(treeWithPlacement, theTree, subTreeNodeToNode)
  #queryNode = treeWithPlacement.find_node_with_label(querySequence)
  #assert not queryNode, "Expected queryNode to not be None type!"
  queryNode = None
  for leaf in treeWithPlacement.leaf_nodes():
    if leaf.taxon.label == querySequence:
      queryNode = leaf
      break
  if queryNode == None:
    print("Expected non null queryNode")


  parentSubTree = queryNode.adjacent_nodes()[0]
  print(f"Parent sub tree = {parentSubTree}")
  for key in subTreeNodeToNode.keys():
    print(f"Key = {key}")
  parent = subTreeNodeToNode[parentSubTree]
  parent.add_child(queryNode)

  #theTree.update_bipartitions()

  # Dump out tree
  theTree.write(file=open("the-result.tre", "w"), schema="newick")


