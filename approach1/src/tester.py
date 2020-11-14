from tree_decomposition import *
from script_executor import *


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
  parentNodeName = ""
  queryNode = None
  branchHit = False
  parentNode = None
  # Likely a much better way to do this...
  print(f"Hunting for query sequence {querySequence}")
  for leaf in treeWithPlacement.leaf_nodes():
    print(f"On leaf {leaf}")
    if leaf.taxon.label == querySequence:
      print(leaf.adjacent_nodes())
      #parentNode, _ = leaf.adjacent_nodes()
      parentNode = leaf.adjacent_nodes()[0]
      #parentNodeName = leaf.adjacent_nodes()[0].taxon.label
      queryNode = leaf
      branchHit = True
      break
  assert branchHit, "Expected to actually find the parent node!"
  
  branchHit = False
  # Clone larger tree, add in sequence
  #backBoneTreeWithPlacement = read_tree(open(prunedTree, "r")).get_tree().clone()
  print(f"Parent node = {parentNode}")
  for node in tree.get_tree().nodes():
    print(f"Node = {node}")
    if node == parentNode:
      node.add_child(queryNode)
      branchHit = True
  assert branchHit, "Expected to actually find the parent node in the backbone tree!"

  tree.update_bipartitions()

  # Dump out tree
  tree.write(file=open("the-result.tre", "w"), schema="newick")


