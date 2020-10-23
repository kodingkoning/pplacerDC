from sepp_utils.alignment import MutableAlignment
from sepp_utils.tree import PhylogeneticTree
import dendropy

def read_alignment_and_tree(alignment_file, tree_file):
    alignment = MutableAlignment()
    alignment.read_file_object(alignment_file)
    tree = PhylogeneticTree(
      dendropy.Tree.get_from_stream(tree_file,
                                    schema="newick",
                                    preserve_underscores=True))
    return (alignment, tree)

if __name__ == "__main__":
  base_dir = "/home/malachi/work/classes/variable-size/data/500/0/"
  tree_size = 500
  tree_file = open(f"{base_dir}/RAxML_result.REF8", "r")
  alignment_file = open(f"{base_dir}/aln_dna.fa", "r")
  alignment, tree = read_alignment_and_tree(alignment_file,
       tree_file)

  maxSize = 100 # e.g.
  decomposed_tree = tree.decompose_tree(maxSize, strategy="centroid", minSize=1)
  assert len(decomposed_tree.keys())*maxSize >= tree_size, "Expected at least 5 sub trees in the decomposition!"
  print(decomposed_tree.keys())
