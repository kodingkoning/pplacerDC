import dendropy
import argparse

def compareTreesFromPath(treePath1, treePath2):
    #print("Comparing {} with {}".format(treePath1, treePath2))
    
    tax = dendropy.TaxonNamespace()
    tr1 = dendropy.Tree.get(path=treePath1,
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax,
                            preserve_underscores=True)
    tr2 = dendropy.Tree.get(path=treePath2,
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax,
                            preserve_underscores=True)
    
    tr1.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    tr2.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    
    return compareDendropyTrees(tr1, tr2)
    #print("RF distance on %d shared leaves: %d" % (nl, fp + fn))
    
def compareDendropyTrees(tr1, tr2):
    from dendropy.calculate.treecompare \
        import false_positives_and_negatives
    
    lb1 = set([l.taxon.label for l in tr1.leaf_nodes()])
    lb2 = set([l.taxon.label for l in tr2.leaf_nodes()])
    
    com = lb1.intersection(lb2)
    if com != lb1 or com != lb2:
        com = list(com)
        tns = dendropy.TaxonNamespace(com)
        
        tr1.retain_taxa_with_labels(com)
        tr1.migrate_taxon_namespace(tns)
        
        tr2.retain_taxa_with_labels(com)
        tr2.migrate_taxon_namespace(tns)
    com = list(com)
    
    tr1.update_bipartitions()
    tr2.update_bipartitions()
    
    nl = len(com)
    ei1 = len(tr1.internal_edges(exclude_seed_edge=True))
    ei2 = len(tr2.internal_edges(exclude_seed_edge=True))
    
    [fp, fn] = false_positives_and_negatives(tr1, tr2)
    rf = float(fp + fn) / (ei1 + ei2)
    
    return (nl, ei1, ei2, fp, fn, rf)
def delta_error(t1,t2, trueTree):
    nl, ei1, ei2, fp, fn1, rf = compareTreesFromPath(trueTree, t1)
    nl, ei1, ei2, fp, fn2, rf = compareTreesFromPath(trueTree, t2)
    return fn1 - fn2

def print_all_errors(t1,t2,trueTree):
    nl_t1, ei1_t1, ei2_t1, fp_t1, fn_t1, rf_t1 = compareTreesFromPath(trueTree, t1)
    nl, ei1, ei2, fp, fn2, rf = compareTreesFromPath(trueTree, t2)
    # print delta error, fp rate, fn rate, and rf rate. Focus on t1, the output tree
    # For binary trees, fn = fp = rf
    #print(abs(fn_t1-fn2), rf_t1)
    fp_rate = fp_t1/ei2_t1
    fn_rate = fn_t1/ei1_t1
    print(abs(fn_t1-fn2), fp_rate, fn_rate, rf_t1)

if __name__ == "__main__":
  #raxml = "RAxML_result.REF"
  #true = "true_topo.tree"
  #approach1Tree = "approach1.tre"
  #approach2Tree = "approach2.tre"
  #applesTree = "apples.tre"
  #de_apples = delta_error(applesTree,raxml, true)
  #de1 = delta_error(approach1Tree,raxml, true)
  #de2 = delta_error(approach2Tree,raxml, true)
  #print(f"Delta error apples = {de_apples}")
  #print(f"Delta error approach1 = {de1}")
  #print(f"Delta error approach2 = {de2}")
  parser = argparse.ArgumentParser(description='Evaluate the delta error.')
  parser.add_argument('trueTreeFile', type=str,
                     help='True backbone tree')
  parser.add_argument('outputTreeFile', type=str,
                     help='Output tree (e.g., from pplacer)')
  parser.add_argument('estimatedTreeFile', type=str,
                     help='Input tree (e.g., RAxML) with the correct placement')
  args = parser.parse_args()
  trueTreeFile = args.trueTreeFile
  outputTreeFile = args.outputTreeFile
  estimatedTreeFile = args.estimatedTreeFile
  print_all_errors(outputTreeFile, estimatedTreeFile, trueTreeFile)
  #de = delta_error(outputTreeFile, estimatedTreeFile, trueTreeFile)
  #print(abs(de))

