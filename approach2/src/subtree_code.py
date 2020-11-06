# Code from Vlad Smirnov

def sampleCompact(leaf, size):
    taxons = [leaf.taxon.label]
    node = leaf
    while len(taxons) < size:
        for bro in node.sibling_nodes():
            taxons.extend(collectSubtreeTaxa(bro, size - len(taxons)))
            if len(taxons) >= size:
                return taxons
        
        node = node.parent_node      
    return taxons

def collectSubtreeTaxa(node, numTaxa):
    if node.is_leaf():
        return [node.taxon.label]
    
    taxons = []
    for child in node.child_nodes():
        taxons.extend(collectSubtreeTaxa(child, numTaxa - len(taxons)))
        if len(taxons) >= numTaxa:
            return taxons
    return taxons