import dendropy
import subtree_code

def runPplacerAPPLES(tree, alignment, query, maxPplacer):
    # TODO: will later need to find out how to do all of this for multiple queries at once

    # Step 1: run APPLES on backbone tree
    # TODO: understand running APPLES and pass query to the tree with the alignment
    apples_result = runAPPLES(tree, alignment, query)

    # Step 2: find the "region" of APPLES' placement
    # TODO: understand how to use dendropy to find the appropriate clade(s) to represent the region of the tree near where APPLES placed the query seq(s)
    tree_region = sampleCompact(apples_result, maxPplacer)

    # Step 3: Run pplacer on the subtree
    # TODO: understand how to use pplacer on that subtree (need to be able to pass that subtree in an appropriate format to pplacer and )
    pplacer_result = runPplacer(tree_region, alignment, query)

    return pplacer_result

def runPplacer(tree, alignment, query):
    # run for comparison with the combined version
    # pplacer -m GTR -s RAxML_info.REF -t backbone.nwk -o query.jplace aln_dna.fa -j 1
    return

def runAPPLES(tree, alignment, query):
    # run for comparison with the combined version
    # python3 ~/apples/run_apples.py -t backbone.nwk -s aln_dna.fa -q query.fa -T 16 -o apples.jplace
    return

def runWithTrainingData():
    train_data = "../../data/train/1000M1" # testing data includes R0..R19 (20 reps)
    print("Running APPLES+pplacer with: " + train_data)

    # for each of the sequences in each of the sets in the training data, perform a take-one-out run to see the placement
    # take-one-out will only show the accuracy, not how well it works with many sequences at once... need to 

    print("Done, but you haven't implemented anything yet")

if __name__ == "__main__":
    # Step 1: run APPLES on full backbone tree

    #TODO: take in pplacer and APPLES paths to where they are locally installed, if they are locally installed
    runWithTrainingData()
