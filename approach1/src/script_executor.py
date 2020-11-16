import subprocess
import os
import tempfile
import uuid
import re
import glob
def field_by_regex(regex,log_file_name, fieldnum = 0):
    with open(log_file_name) as log_file:
        content =log_file.readlines()
    content = [x.strip() for x in content]
    field=[]
    for line in content:
        m = re.search(regex,line)
        if m:
            field.append(float(m.groups()[fieldnum]))
    return field
def generate_fasta_file(subtree, querySequence, referenceFastaFile, outputReferenceFile, debugOutput=False):
    concatSequences = ""
    leaves = subtree.leaf_node_names()
    for leaf in leaves:
      concatSequences += f" {leaf}"
    concatSequences += f" {querySequence}"
    command = "python3 faSomeRecords.py"
    command += " --records " + concatSequences
    command += " --fasta " + referenceFastaFile
    command += " --outfile " + outputReferenceFile
    if debugOutput:
        print(f"Command =\n{command}")
    subprocess.call([command], shell=True)

def place_sequence_in_subtree(pplacerOutput, outputTreeFileName):
    subprocess.call(["guppy","tog",
        "-o", outputTreeFileName,
        pplacerOutput
        ])

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
def score_raxml(treeFile, referenceAln):
    """
    Run RAxML-ng in fixed tree mode to determine the best LogLikelihood score possible
    treeFile: the file where the tree is located
    referenceAln: the file where the reference alignment is located
    """
    tmpFile = str(uuid.uuid4())
    tmpFileHandle = open(tmpFile,"w")
    # generate temporary file handle
    subprocess.call(["raxml-ng",
                     "--msa", referenceAln,
                     "--model", "GTR+G",
                     "--tree", treeFile,
                     "--threads", "1", # run in serial
                     "--opt-branches", "off", # do not optimize branch lengths
                     "--opt-model", "off", # do not optimize model conditions
                     "--evaluate"   # fixed-tree evaluation
                     ],
                     stdout=tmpFileHandle
                     )
    # parse the output for the score
    regex = "Final LogLikelihood: (.+)"
    score = field_by_regex(regex, tmpFile)[0]
    # delete temporary file
    tmpFileHandle.close()
    os.remove(tmpFile)

    # delete raxml files
    for fl in glob.glob(f"{referenceAln}.raxml*"):
        os.remove(fl)
    return score
