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
    command = "faSomeRecords.py"
    command += " --records " + concatSequences
    command += " --fasta " + referenceFastaFile
    command += " --outfile " + outputReferenceFile
    if debugOutput:
        print(f"Command =\n{command}")
    tmpFile = str(uuid.uuid4())
    tmpFileHandle = open(tmpFile,"w")
    ret = subprocess.call([command], shell=True, stdout=tmpFileHandle, stderr=tmpFileHandle)
    if ret != 0:
        faSomeRecords_error = tmpFileHandle.readlines()
        print(f"Failed to run fasomeRecords.py, it said:\n{faSomeRecords_error}")
    tmpFileHandle.close()
    if ret != 0:
        exit(-1)


def place_sequence_in_subtree(pplacerOutput, outputTreeFileName):
    # Ubuntu 18.0.4 workaround for a bad assertion in loadLocale.c:129
    # THIS IS A BODGE!!!
    # see: https://askubuntu.com/questions/1081901/what-is-the-correct-way-to-fix-an-assertion-in-loadlocale-c
    # for more details
    os.environ["LC_ALL"] = "C"
    subprocess.call(["guppy","tog",
        "-o", outputTreeFileName,
        pplacerOutput
        ])

def run_pplacer(raxml_info_file, backbone_tree, queries, output):
    # Ubuntu 18.0.4 workaround for a bad assertion in loadLocale.c:129
    # THIS IS A BODGE!!!
    # see: https://askubuntu.com/questions/1081901/what-is-the-correct-way-to-fix-an-assertion-in-loadlocale-c
    # for more details
    os.environ["LC_ALL"] = "C"
    tmpFile = str(uuid.uuid4())
    tmpFileHandle = open(tmpFile,"w+")
    ret = subprocess.call(["pplacer",
                       "-m", "GTR", # model
                       "-s", raxml_info_file, # raxml info file location
                       "-t", backbone_tree, # backbone tree
                       "-o", f"{output}", # output location
                       "-j", "1", # Run on single thread
                       queries    # Name of file containing query sequences
                       ],
                       stdout=tmpFileHandle,
                       stderr=tmpFileHandle
                       )
    if ret != 0:
        pplacer_error = tmpFileHandle.readlines()
        print(f"Failed to run pplacer, it said:\n{pplacer_error}")

    tmpFileHandle.close()
    if ret != 0:
      exit(-1)

def score_raxml(treeFile, referenceAln, tid, trial=0):
    """
    Run RAxML-ng in fixed tree mode to determine the best LogLikelihood score possible
    treeFile: the file where the tree is located
    referenceAln: the file where the reference alignment is located
    """
    maxTrials = 2
    if trial > maxTrials:
      # I give up.
      return
    #random_prefix = str(uuid.uuid4())
    prefix = f"raxml-prefix-{tid}.score"
    tmpFileHandle = open(prefix,"w")
    # generate temporary file handle
    ret = subprocess.call(["raxml-ng",
                     "--msa", referenceAln,
                     "--model", "GTR+G",
                     "--tree", treeFile,
                     #"--threads", "1", # run in serial
                     "--opt-branches", "off", # do not optimize branch lengths
                     "--opt-model", "off", # do not optimize model conditions
                     "--evaluate",   # fixed-tree evaluation
                     "--nofiles", # do not generate files
                     "--log", "result", # do not log anything
                     ],
                     stdout=tmpFileHandle
                     )
    # parse the output for the score
    tmpFileHandle.close()
    regex = "Final LogLikelihood: (.+)"
    try:
      score = field_by_regex(regex, prefix)[0]
    except:
      # there are spurious IO errors, so we can just re-run this
      trial += 1
      score_raxml(treeFile, referenceAln, tid, trial)
    if ret != 0:
      raxml_error = tmpFileHandle.readlines()
      print(f"Failed to run raxml-ng, it said:\n{raxml_error}")

    if ret != 0:
      exit(-1)
    return
