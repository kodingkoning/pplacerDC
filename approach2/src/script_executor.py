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

def generate_fasta_file_apples(subtree, querySequence, referenceFastaFile, outputReferenceFile, debugOutput=False):
    concatSequences = ""
    #leaves = [node.taxon.label for node in subtree.leaf_nodes()]
    #for leaf in leaves:
    #  concatSequences += f" {leaf}"
    concatSequences += f" {querySequence}"
    command = "faSomeRecords.py"
    command += " --records " + concatSequences
    command += " --fasta " + referenceFastaFile
    command += " --outfile " + outputReferenceFile
    if debugOutput:
        print(f"Command =\n{command}")
    tmpFile = str(uuid.uuid4())
    tmpFileHandle = open(tmpFile,"w+")
    ret = subprocess.call([command], shell=True, stdout=tmpFileHandle)
    if ret != 0:
        faSomeRecords_error = tmpFileHandle.readlines()
        print(f"Failed to run fasomeRecords.py, it said:\n{faSomeRecords_error}")
    tmpFileHandle.close()
    os.remove(tmpFile)
    if ret != 0:
        exit(-1)

def generate_fasta_file(subtree, querySequence, referenceFastaFile, outputReferenceFile, debugOutput=False):
    concatSequences = ""
    leaves = [node.taxon.label for node in subtree.leaf_nodes()]
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
    ret = subprocess.call([command], shell=True, stdout=tmpFileHandle)
    if ret != 0:
        faSomeRecords_error = tmpFileHandle.readlines()
        print(f"Failed to run fasomeRecords.py, it said:\n{faSomeRecords_error}")
    tmpFileHandle.close()
    os.remove(tmpFile)
    if ret != 0:
        exit(-1)


def place_sequence_in_subtree(pplacerOutput, outputTreeFileName):
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
    tmpFileHandle = open(tmpFile,"w")
    ret = subprocess.call(["pplacer",
                       "-m", "GTR", # model
                       "-s", raxml_info_file, # raxml info file location
                       "-t", backbone_tree, # backbone tree
                       "-o", f"{output}", # output location
                       "-j", "1", # Run on single thread
                       queries    # Name of file containing query sequences
                       ],
                       stdout=tmpFileHandle
                       )
    if ret != 0:
        pplacer_error = tmpFileHandle.readlines()
        print(f"Failed to run pplacer, it said:\n{pplacer_error}")

    tmpFileHandle.close()
    os.remove(tmpFile)
    if ret != 0:
      exit(-1)
def run_apples(alignmentFile, backBoneTreeFile, querySequenceAlignmentFile, output, threads):
    tmpFile = str(uuid.uuid4())
    tmpFileHandle = open(tmpFile,"w")
    ret = subprocess.call(["run_apples.py",
                       "-s", alignmentFile, # reference alignment file
                       "-t", backBoneTreeFile, # backbone tree
                       "-o", f"{output}", # output location
                       "--threads", str(threads), # num threads to run
                       "-q", querySequenceAlignmentFile,
                       ],
                       #shell=True,
                       stdout=tmpFileHandle
                       )
    tmpFileHandle.close()
    if ret != 0:
        reader = open(tmpFile, "r")
        apples_error = reader.readlines()
        print(f"Failed to run apples, it said:\n{apples_error}")

    os.remove(tmpFile)
    if ret != 0:
      exit(-1)

def score_raxml(treeFile, referenceAln):
    """
    Run RAxML-ng in fixed tree mode to determine the best LogLikelihood score possible
    treeFile: the file where the tree is located
    referenceAln: the file where the reference alignment is located
    """
    tmpFile = str(uuid.uuid4())
    tmpFileHandle = open(tmpFile,"w")
    random_prefix = str(uuid.uuid4())
    # generate temporary file handle
    ret = subprocess.call(["raxml-ng",
                     "--msa", referenceAln,
                     "--model", "GTR+G",
                     "--tree", treeFile,
                     "--threads", "1", # run in serial
                     "--opt-branches", "off", # do not optimize branch lengths
                     "--opt-model", "off", # do not optimize model conditions
                     "--prefix", random_prefix, # do not optimize model conditions
                     "--evaluate"   # fixed-tree evaluation
                     ],
                     stdout=tmpFileHandle,
                     stderr=tmpFileHandle
                     )
    # parse the output for the score
    regex = "Final LogLikelihood: (.+)"
    # fields = field_by_regex(regex, tmpFile)
    # score = 0
    # if len(fields):
    #     score = fields[0]
    score = field_by_regex(regex, tmpFile)[0]
    if ret != 0:
      print(ret)
      raxml_error = tmpFileHandle.readlines()
      print(f"Failed to run raxml-ng, it said:\n{raxml_error}")

    # delete temporary file
    tmpFileHandle.close()
    os.remove(tmpFile)

    if ret != 0:
      exit(-1)
    return score
