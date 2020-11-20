import matplotlib.pyplot as plt
import numpy as np

# just take R0-treenwk.output and make a graph of the results from R0 and the tree.nwk vs the predicted placement

# ppa_results = open("R0-treenwk.output").readlines()
# app_results = open("R0-apples-output.txt").readlines()
# pp_results  = open("R0-pp-output.txt").readlines()

# get list of the error (FN, FP, RF) over all of the placements

# output from treecompare.py is: nl, ei1, ei2, fp, fn, rf
def get_error(file_name):
    FN = []
    FP = []
    RF = []

    for line in open(file_name).readlines():
        if line.startswith("("):
            line_list = line.split(",")
            FP.append(int(line_list[3])/int(line_list[2]))
            FN.append(int(line_list[4])/int(line_list[1]))
            RF.append(float(line_list[5][0:-2]))
    return FN, FP, RF

ppa_FN, ppa_FP, ppa_RF = get_error("R0-treenwk.output")
app_FN, app_FP, app_RF = get_error("R0-apples-output.txt")
pp_FN,  pp_FP,  pp_RF  = get_error("R0-pp-output-results.txt")

# make plot of FN rate
fig, ax = plt.subplots()
ax.boxplot([ppa_FN, app_FN, pp_FN], labels=["Combined", "APPLES", "pplacer"], positions=[0,1,2])
plt.title("FN rates for R0 of 1000M1")
plt.tight_layout()
plt.show()

# make plot of FP rate
fig, ax = plt.subplots()
ax.boxplot([ppa_FP, app_FP, pp_FP], labels=["Combined", "APPLES", "pplacer"], positions=[0,1,2])
plt.title("FP rates for R0 of 1000M1")
plt.tight_layout()
plt.show()

# make plot of RF rate
fig, ax = plt.subplots()
ax.boxplot([ppa_RF, app_RF, pp_RF], labels=["Combined", "APPLES", "pplacer"], positions=[0,1,2])
plt.title("RF rates for R0 of 1000M1")
plt.tight_layout()
plt.show()


# make these box plots
# fig, ax = plt.subplots()
# ax.boxplot(FP, labels=["FP"], positions=[0])
# ax.boxplot(FN, labels=["FN"], positions=[1])
# ax.boxplot(RF, labels=["RF"], positions=[2])
# plt.show()
# plt.tight_layout()
