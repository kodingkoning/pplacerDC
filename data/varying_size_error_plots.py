import numpy as np
import re
import matplotlib.pyplot as plt
def read_list(fileName):
    with open(fileName) as fileHandle:
        content = fileHandle.readlines()
    content=[float(x.strip()) for x in content]
    return np.array(content)
def get_errs(method, size):
  return np.sort(read_list(f"vs-data/{size}/delta_error_{method}.txt"))

problemSize = 5000
testedSizes = [100, 500, 1000]
data=[]
for size in testedSizes:
  data.append(get_errs(f"approach1_size_{size}", problemSize))
data.append(np.sort(read_list("BW-results/5000/delta_error_approach1.txt")))
# TODO: append data for 2500
fig, ax = plt.subplots()
ax.boxplot(data)
labels = [100, 500, 1000, 2500]
ax.set_xticklabels(labels)
plt.title(f"Delta Error for DC-pplacer v. Subproblem Size")
plt.ylabel("Delta Error")
plt.xlabel("Subproblem Size")
plt.xticks(rotation=45)
plt.tight_layout()
#plt.show()
plt.savefig(f"../writeup/Figs/varying-subproblem-size.png",dpi=150)
