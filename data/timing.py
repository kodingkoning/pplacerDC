import numpy as np
import re
import matplotlib.pyplot as plt
def field_by_regex(regex,log_file_name, fieldnum = 0):
    with open(log_file_name) as log_file:
        content =log_file.readlines()
    content = [x.strip() for x in content]
    field=[]
    for line in content:
        m = re.search(regex,line)
        if m:
            field.append(float(m.groups()[fieldnum]))
    return np.array(field)

def field_by_name(log_file_name):
    with open(log_file_name) as log_file:
        content =log_file.readlines()
    content = [x.strip() for x in content]
    field=[]
    for line in content:
        field.append(float(line))
    return np.array(field)

def read_list(fileName):
    with open(fileName) as fileHandle:
        content = fileHandle.readlines()
    content=[float(x.strip()) for x in content]
    return np.array(content)

sizes = np.array([500, 1000, 5000, 10000, 50000])
def get_time_by_regex(filename, regex, sizes):
  avg_times = np.zeros(len(sizes))
  std_times = np.zeros(len(sizes))
  for i, size in enumerate(sizes):
    outputFile = f"RNASim-VS-results/variable-size/data/{size}/{filename}"
    times = field_by_regex(regex, outputFile)
    avg_times[i] = np.average(times)
    std_times[i] = np.std(times)
  return avg_times, std_times

def get_time_by_list(filename,regex, sizes):
  avg_times = np.zeros(len(sizes))
  std_times = np.zeros(len(sizes))
  for i, size in enumerate(sizes):
    outputFile = f"test/RNASim-VS-results/variable-size/data/{size}/{filename}"
    times = field_by_name(outputFile)
    avg_times[i] = np.average(times)
    std_times[i] = np.std(times)
  return avg_times, std_times


regex = "Program execution took (.+?) s"
avg_times, std_times = get_time_by_list("time_approach1.txt",regex,sizes)
fig = plt.figure()
ax = plt.axes()
ax.set_xscale("log")
ax.set_yscale("log")
ax.errorbar(sizes,avg_times,yerr=std_times, label="pplacerDC", marker="o")
ax.set_xlabel("Number of Taxa")
ax.set_ylabel("Time (s)")
ax.set_title("Time v. Number of Taxa")

regex = "execution took (.+?) s"
avg_times, std_times = get_time_by_list("time_approach2.txt", regex, sizes)
ax.errorbar(sizes,avg_times,yerr=std_times, label="pplacerAPPLES*", marker="o")

# Grab some pplacer timing data on 500, 1000 datasets
pplacer_regex = "Running pplacer... took (.+?) s"
pplacer_sizes = np.array([500,1000])
avg_t_pplacer, std_t_pplacer = get_time_by_list("time_pplacer.txt", pplacer_regex, pplacer_sizes)
ax.errorbar(pplacer_sizes,avg_t_pplacer,yerr=std_t_pplacer, label="pplacer", marker="x")

# Grab some pplacer timing data on 500, 1000 datasets
apples_regex = "apples took (.*) s"
avg_t_apples, std_t_apples = get_time_by_list("time_apples.txt", apples_regex, sizes)
ax.errorbar(sizes,avg_t_apples,yerr=std_t_apples, label="APPLES*", marker="x")

plt.xticks([500,1000,5000,10000,50000])
ax.set_xticklabels([500,1000,5000,10000,50000])

ax.legend()
#plt.show()

plt.tight_layout()
plt.savefig("../writeup/Figs/VS-timing-results-BW.png", dpi=150)
