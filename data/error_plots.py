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
def read_list(fileName):
    with open(fileName) as fileHandle:
        content = fileHandle.readlines()
    content=[float(x.strip()) for x in content]
    return np.array(content)
def get_errs(method, size):
  return np.sort(read_list(f"vs-data/{size}/delta_error_{method}.txt"))
def get_errs(method, size):
  return np.sort(read_list(f"vs-data/{size}/delta_error_{method}.txt"))

sizes = [500,1000,5000,10000]
pplacer_sizes = [500,1000]

def make_figure(size):
  data=[]
  methods=[]
  if size in pplacer_sizes:
    methods=["DC-pplacer", "pplacer", "pplacer+APPLES*", "APPLES*"]
  else:
    methods=["DC-pplacer", "pplacer+APPLES*", "APPLES*"]
  data.append(get_errs("approach1", size))
  if size in pplacer_sizes:
    data.append(get_errs("approach1", size))
  data.append(get_errs("approach2", size))
  data.append(get_errs("apples", size))
  fig, ax = plt.subplots()
  ax.boxplot(data)
  ax.set_xticklabels(methods)
  plt.title(f"Delta Error for Various Methods, VS-{size}")
  plt.ylabel("Delta Error")
  plt.xticks(rotation=45)
  plt.tight_layout()
  plt.savefig(f"../writeup/Figs/VS-delta-error-{size}.png",dpi=150)

for size in sizes:
  make_figure(size)

def make_everything():
  data=[[],[],[],[]]
  methods=["DC-pplacer", "pplacer", "pplacer+APPLES*", "APPLES*"]
  for size in sizes:
    data[0].extend(get_errs("approach1", size))
    if size in pplacer_sizes:
      data[1].extend(get_errs("approach1", size))
    data[2].extend(get_errs("approach2", size))
    data[3].extend(get_errs("apples", size))
  fig, ax = plt.subplots()
  ax.boxplot(data)
  ax.set_xticklabels(methods)
  plt.title(f"Delta Error for Various Methods, All sizes")
  plt.ylabel("Delta Error")
  plt.xticks(rotation=45)
  plt.tight_layout()
  plt.savefig(f"../writeup/Figs/VS-delta-error-all-sizes.png",dpi=150)
make_everything()


if 0:
  plt.step(data[0], np.arange(data[0].size), label="DQ-pplacer")  # From 0 to the number of data points-1
  plt.step(data[1], np.arange(data[1].size), label="pplacer")  # From 0 to the number of data points-1
  plt.legend()
  plt.show()
