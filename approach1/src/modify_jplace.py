#!/usr/bin/env python3
import json
import sys

filenameIn  = sys.argv[1]
filenameOut = sys.argv[2]
newIndex = int(sys.argv[3])
newTreeFile = sys.argv[4]
with open(filenameIn,"r") as f:
    data = json.load(f)
placement_p = data['placements'][0]['p'][0]
data['placements'][0]['p'] = [placement_p]
data['placements'][0]['p'][0][1] = newIndex

with open(newTreeFile,"r") as f:
    data['tree'] = f.readlines()[0]

with open(filenameOut,"w") as f:
    json.dump(data,f,indent=1)

