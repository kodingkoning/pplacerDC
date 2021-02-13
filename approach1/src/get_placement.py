#!/usr/bin/env python3
import json
import sys

filename = sys.argv[1]
f = open(filename,)
data = json.load(f)
print(data['placements'][0]['p'][0][1])

