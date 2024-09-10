#!/usr/bin/env python

import sys
import json

f = open(sys.argv[1])
map = json.load(f)
f.close()

outfile = open("pu_file-new.txt", "w")

for run, lumis in map.items():
   for lumi in lumis:
				print(run, " ".join([str(num) for num in lumi]))
				outfile.write(run+" ".join([str(num) for num in lumi]))

outfile.close()
