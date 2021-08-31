import namespace as ns
import numpy as np
import os

mod = 82503
dist = 10
#nuMass = 0.1
nuType = -1
chaname = 1
path = "/junofs/users/miaoyu/supernova/wenlj"
dataMass = 0.0
dataMH = 2
group = 999

for num in range(20, 30, 1):
    nuMass = num/10.
    print("Submit jobs from nuMass : %.1f" %nuMass)
    for imh in range(1, 3, 1):
        print("mass hierarchy: %d" %imh)
        #for no in range(17, len(Evismin), 1):
        shname = ns.fitGenShName(mod, chaname, nuType, dataMass, dataMH, nuMass, imh, dist, group)
        with open( shname, "w") as f:
            f.write("#!/bin/bash")
            f.write("\n")
            f.write("cd %s" %path)
            f.write("\n")
            f.write("python nllFitnew.py %d %d %d %.1f %d %.1f %d %.1f %.1f %d &> %s" %(mod, chaname, nuType, dataMass, dataMH, nuMass, imh, dist, 0.1, group, ns.fitGenLogName(mod, chaname, nuType, dataMass, dataMH, nuMass, imh, dist, group)))
    
            os.system("chmod +x %s" %shname)
            os.system("hep_sub %s" %shname)

 
