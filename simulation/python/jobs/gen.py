import os
import sys

N_Emax = 3
N_time = 6
N = N_Emax * N_time

c = 0
Emin = 0.10
for i in range(N_Emax):
    Emax = 3 + i
    for j in range(N_time):
        tmin = -0.03 + 0.01 * j
        tmax = tmin + 0.01

        filename = f"./job{c}.sh"

        with open(filename, "w") as f:
            f.write("#!/bin/bash")
            f.write("\n")
            f.write(f"python /junofs/users/miaoyu/supernova/simulation/main.py {tmin} {tmax} {Emax}")


        os.system(f"chmod +x {filename}")
        os.system(f"hep_sub {filename}")

        c += 1
