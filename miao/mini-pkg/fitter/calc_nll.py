import numpy as np
import matplotlib.pyplot as plt

nevt = 304

nuTime, deltaNuTime, prob = [], [], []
tmp_prob = 0
with open("log") as f:
    for lines in f.readlines():
        line = lines.strip("\n")
        data = line.split(" ")
        if int(data[0]) == 0:
            print(tmp_prob)
            tmp_prob = 0
        if int(data[0]) < nevt:
            tmp_prob += -np.log(float(data[4]))




