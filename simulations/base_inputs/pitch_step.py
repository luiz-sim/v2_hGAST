import os, math
import numpy as np

with open('./qsdof.dat','r') as f:
    data = [i.split() for i in f]

time = [float(data[i][0]) for i in range(0,len(data))]
pitch = [float(data[i][4]) for i in range(0,len(data))]


PITCH = np.average(pitch)
print(PITCH)
