import os, math
import numpy as np

with open('./LOADS_aer.dat','r') as f:
    data = [i.split() for i in f]

cp = [float(data[i][3]) for i in range(0,len(data))]
lamda = [float(data[i][4]) for i in range(0,len(data))]

cpmax = max(cp)
lopt = lamda[cp.index(max(cp))]

print(cpmax)
print(lopt)

R = 121.119
dens = 1.225
omega = 0.7854


M = cpmax*0.5*dens*math.pi*R**5 * omega**2/lopt**3
print(M)
