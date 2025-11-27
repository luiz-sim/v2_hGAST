import os,math
import numpy as np

with open('truss.inp','r') as f:
	data = [i.split() for i in f]
L,dens= [float(data[i][3]) for i in range(71,119)],[float(data[i][6]) for i in range(71,119)]  #dens N/m


weight=[L[i]*dens[i] for i in range(5,15)]

for i in range(21,31):
    weight.append(L[i]*dens[i])

for i in range(37,len(dens)):
    weight.append(L[i]*dens[i])

weight=np.array(weight)
WEIGHT=np.sum(weight)
print(round(WEIGHT,2))
