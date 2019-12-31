import numpy as np
import matplotlib.pyplot as plt
import csv
from math import *

# orbit class
class orbit:
    def __init__(self , m , C , x , y , xp , yp):
        self.m = m
        self.C = C
        self.x = x
        self.y = y
        self.xp = xp
        self.yp = yp


# we read the orbits
orbits = []
filename = "../cmake-build-debug/Ass11.out"
with open(filename, 'r') as csvfile:
    myreader = csv.reader(csvfile, delimiter=' ', quotechar='|', skipinitialspace=True)

    row = next(myreader)
    n = int(row[0]) # number of orbits

    for i in range(n):
        row = next(myreader)
        C = float(row[0]) # Value of C for this orbit

        row = next(myreader)
        m = int(row[0]) # Number of points in the orbit

        t = np.zeros(m)
        x = np.zeros(m)
        y = np.zeros(m)
        xp = np.zeros(m)
        yp = np.zeros(m)

        for j in range(m):
            row = next(myreader)
            t[j] = float(row[0])
            x[j] = float(row[1])
            y[j] = float(row[2])
            xp[j] = float(row[3])
            yp[j] = float(row[4])

        orbits.append(orbit(m,C,x,y,xp,yp))


# PLOT OF (C,x)
C = np.zeros(n)
x = np.zeros(n)

for i in range(n):
    C[i] = orbits[i].C
    x[i] = orbits[i].x[0]

fig1 = plt.figure()
ax1 = fig1.add_subplot(1 , 1 , 1)
ax1.plot(C , x)
ax1.set_xlabel("C")
ax1.set_ylabel("x")

# PLOT OF SOME ORBITS
targets = [2.1 , 2.5 , 3.15]

fig2 = plt.figure()
ax2 = fig2.add_subplot(1 , 1 , 1)
for i in range(3):
    index = np.where(np.amin(np.abs(C - targets[i])) == np.abs(C - targets[i]))[0][0]
    ax2.plot(orbits[index].x , orbits[index].y)
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.legend(["C=2.1" , "C=2.5" , "C=3.15"])

plt.show()