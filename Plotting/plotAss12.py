import numpy as np
import matplotlib.pyplot as plt
import csv
from math import *

# orbit class
class orbit:
    def __init__(self , m , x , y , xp , yp):
        self.m = m
        self.x = x
        self.y = y
        self.xp = xp
        self.yp = yp


# we read the orbits
orbits = []
filename = "../cmake-build-debug/Ass12.out"
with open(filename, 'r') as csvfile:
    myreader = csv.reader(csvfile, delimiter=' ', quotechar='|', skipinitialspace=True)

    row = next(myreader)
    n = int(row[0]) # number of orbits

    for i in range(4*n):
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

        orbits.append(orbit(m,x,y,xp,yp))


# PLOT OF THE UNSTABLE MANIFOLDS, (x,y) projection
fig1 = plt.figure()
ax1 = fig1.add_subplot(1 , 1 , 1)
for i in range(n):
    if (i%6 == 0):
        ax1.plot(orbits[i].x , orbits[i].y , 'r')
        ax1.plot(orbits[n + i].x , orbits[n + i].y , 'r')

ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_title("Plot of unstable manifolds, (x,y) projection")


# PLOT OF THE STABLE MANIFOLDS, (x,y) projection
fig2 = plt.figure()
ax2 = fig2.add_subplot(1 , 1 , 1)
for i in range(n):
    if (i%6 == 0):
        ax2.plot(orbits[2*n + i].x , orbits[2*n + i].y , 'b')
        ax2.plot(orbits[3*n + i].x , orbits[3*n + i].y , 'b')

ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_title("Plot of stable manifolds, (x,y) projection")


# PLOT OF THE STABLE & UNSTABLE MANIFOLDS, (x,y) projection
fig3 = plt.figure()
ax3 = fig3.add_subplot(1 , 1 , 1)
for i in range(n):
    if (i%6 == 0):
        ax3.plot(orbits[i].x , orbits[i].y , 'r')
        ax3.plot(orbits[n + i].x , orbits[n + i].y , 'r')
        ax3.plot(orbits[2*n + i].x , orbits[2*n + i].y , 'b')
        ax3.plot(orbits[3*n + i].x , orbits[3*n + i].y , 'b')

ax3.set_xlabel("x")
ax3.set_ylabel("y")
ax3.set_title("Plot of unstable (red) and stable (blue) manifolds, (x,y) projection")


# PLOT OF THE STABLE & UNSTABLE MANIFOLDS, (x,y) projection
fig4 = plt.figure()
ax4 = fig4.add_subplot(1 , 1 , 1)
fig5 = plt.figure()
ax5 = fig5.add_subplot(1 , 1 , 1)

y = np.zeros(2*n)
xp = np.zeros(2*n)
for i in range(n):
    m = orbits[i].m
    y[i] = orbits[i].y[m-1]
    xp[i] = orbits[i].xp[m-1]

    m = orbits[n + i].m
    y[n + i] = orbits[n + i].y[m-1]
    xp[n + i] = orbits[n + i].xp[m-1]
ax4.plot(y[0:n] , xp[0:n] , 'r')
ax4.plot(y[n:2*n] , xp[n:2*n] , 'r')
ax5.plot(y[0:n] , xp[0:n] , 'r')

for i in range(n):
    m = orbits[2*n + i].m
    y[i] = orbits[2*n + i].y[m-1]
    xp[i] = orbits[2*n + i].xp[m-1]

    m = orbits[3*n + i].m
    y[n + i] = orbits[3*n + i].y[m-1]
    xp[n + i] = orbits[3*n + i].xp[m-1]
ax4.plot(y[0:n] , xp[0:n] , 'b')
ax4.plot(y[n:2*n] , xp[n:2*n] , 'b')

ax4.set_xlabel("y")
ax4.set_ylabel("x'")
ax4.set_title("Plot of unstable (red) and stable (blue) manifolds, (y,x') projection")

ax5.set_xlabel("y")
ax5.set_ylabel("x'")
ax5.set_title("Plot of half of the unstable manifolds, (y,x') projection")

plt.show()