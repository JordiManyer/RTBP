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


def readOrbits(filename):
    orbits = []
    with open(filename, 'r') as csvfile:
        myreader = csv.reader(csvfile, delimiter=' ', quotechar='|', skipinitialspace=True)

        row = next(myreader)
        n = int(row[0]) # number of orbits

        for i in range(n):
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
    return n , orbits


# PC Orbits with y' < 0
filename = "../cmake-build-debug/Ass13.out1"
n1 , orbits1 = readOrbits(filename)

# PC Orbits with y' > 0
filename = "../cmake-build-debug/Ass13.out2"
n2 , orbits2 = readOrbits(filename)

yp1 = np.linspace(0.001 , 5.0 , n1)
yp2 = np.linspace(0.001 , 5.0 , n2)
print(yp1)
print(yp2)

nstart = 3
fig1 = plt.figure()
ax1 = fig1.add_subplot(1,1,1)
for i in range(0,n1):
    ax1.plot(orbits1[i].x[nstart:] , orbits1[i].y[nstart:] , '.r' , markersize=0.5)

fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)
for i in range(0,n2):
    ax2.plot(orbits2[i].x[nstart:] , orbits2[i].y[nstart:] , '.r' , markersize=0.5)


# nPlots = 4
# skipPlots = 20
# fig1 = plt.figure()
# axes1 = []
# for i in range(0,nPlots):
#     j = i * skipPlots
#     axes1.append(fig1.add_subplot(nPlots ,1 ,i+1))
#     axes1[i].plot(orbits1[j].x , orbits1[j].y)
#
# fig2 = plt.figure()
# axes2 = []
# for i in range(0,nPlots):
#     j = i * skipPlots
#     axes2.append(fig2.add_subplot(nPlots ,1 ,i+1))
#     axes2[i].plot(orbits2[j].x , orbits2[j].y)




plt.show()