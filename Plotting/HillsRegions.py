import numpy as np

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

from math import *

from skimage import measure

######################################################################33
# Functions to find the equilibrium points
def F1(x,mu):
    y = ( mu * (1 - x)**2 / (3 - 2*mu - x*(3 - mu - x)))**(1/3)
    return y

def F2(x,mu):
    y = ( mu * (1 + x)**2 / (3 - 2*mu + x*(3 - mu + x)))**(1/3)
    return y

def F3(x,mu):
    y = ( (1 - mu) * (1 + x)**2 / (1 + 2*mu + x*(2 + mu + x)))**(1/3)  
    return y 

def solveL123(mu , tol):
    x02 = (mu / (3 *(1-mu)))**(1/3)
    x01 = (mu / (3 *(1-mu)))**(1/3)
    x03 = 1 - 7/12 * mu

    x1 = x01
    while abs(x1-F1(x1,mu)) > tol:
        x1 = F1(x1,mu)
    x2 = x02
    while abs(x2-F2(x2,mu)) > tol:
        x2 = F2(x2,mu)
    x3 = x03
    while abs(x3-F3(x3,mu)) > tol:
        x3 = F3(x3,mu)

    L1 = mu - 1 - x1
    L2 = mu - 1 + x2
    L3 = mu + x3 
    return L1 , L2 , L3

######################################################################

# Jacobi constant 
def Jacobi(x,mu):
    r1 = sqrt((x[0]-mu)**2 + x[1]**2)
    r2 = sqrt((x[0]-mu+1)**2 + x[1]**2)
    Omega = 0.5*(x[0]**2 + x[1]**2) + (1-mu)/r1 + mu/r2 + 0.5*mu*(1-mu)

    C = 2.0*Omega - (x[2]**2 + x[3]**2)
    return C

# Hills function
def hillsF(x,y,C,mu):
    r1 = np.sqrt( np.power(x-mu , 2) + np.power(y,2))
    r2 = np.sqrt(np.power(x-mu+1,2) + np.power(y,2))
    Omega = 0.5*(x**2 + y**2) + (1-mu)/r1 + mu/r2 + 0.5*mu*(1-mu)
    F = 2.0 * Omega - C
    return F

########################################################################
########################################################################
# MAIN PROGRAM 
mu = 0.1

# Get equilibrium points and their Jacobi constants
tol = 1.e-15
L1,L2,L3 = solveL123(mu,tol)
print(L1,L2,L3)

xLvec = [L1 , L2 , L3 , mu-1.0/2.0 , mu-1.0/2.0]
yLvec = [ 0 , 0 , 0 , sqrt(3)/2.0 , -sqrt(3)/2.0]

xMvec = [mu , mu-1]
yMvec = [0 ,0]


L1 = np.array([L1 , 0 , 0 , 0])
L2 = np.array([L2 , 0 , 0 , 0])
L3 = np.array([L3 , 0 , 0 , 0])
L4 = np.array([mu-1.0/2.0 , sqrt(3)/2.0 , 0 , 0])



C1 = Jacobi(L2,mu)
C2 = Jacobi(L1,mu)
C3 = Jacobi(L3,mu)
C4 = Jacobi(L4,mu)


# Base grid we will use to draw contour plots
n = 6000
xylim = 3
elemsize = 2*xylim/(n-1)
x = np.linspace(-xylim , xylim , n)
y = np.linspace(-xylim , xylim , n)
X , Y = np.meshgrid(x, y, sparse=False, indexing='ij')
F = hillsF(X,Y,0.0,mu)

# Computing hills function for different cases 
Cvec = [C1+1 , C1 , (C1+C2)/2 , C2 , (C2+C3)/2 , C3 , (C4+C3)/2]
colors = [ 'red' , 'darkorange' , 'green' , 'navy' , 'royalblue' , 'purple' , 'saddlebrown']


fig1 = plt.figure()
ax = fig1.add_subplot(1,1,1)
ax.plot(xLvec , yLvec , '+k')
ax.plot(xMvec , yMvec , '*k')

fig2 = plt.figure()
axes = []

leg = ['C > C1' , 'C = C1' , 'C2 < C < C1' , 'C = C2' , 'C3 < C < C2' , 'C = C3' , 'C4 < C < C3']
for i in range(0 , 7):
    axes.append(fig2.add_subplot(4,2,i+1))
    C = Cvec[i]
    contours = measure.find_contours(F, C)

    axes[i].plot(xLvec , yLvec , '+k')
    axes[i].plot(xMvec , yMvec , '*k')
    axes[i].set_title(leg[i])

    for n2, contour in enumerate(contours):
        contour = elemsize*(contour - n/2)
        ax.plot(contour[:, 0], contour[:, 1], linewidth=2 , color=colors[i])
        axes[i].plot(contour[:, 0], contour[:, 1], linewidth=2 , color=colors[i])



plt.show()