import numpy as np

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

from math import *


def hillsF(x,y,C,mu):
    r1 = sqrt( (x-mu)**2 + y**2)
    r2 = sqrt((x-mu+1)**2 + y**2)
    Omega = 0.5*(x**2 + y**2) + (1-mu)/r1 + mu/r2 + 0.5*mu*(1-mu)
    F = 2.0 * Omega - C
    return F

def hillsFx(x,y,C,mu):
    r1 = sqrt( (x[0]-mu)**2 + x[1]**2)
    r2 = sqrt((x[0]-mu+1)**2 + x[1]**2)
    Ox = x[0] - ((1.0-mu)*(x[0]-mu))/(r1**3) - (mu*(x[0]-mu+1))/(r2**3)
    Fx = 2.0 * Ox
    return Fx

def hillsFy(x,y,C,mu):
    r1 = sqrt( (x[0]-mu)**2 + x[1]**2)
    r2 = sqrt((x[0]-mu+1)**2 + x[1]**2)
    Oy = x[1]*( 1 - (1-mu)/(r1**3) - mu/(r2**3) )
    Fy = 2.0 * Oy
    return Fy


mu = 0.3
C = 1000


n = 1000
x = np.linspace(-1 , 1 , n)
y = np.linspace(-1 , 1 , n)
X , Y = np.meshgrid(x, y, sparse=False, indexing='ij')

F = np.zeros((n,n))
for i in range(0,n):
	for j in range(0,n):
		F[i,j] = hillsF(X[i,j] , Y[i,j] , C , mu)

Z = F >= 0
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,linewidth=0, antialiased=False)


plt.show()