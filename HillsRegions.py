import numpy as np

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

from math import *



def Jacobi(x,mu):
    r1 = math.sqrt((x[0]-mu)**2 + x[1]**2)
    r2 = math.sqrt((x[0]-mu+1)**2 + x[1]**2)
    Omega = 0.5*(x[0]**2 + x[1]**2) + (1-mu)/r1 + mu/r2 + 0.5*mu*(1-mu)

    C = 2.0*Omega - (x[2]**2 + x[3]**2)
    return C

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



mu = 0.1
tol = 1.e-15
ncross = 1

L1,L2,L3 = solveL123(mu,tol)
print(L1,L2,L3)
L1 = np.array([L1 , 0 , 0 , 0])
L2 = np.array([L2 , 0 , 0 , 0])
L3 = np.array([L3 , 0 , 0 , 0])
L4 = np.array([mu-1.0/2.0 , math.sqrt(3)/2.0 , 0 , 0])


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