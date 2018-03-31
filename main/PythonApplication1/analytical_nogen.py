# This script gives the analytical solution of 2D axisymmetric heat conduction of a cylinder (r,z)

import numpy as np
import scipy
from numpy import pi
from scipy import special
#H=1
#R=1
#r=0.5
#z=1.2
#alpha=1
#t=1
#hzH=10
#kz=0.25
def calc_nogen(H,R,z,r,t,alpha,hzH,kz):
    
    # Planar part of the solution (z direction)
    L=H/2
    xstar=z-L
    Bi=hzH*L/kz
    Fo=alpha*t/(L**2)
    zeta1=0.8633
    P_C1=4*np.sin(zeta1)/(2*zeta1+np.sin(2*zeta1))
    P=P_C1*np.exp(-(zeta1**2)*Fo)*np.cos(zeta1*xstar)
    print(P)

    # Cylindrical part of the solution (r direction)
    rstar=r/R
    zeta1=0.86033
    C_C1=(2/zeta1)*scipy.special.j1(zeta1)/(scipy.special.j0(zeta1)**2+scipy.special.j1(zeta1)**2)
    C=C_C1*np.exp(-(zeta1**2)*Fo)*scipy.special.j0(zeta1*rstar)
    print(C)

    # Overall solution is the combination of planar and cylindrical
    Theta=P*C
    return Theta
