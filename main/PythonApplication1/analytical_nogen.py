# This script gives the analytical solution of 2D axisymmetric heat conduction of a cylinder (r,z)

import numpy as np
import scipy
from numpy import pi
from scipy import special
from scipy import optimize
#H=1
#R=1
#r=0.5
#z=1.2
#alpha=1
#t=1
#hzH=10
#hR=10
#kz=0.25
#kr=0.25
def calc_nogen(H,R,z,r,t,rho,cp,hzH,hR,kz,kr):
    
    # Planar part of the solution (z direction)
    L=H/2
    xstar=z-L
    Bi=hzH*L/kz
    alpha=kz/rho/cp
    print(alpha)
    Fo=alpha*t/(L**2)
    zeta1=scipy.optimize.fsolve(eigfcnplanar,1,args=(Bi))
    print("Bi planar is %2.2f" %(Bi))
    print(zeta1)
    P_C1=4*np.sin(zeta1)/(2*zeta1+np.sin(2*zeta1))
    print("P_C1 is %3.3f" %(P_C1))
    P=P_C1*np.exp(-(zeta1**2)*Fo)*np.cos(zeta1*xstar)
    print("P is %4.4f" %(P))

    # Cylindrical part of the solution (r direction)
    rstar=r/R
    Bi=hR*R/kr
    print("Bi circular is %2.0f" %(Bi))
    alpha=kr/rho/cp
    print(alpha)
    Fo=alpha*t/(R**2)
    print("Fo=%2.2f" %(Fo))
    zeta1=scipy.optimize.fsolve(eigfcncircular,1,args=(Bi))
    print(zeta1)
    C_C1=(2/zeta1)*scipy.special.j1(zeta1)/(scipy.special.j0(zeta1)**2+scipy.special.j1(zeta1)**2)
    print("C_C1 is %3.3f" %(C_C1))
    C=C_C1*np.exp(-(zeta1**2)*Fo)*scipy.special.j0(zeta1*rstar)
    print("C is %4.4f" %(C))

    # Overall solution is the combination of planar and cylindrical
    Theta=P*C
    return Theta

def eigfcnplanar(zeta,Bi):
    return (zeta*np.tan(zeta)-Bi)
def eigfcncircular(zeta,Bi):
    return (zeta*scipy.special.j1(zeta)/scipy.special.j0(zeta)-Bi)

#calc_nogen(H,R,z,r,t,rho,cp,hzH,hR,kz,kr)