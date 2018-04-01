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
    xstar=(z-L)/L
    Bi=hzH*L/kz
    print("Bi planar is %2.2f" %(Bi))
    alpha=kz/rho/cp
    #print("Alpha is %4.4f" %(alpha))
    Fo=alpha*t/(L**2)
    print("Fo is %4.4f" %(Fo))
    if Fo<0.2:
        zeta1=scipy.optimize.bisect(eigfcnplanar,0,1.57,args=(Bi,))
        zeta2=scipy.optimize.bisect(eigfcnplanar,3.1415,4.71,args=(Bi,))
        zeta3=scipy.optimize.bisect(eigfcnplanar,6.2832,7.85,args=(Bi,))
        zeta4=scipy.optimize.bisect(eigfcnplanar,9.4248,10.99,args=(Bi,))
        #zeta1=scipy.optimize.fsolve(eigfcnplanar,1,args=(Bi,))
        #zeta2=scipy.optimize.fsolve(eigfcnplanar,4,args=(Bi,))
        #zeta3=scipy.optimize.fsolve(eigfcnplanar,7,args=(Bi,))
        #zeta4=scipy.optimize.fsolve(eigfcnplanar,10.2,args=(Bi,))
        print("zeta1 is %4.4f" %(zeta1))
        print("zeta2 is %4.4f" %(zeta2))
        print("zeta3 is %4.4f" %(zeta3))
        print("zeta4 is %4.4f" %(zeta4))
        P_C1=P_Cn(zeta1)
        P_C2=P_Cn(zeta2)
        P_C3=P_Cn(zeta3)
        P_C4=P_Cn(zeta4)
        #P_C1=4*np.sin(zeta1)/(2*zeta1+np.sin(2*zeta1))
        #print("P_C1 is %3.3f" %(P_C1))
        #print("P_C2 is %3.3f" %(P_C2))
        #print("P_C3 is %3.3f" %(P_C3))
        #print("P_C4 is %3.3f" %(P_C4))
        P=P_C1*np.exp(-(zeta1**2)*Fo)*np.cos(zeta1*xstar)+P_C2*np.exp(-(zeta2**2)*Fo)*np.cos(zeta2*xstar)+P_C3*np.exp(-(zeta3**2)*Fo)*np.cos(zeta3*xstar)+P_C4*np.exp(-(zeta4**2)*Fo)*np.cos(zeta4*xstar)
        print("P is %4.4f" %(P))
        if P>1:
            print("Unrealistic P is found!")
    else:
        zeta1=scipy.optimize.fsolve(eigfcnplanar,1,args=(Bi))
        print("zeta1 is %4.4f" %(zeta1))
        P_C1=P_Cn(zeta1)
        #P_C1=4*np.sin(zeta1)/(2*zeta1+np.sin(2*zeta1))
        print("P_C1 is %3.3f" %(P_C1))
        P=P_C1*np.exp(-(zeta1**2)*Fo)*np.cos(zeta1*xstar)
        print("P is %4.4f" %(P))
        if P>1:
            print("Unrealistic P is found!")
   
    # Cylindrical part of the solution (r direction)
    rstar=r/R
    Bi=hR*R/kr
    print("Bi circular is %2.2f" %(Bi))
    alpha=kr/rho/cp
    #print("Alpha is %4.4f" %(alpha))
    Fo=alpha*t/(R**2)
    print("Fo=%4.4f" %(Fo))
    if Fo<0.2:
        zeta1=scipy.optimize.bisect(eigfcncircular,0,2.40,args=(Bi,))
        zeta2=scipy.optimize.bisect(eigfcncircular,3.8317,5.52,args=(Bi,))
        zeta3=scipy.optimize.bisect(eigfcncircular,7.0156,8.65,args=(Bi,))
        zeta4=scipy.optimize.bisect(eigfcncircular,10.1735,11.79,args=(Bi,))
        zeta5=scipy.optimize.bisect(eigfcncircular,13.3237,14.93,args=(Bi,))
        zeta6=scipy.optimize.bisect(eigfcncircular,16.470,18.07,args=(Bi,))
        #print("zeta1 is %4.4f" %(zeta1))
        #print("zeta2 is %4.4f" %(zeta2))
        #print("zeta3 is %4.4f" %(zeta3))
        #print("zeta4 is %4.4f" %(zeta4))
        #print("zeta5 is %4.4f" %(zeta5))
        #print("zeta6 is %4.4f" %(zeta6))
        #zeta1=scipy.optimize.fsolve(eigfcncircular,2,args=(Bi),xtol=1.0e-12)
        #zeta1=scipy.optimize.fmin_tnc(eigfcncircular,(1,),args=(Bi,),bounds=((0.1,2.41),))
        #zeta1=scipy.optimize.fminbound(eigfcncircular,0,2.41,args=(Bi,),xtol=1e-08)
        C_C1=C_Cn(zeta1)
        C_C2=C_Cn(zeta2)
        C_C3=C_Cn(zeta3)
        C_C4=C_Cn(zeta4)
        C_C5=C_Cn(zeta5)
        C_C6=C_Cn(zeta6)
        #C_C1=(2/zeta1)*scipy.special.j1(zeta1)/(scipy.special.j0(zeta1)**2+scipy.special.j1(zeta1)**2)
        #print("C_C1 is %3.3f" %(C_C1))
        #print("C_C2 is %3.3f" %(C_C2))
        #print("C_C3 is %3.3f" %(C_C3))
        #print("C_C4 is %3.3f" %(C_C4))
        #print("C_C5 is %3.3f" %(C_C5))
        #print("C_C6 is %3.3f" %(C_C6))
        C=C_C1*np.exp(-(zeta1**2)*Fo)*scipy.special.j0(zeta1*rstar)+C_C2*np.exp(-(zeta2**2)*Fo)*scipy.special.j0(zeta2*rstar)+C_C3*np.exp(-(zeta3**2)*Fo)*scipy.special.j0(zeta3*rstar)+C_C4*np.exp(-(zeta4**2)*Fo)*scipy.special.j0(zeta4*rstar)+C_C5*np.exp(-(zeta5**2)*Fo)*scipy.special.j0(zeta5*rstar)+C_C6*np.exp(-(zeta6**2)*Fo)*scipy.special.j0(zeta6*rstar)
        print("C is %4.4f" %(C))
        if C>1:
            print("Unrealistic C is found!")
    else:
        zeta1=1.9898
        print("zeta1 is %4.4f" %(zeta1))
        #zeta1=scipy.optimize.fsolve(eigfcncircular,1,args=(Bi),xtol=1.0e-12)
        zeta1=scipy.optimize.fmin_tnc(eigfcncircular,1,args=(Bi),bounds=((0.0,2.41),))
        if zeta1>2.41:
            print("Unrealistic zeta1 is found!")
        C_C1=C_Cn(zeta1)
        #C_C1=(2/zeta1)*scipy.special.j1(zeta1)/(scipy.special.j0(zeta1)**2+scipy.special.j1(zeta1)**2)
        print("C_C1 is %3.3f" %(C_C1))
        C=C_C1*np.exp(-(zeta1**2)*Fo)*scipy.special.j0(zeta1*rstar)
        print("C is %4.4f" %(C))
        if C>1:
            print("Unrealistic C is found!")
    # Overall solution is the combination of planar and cylindrical
    Theta=P*C
    return Theta

def eigfcnplanar(zeta,Bi):
    return ((zeta*np.tan(zeta))-Bi)
def eigfcncircular(zeta,Bi):
    return (zeta*scipy.special.j1(zeta)/scipy.special.j0(zeta)-Bi)
def P_Cn(zeta):
    return (4*np.sin(zeta)/(2*zeta+np.sin(2*zeta)))
def C_Cn(zeta):
    return ((2/zeta)*scipy.special.j1(zeta)/(scipy.special.j0(zeta)**2+scipy.special.j1(zeta)**2))
#calc_nogen(H,R,z,r,t,rho,cp,hzH,hR,kz,kr)