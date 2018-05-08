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
def calc_probe(H,R,z,r,t,rho,cp,hzH,hR,kz,kr,Tinit,Tamb):
    Tavg_all=np.empty((tsample.shape[0],))*np.nan
    Err_all=np.empty((tsample.shape[0],3))*np.nan
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
        zeta5=scipy.optimize.bisect(eigfcnplanar,12.5664,14.137,args=(Bi,))
        zeta6=scipy.optimize.bisect(eigfcnplanar,15.7080,17.278,args=(Bi,))
        zeta7=scipy.optimize.bisect(eigfcnplanar,18.8497,20.41,args=(Bi,))
        #zeta1=scipy.optimize.fsolve(eigfcnplanar,1,args=(Bi,))
        #zeta2=scipy.optimize.fsolve(eigfcnplanar,4,args=(Bi,))
        #zeta3=scipy.optimize.fsolve(eigfcnplanar,7,args=(Bi,))
        #zeta4=scipy.optimize.fsolve(eigfcnplanar,10.2,args=(Bi,))
        #print("zeta1 is %4.4f" %(zeta1))
        #print("zeta2 is %4.4f" %(zeta2))
        #print("zeta3 is %4.4f" %(zeta3))
        #print("zeta4 is %4.4f" %(zeta4))
        P_C1=P_Cn(zeta1)
        P_C2=P_Cn(zeta2)
        P_C3=P_Cn(zeta3)
        P_C4=P_Cn(zeta4)
        P_C5=P_Cn(zeta5)
        P_C6=P_Cn(zeta6)
        P_C7=P_Cn(zeta7)
        #P_C1=4*np.sin(zeta1)/(2*zeta1+np.sin(2*zeta1))
        #print("P_C1 is %3.3f" %(P_C1))
        #print("P_C2 is %3.3f" %(P_C2))
        #print("P_C3 is %3.3f" %(P_C3))
        #print("P_C4 is %3.3f" %(P_C4))
        P6=P_C1*np.exp(-(zeta1**2)*Fo)*np.cos(zeta1*xstar)+P_C2*np.exp(-(zeta2**2)*Fo)*np.cos(zeta2*xstar)+P_C3*np.exp(-(zeta3**2)*Fo)*np.cos(zeta3*xstar)+P_C4*np.exp(-(zeta4**2)*Fo)*np.cos(zeta4*xstar)+P_C5*np.exp(-(zeta5**2)*Fo)*np.cos(zeta5*xstar)+P_C6*np.exp(-(zeta6**2)*Fo)*np.cos(zeta6*xstar)   
        P7=P_C1*np.exp(-(zeta1**2)*Fo)*np.cos(zeta1*xstar)+P_C2*np.exp(-(zeta2**2)*Fo)*np.cos(zeta2*xstar)+P_C3*np.exp(-(zeta3**2)*Fo)*np.cos(zeta3*xstar)+P_C4*np.exp(-(zeta4**2)*Fo)*np.cos(zeta4*xstar)+P_C5*np.exp(-(zeta5**2)*Fo)*np.cos(zeta5*xstar)+P_C6*np.exp(-(zeta6**2)*Fo)*np.cos(zeta6*xstar)+P_C7*np.exp(-(zeta7**2)*Fo)*np.cos(zeta7*xstar)   
        Err6_7=np.abs(Tinit-Tamb)*np.abs(P6-P7)/P7
        print(Err6_7)
        if P>1:
            print("P is %4.4f" %(P))
            print("Unrealistic P is found!")
    else:
        zeta1=scipy.optimize.fsolve(eigfcnplanar,1,args=(Bi))
        #print("zeta1 is %4.4f" %(zeta1))
        P_C1=P_Cn(zeta1)
        #print("P_C1 is %3.3f" %(P_C1))
        P=P_C1*np.exp(-(zeta1**2)*Fo)*np.cos(zeta1*xstar)    
        if P>1:
            print("P is %4.4f" %(P))
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
        #print("C_C1 is %3.3f" %(C_C1))
        #print("C_C2 is %3.3f" %(C_C2))
        #print("C_C3 is %3.3f" %(C_C3))
        #print("C_C4 is %3.3f" %(C_C4))
        #print("C_C5 is %3.3f" %(C_C5))
        #print("C_C6 is %3.3f" %(C_C6))
        C=C_C1*np.exp(-(zeta1**2)*Fo)*scipy.special.j0(zeta1*rstar)+C_C2*np.exp(-(zeta2**2)*Fo)*scipy.special.j0(zeta2*rstar)+C_C3*np.exp(-(zeta3**2)*Fo)*scipy.special.j0(zeta3*rstar)+C_C4*np.exp(-(zeta4**2)*Fo)*scipy.special.j0(zeta4*rstar)+C_C5*np.exp(-(zeta5**2)*Fo)*scipy.special.j0(zeta5*rstar)+C_C6*np.exp(-(zeta6**2)*Fo)*scipy.special.j0(zeta6*rstar)     
        if C>1:
            print("C is %4.4f" %(C))
            print("Unrealistic C is found!")
    else:
        zeta1=scipy.optimize.bisect(eigfcncircular,0,2.40,args=(Bi,))
        #print("zeta1 is %4.4f" %(zeta1))
        #zeta1=scipy.optimize.fsolve(eigfcncircular,1,args=(Bi),xtol=1.0e-12)
        #zeta1=scipy.optimize.fmin_tnc(eigfcncircular,1,args=(Bi),bounds=((0.0,2.41),))
        if zeta1>2.41:
            print("Unrealistic zeta1 is found!")
        C_C1=C_Cn(zeta1)
        #print("C_C1 is %3.3f" %(C_C1))
        C=C_C1*np.exp(-(zeta1**2)*Fo)*scipy.special.j0(zeta1*rstar)
        
        if C>1:
            print("C is %4.4f" %(C))
            print("Unrealistic C is found!")
    # Overall solution is the combination of planar and cylindrical
    Theta=P*C
    Temp_analytical=Theta*(Tinit-Tamb)+Tamb
    return Temp_analytical

def calc_whole(H,R,zLoc,rLoc,tsample,rho,cp,hzH,hR,kz,kr,Tinit,Tamb,GridMap,Ncellr,Ncellz):
    Tavg_all=np.empty((tsample.shape[0],))*np.nan
    Err_all=np.empty((tsample.shape[0],3))*np.nan
    Tavg_all[0]=Tinit
    T_analy_all=np.empty((Ncellr*Ncellz,tsample.shape[0]))
    # Calculate zeta_n and C_n for planar
    L=H/2
    Bi_P=hzH*L/kz
    print("Bi planar is %2.2f" %(Bi_P))
    alpha_P=kz/rho/cp
    #print("Alpha is %4.4f" %(alpha_P))
    
    P_zeta1=scipy.optimize.bisect(eigfcnplanar,0,1.57,args=(Bi_P,))
    P_zeta2=scipy.optimize.bisect(eigfcnplanar,3.1415,4.71,args=(Bi_P,))
    P_zeta3=scipy.optimize.bisect(eigfcnplanar,6.2832,7.85,args=(Bi_P,))
    P_zeta4=scipy.optimize.bisect(eigfcnplanar,9.4248,10.99,args=(Bi_P,))
    P_zeta5=scipy.optimize.bisect(eigfcnplanar,12.5664,14.137,args=(Bi_P,))
    P_zeta6=scipy.optimize.bisect(eigfcnplanar,15.7080,17.278,args=(Bi_P,))
    P_zeta7=scipy.optimize.bisect(eigfcnplanar,18.8497,20.41,args=(Bi_P,))
    P_zeta8=scipy.optimize.bisect(eigfcnplanar,21.9912,23.56,args=(Bi_P,))
    P_zeta9=scipy.optimize.bisect(eigfcnplanar,25.1328,26.7,args=(Bi_P,))
    #print("zeta1 is %4.4f" %(P_zeta1))
    #print("zeta2 is %4.4f" %(P_zeta2))
    #print("zeta3 is %4.4f" %(P_zeta3))
    #print("zeta4 is %4.4f" %(P_zeta4))
    #print("zeta5 is %4.4f" %(P_zeta5))
    #print("zeta6 is %4.4f" %(P_zeta6))
    P_C1=P_Cn(P_zeta1)
    P_C2=P_Cn(P_zeta2)
    P_C3=P_Cn(P_zeta3)
    P_C4=P_Cn(P_zeta4)
    P_C5=P_Cn(P_zeta5)
    P_C6=P_Cn(P_zeta6)
    P_C7=P_Cn(P_zeta7)
    P_C8=P_Cn(P_zeta8)
    P_C9=P_Cn(P_zeta9)
        #print("P_C1 is %3.3f" %(P_C1))
        #print("P_C2 is %3.3f" %(P_C2))
        #print("P_C3 is %3.3f" %(P_C3))
        #print("P_C4 is %3.3f" %(P_C4))
    
    # Calculate zeta_n and C_n for cylindrical
    Bi_C=hR*R/kr
    print("Bi circular is %2.2f" %(Bi_C))
    alpha_C=kr/rho/cp
    #print("Alpha is %4.4f" %(alpha))

    C_zeta1=scipy.optimize.bisect(eigfcncircular,0,2.40,args=(Bi_C,))
    C_zeta2=scipy.optimize.bisect(eigfcncircular,3.8317,5.52,args=(Bi_C,))
    C_zeta3=scipy.optimize.bisect(eigfcncircular,7.0156,8.65,args=(Bi_C,))
    C_zeta4=scipy.optimize.bisect(eigfcncircular,10.1735,11.79,args=(Bi_C,))
    C_zeta5=scipy.optimize.bisect(eigfcncircular,13.3237,14.93,args=(Bi_C,))
    C_zeta6=scipy.optimize.bisect(eigfcncircular,16.470,18.07,args=(Bi_C,))
    C_zeta7=scipy.optimize.bisect(eigfcncircular,19.6159,21.21,args=(Bi_C,))
    C_zeta8=scipy.optimize.bisect(eigfcncircular,22.7601,24.35,args=(Bi_C,))
    C_zeta9=scipy.optimize.bisect(eigfcncircular,25.9037,27.49,args=(Bi_C,))
    #print("zeta1 is %4.4f" %(C_zeta1))
    #print("zeta2 is %4.4f" %(C_zeta2))
    #print("zeta3 is %4.4f" %(C_zeta3))
    #print("zeta4 is %4.4f" %(C_zeta4))
    #print("zeta5 is %4.4f" %(C_zeta5))
    #print("zeta6 is %4.4f" %(C_zeta6))
    C_C1=C_Cn(C_zeta1)
    C_C2=C_Cn(C_zeta2)
    C_C3=C_Cn(C_zeta3)
    C_C4=C_Cn(C_zeta4)
    C_C5=C_Cn(C_zeta5)
    C_C6=C_Cn(C_zeta6)
    C_C7=C_Cn(C_zeta7)
    C_C8=C_Cn(C_zeta8)
    C_C9=C_Cn(C_zeta9)

        #print("C_C1 is %3.3f" %(C_C1))
        #print("C_C2 is %3.3f" %(C_C2))
        #print("C_C3 is %3.3f" %(C_C3))
        #print("C_C4 is %3.3f" %(C_C4))
        #print("C_C5 is %3.3f" %(C_C5))
        #print("C_C6 is %3.3f" %(C_C6))
    
        # Here time loop starts 
    T=Tinit*np.ones([Ncellr*Ncellz,1])
    T_analy_all[:,0]=T.reshape(Ncellr*Ncellz,)
    for i in range(1,tsample.shape[0]):
        t=tsample[i]
        print("t_analy=%4.4f" %(t))
        Fo_P=alpha_P*t/(L**2)
        #print("Fo planar is %4.4f" %(Fo_P))
        Fo_C=alpha_C*t/(R**2)
        #print("Fo cylindrical is %4.4f" %(Fo_C))
        count=0
        for m in range(0,Ncellz):
            for n in range(0,Ncellr):
    # PLanar part of the solution (z direction)
                xstar=np.abs((zLoc[m,n]-L)/L)
                if xstar>0:
                    if Fo_P<0.2:
                        P6=P_C1*np.exp(-(P_zeta1**2)*Fo_P)*np.cos(P_zeta1*xstar)+P_C2*np.exp(-(P_zeta2**2)*Fo_P)*np.cos(P_zeta2*xstar)+P_C3*np.exp(-(P_zeta3**2)*Fo_P)*np.cos(P_zeta3*xstar)+P_C4*np.exp(-(P_zeta4**2)*Fo_P)*np.cos(P_zeta4*xstar)+P_C5*np.exp(-(P_zeta5**2)*Fo_P)*np.cos(P_zeta5*xstar)+P_C6*np.exp(-(P_zeta6**2)*Fo_P)*np.cos(P_zeta6*xstar)
                        P7=P_C1*np.exp(-(P_zeta1**2)*Fo_P)*np.cos(P_zeta1*xstar)+P_C2*np.exp(-(P_zeta2**2)*Fo_P)*np.cos(P_zeta2*xstar)+P_C3*np.exp(-(P_zeta3**2)*Fo_P)*np.cos(P_zeta3*xstar)+P_C4*np.exp(-(P_zeta4**2)*Fo_P)*np.cos(P_zeta4*xstar)+P_C5*np.exp(-(P_zeta5**2)*Fo_P)*np.cos(P_zeta5*xstar)+P_C6*np.exp(-(P_zeta6**2)*Fo_P)*np.cos(P_zeta6*xstar)+P_C7*np.exp(-(P_zeta7**2)*Fo_P)*np.cos(P_zeta7*xstar)
                        P8=P_C1*np.exp(-(P_zeta1**2)*Fo_P)*np.cos(P_zeta1*xstar)+P_C2*np.exp(-(P_zeta2**2)*Fo_P)*np.cos(P_zeta2*xstar)+P_C3*np.exp(-(P_zeta3**2)*Fo_P)*np.cos(P_zeta3*xstar)+P_C4*np.exp(-(P_zeta4**2)*Fo_P)*np.cos(P_zeta4*xstar)+P_C5*np.exp(-(P_zeta5**2)*Fo_P)*np.cos(P_zeta5*xstar)+P_C6*np.exp(-(P_zeta6**2)*Fo_P)*np.cos(P_zeta6*xstar)+P_C7*np.exp(-(P_zeta7**2)*Fo_P)*np.cos(P_zeta7*xstar)+P_C8*np.exp(-(P_zeta8**2)*Fo_P)*np.cos(P_zeta8*xstar)
                        P9=P_C1*np.exp(-(P_zeta1**2)*Fo_P)*np.cos(P_zeta1*xstar)+P_C2*np.exp(-(P_zeta2**2)*Fo_P)*np.cos(P_zeta2*xstar)+P_C3*np.exp(-(P_zeta3**2)*Fo_P)*np.cos(P_zeta3*xstar)+P_C4*np.exp(-(P_zeta4**2)*Fo_P)*np.cos(P_zeta4*xstar)+P_C5*np.exp(-(P_zeta5**2)*Fo_P)*np.cos(P_zeta5*xstar)+P_C6*np.exp(-(P_zeta6**2)*Fo_P)*np.cos(P_zeta6*xstar)+P_C7*np.exp(-(P_zeta7**2)*Fo_P)*np.cos(P_zeta7*xstar)+P_C8*np.exp(-(P_zeta8**2)*Fo_P)*np.cos(P_zeta8*xstar)+P_C9*np.exp(-(P_zeta9**2)*Fo_P)*np.cos(P_zeta9*xstar)
                        P=P9
                    else:
                        P6=P_C1*np.exp(-(P_zeta1**2)*Fo_P)*np.cos(P_zeta1*xstar)+P_C2*np.exp(-(P_zeta2**2)*Fo_P)*np.cos(P_zeta2*xstar)+P_C3*np.exp(-(P_zeta3**2)*Fo_P)*np.cos(P_zeta3*xstar)+P_C4*np.exp(-(P_zeta4**2)*Fo_P)*np.cos(P_zeta4*xstar)+P_C5*np.exp(-(P_zeta5**2)*Fo_P)*np.cos(P_zeta5*xstar)+P_C6*np.exp(-(P_zeta6**2)*Fo_P)*np.cos(P_zeta6*xstar)
                        P7=P6
                        P8=P7
                        P9=P8
                        P=P6
                    #Err6_7=np.abs(Tinit-Tamb)*np.abs(P6-P7)/P7
                    #print(Err6_7)
                    #if P>1:
                    #    print("P is %4.4f" %(P))
                    #    print("Unrealistic P is found!")
                    #else:
                    #    print("P is %4.4f" %(P))
                    #print("P is %4.4f" %(P))

    # Cylindrical part of the solution (r direction)
                rstar=rLoc[m,n]/R

                if Fo_C<0.2:
                    C6=C_C1*np.exp(-(C_zeta1**2)*Fo_C)*scipy.special.j0(C_zeta1*rstar)+C_C2*np.exp(-(C_zeta2**2)*Fo_C)*scipy.special.j0(C_zeta2*rstar)+C_C3*np.exp(-(C_zeta3**2)*Fo_C)*scipy.special.j0(C_zeta3*rstar)+C_C4*np.exp(-(C_zeta4**2)*Fo_C)*scipy.special.j0(C_zeta4*rstar)+C_C5*np.exp(-(C_zeta5**2)*Fo_C)*scipy.special.j0(C_zeta5*rstar)+C_C6*np.exp(-(C_zeta6**2)*Fo_C)*scipy.special.j0(C_zeta6*rstar)
                    C7=C_C1*np.exp(-(C_zeta1**2)*Fo_C)*scipy.special.j0(C_zeta1*rstar)+C_C2*np.exp(-(C_zeta2**2)*Fo_C)*scipy.special.j0(C_zeta2*rstar)+C_C3*np.exp(-(C_zeta3**2)*Fo_C)*scipy.special.j0(C_zeta3*rstar)+C_C4*np.exp(-(C_zeta4**2)*Fo_C)*scipy.special.j0(C_zeta4*rstar)+C_C5*np.exp(-(C_zeta5**2)*Fo_C)*scipy.special.j0(C_zeta5*rstar)+C_C6*np.exp(-(C_zeta6**2)*Fo_C)*scipy.special.j0(C_zeta6*rstar)+C_C7*np.exp(-(C_zeta7**2)*Fo_C)*scipy.special.j0(C_zeta7*rstar)
                    C8=C_C1*np.exp(-(C_zeta1**2)*Fo_C)*scipy.special.j0(C_zeta1*rstar)+C_C2*np.exp(-(C_zeta2**2)*Fo_C)*scipy.special.j0(C_zeta2*rstar)+C_C3*np.exp(-(C_zeta3**2)*Fo_C)*scipy.special.j0(C_zeta3*rstar)+C_C4*np.exp(-(C_zeta4**2)*Fo_C)*scipy.special.j0(C_zeta4*rstar)+C_C5*np.exp(-(C_zeta5**2)*Fo_C)*scipy.special.j0(C_zeta5*rstar)+C_C6*np.exp(-(C_zeta6**2)*Fo_C)*scipy.special.j0(C_zeta6*rstar)+C_C7*np.exp(-(C_zeta7**2)*Fo_C)*scipy.special.j0(C_zeta7*rstar)+C_C8*np.exp(-(C_zeta8**2)*Fo_C)*scipy.special.j0(C_zeta8*rstar)
                    C9=C_C1*np.exp(-(C_zeta1**2)*Fo_C)*scipy.special.j0(C_zeta1*rstar)+C_C2*np.exp(-(C_zeta2**2)*Fo_C)*scipy.special.j0(C_zeta2*rstar)+C_C3*np.exp(-(C_zeta3**2)*Fo_C)*scipy.special.j0(C_zeta3*rstar)+C_C4*np.exp(-(C_zeta4**2)*Fo_C)*scipy.special.j0(C_zeta4*rstar)+C_C5*np.exp(-(C_zeta5**2)*Fo_C)*scipy.special.j0(C_zeta5*rstar)+C_C6*np.exp(-(C_zeta6**2)*Fo_C)*scipy.special.j0(C_zeta6*rstar)+C_C7*np.exp(-(C_zeta7**2)*Fo_C)*scipy.special.j0(C_zeta7*rstar)+C_C8*np.exp(-(C_zeta8**2)*Fo_C)*scipy.special.j0(C_zeta8*rstar)+C_C9*np.exp(-(C_zeta9**2)*Fo_C)*scipy.special.j0(C_zeta9*rstar)
                    C=C9
                else:
                    C6=C_C1*np.exp(-(C_zeta1**2)*Fo_C)*scipy.special.j0(C_zeta1*rstar)+C_C2*np.exp(-(C_zeta2**2)*Fo_C)*scipy.special.j0(C_zeta2*rstar)+C_C3*np.exp(-(C_zeta3**2)*Fo_C)*scipy.special.j0(C_zeta3*rstar)+C_C4*np.exp(-(C_zeta4**2)*Fo_C)*scipy.special.j0(C_zeta4*rstar)+C_C5*np.exp(-(C_zeta5**2)*Fo_C)*scipy.special.j0(C_zeta5*rstar)+C_C6*np.exp(-(C_zeta6**2)*Fo_C)*scipy.special.j0(C_zeta6*rstar)
                    C7=C6
                    C8=C7
                    C9=C8
                    C=C6
                #if C>1:
                #    print("C is %4.4f" %(C))
                #    print("Unrealistic C is found!")
                Temp6=P6*C6*(Tinit-Tamb)+Tamb
                Temp7=P7*C7*(Tinit-Tamb)+Tamb
                Temp8=P8*C8*(Tinit-Tamb)+Tamb
                Temp9=P9*C9*(Tinit-Tamb)+Tamb
                Err6_9=np.abs(Temp6-Temp9)*100/Temp9
                Err7_9=np.abs(Temp7-Temp9)*100/Temp9
                Err8_9=np.abs(Temp8-Temp9)*100/Temp9
                #print(Err6_9)
                #print(Err7_9)
                #print(Err8_9)
                #print("C is %4.4f" %(C))
    # Overall solution is the combination of planar and cylindrical
                #Theta=P*C
                #print("For xstar=%4.4f and rstar=%4.4f, P is %4.4f and C is %4.4f " %(xstar,rstar,P,C))
                #Temp_analytical=Theta*(Tinit-Tamb)+Tamb
                Temp_analytical=Temp9
                T_analy_all[count,i]=Temp_analytical
                if i==1:
                    Err_all[i,:]=[Err6_9,Err7_9,Err8_9]
                if Err_all[i,0]<Err6_9:
                    Err_all[i,0]=Err6_9
                if Err_all[i,1]<Err7_9:
                    Err_all[i,1]=Err7_9
                if Err_all[i,2]<Err8_9:
                    Err_all[i,2]=Err8_9
                count=count+1
        
        Tavg=np.mean(T_analy_all[:,i],axis=0)
        print(Tavg)
        Tavg_all[i]=Tavg

    return (T_analy_all,Tavg_all,Err_all)


def eigfcnplanar(zeta,Bi):
    return ((zeta*np.tan(zeta))-Bi)
def eigfcncircular(zeta,Bi):
    return (zeta*scipy.special.j1(zeta)/scipy.special.j0(zeta)-Bi)
def P_Cn(zeta):
    return (4*np.sin(zeta)/(2*zeta+np.sin(2*zeta)))
def C_Cn(zeta):
    return ((2/zeta)*scipy.special.j1(zeta)/(scipy.special.j0(zeta)**2+scipy.special.j1(zeta)**2))
