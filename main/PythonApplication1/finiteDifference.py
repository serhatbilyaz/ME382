#This script contains all subroutines that solve the problem with finite difference approach

import numpy as np
#import scipy
from numpy import pi
#from scipy import special

def uniformgen(GridMap,kr,kz,hR,hz0,hzH,rho,cp,R,H,Ncellr,Ncellz,delr,delz,delt,Tinit,Tamb,tsample,qgen):
    Coeff=np.zeros((Ncellr*Ncellz,Ncellr*Ncellz))
    b=np.zeros((Ncellr*Ncellz,1))
    T_all=np.empty((Ncellr*Ncellz,tsample.shape[0]))
    T=Tinit*np.ones([Ncellr*Ncellz,1])
    T_all[:,0]=T.reshape(Ncellr*Ncellz,)
    for i in range(1,tsample.shape[0]):
        t=tsample[i]
        print("t_FD=%4.4f" %(t))  
# Coefficient and forcing matrix for interior nodes
        for m in range(2,Ncellr):
            qwestcoeff=kr*2*np.pi*(m-0.5)*delz
            qeastcoeff=kr*2*np.pi*(m+0.5)*delz
            qnorthsouthcoeff=kz*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)/delz
            Estcoeff=rho*cp*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz/delt
            Qgencoeff=qgen*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz
            for n in range(2,Ncellz):
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qeastcoeff-qwestcoeff-2*qnorthsouthcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m-2,n-1]]=qwestcoeff  # Coeff for m-1,n
                Coeff[GridMap[m-1,n-1],GridMap[m,n-1]]=qeastcoeff  # Coeff for m+1,n
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-2]]=qnorthsouthcoeff  # Coeff for m,n-1
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n]]=qnorthsouthcoeff  # Coeff for m,n+1
                b[GridMap[m-1,n-1],0]=-Estcoeff*T[GridMap[m-1,n-1]]-Qgencoeff

# Coefficient and forcing matrix for boundary nodes
# r=0 boundary
        m=1
        qeastcoeff=kr*2*np.pi*(m+0.5)*delz
        qnorthsouthcoeff=kz*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)/delz
        Estcoeff=rho*cp*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz/delt
        Qgencoeff=qgen*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz
        for n in range(1,Ncellz+1):
            Coeff[GridMap[m-1,n-1],GridMap[m,n-1]]=qeastcoeff  # Coeff for m+1,n
            b[GridMap[m-1,n-1],0]=-Estcoeff*T[GridMap[m-1,n-1]]-Qgencoeff
            if n==1:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qeastcoeff-qnorthsouthcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n]]=qnorthsouthcoeff  # Coeff for m,n+1
            elif n==Ncellz:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qeastcoeff-qnorthsouthcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-2]]=qnorthsouthcoeff  # Coeff for m,n-1
            else:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qeastcoeff-2*qnorthsouthcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-2]]=qnorthsouthcoeff  # Coeff for m,n-1
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n]]=qnorthsouthcoeff  # Coeff for m,n+1

        # r=R boundary
        m=Ncellr
        qwestcoeff=kr*2*np.pi*(m-0.5)*delz
        qnorthsouthcoeff=kz*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)/delz
        Estcoeff=rho*cp*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz/delt
        Convcoeff=hR*2*np.pi*R*delz
        Qgencoeff=qgen*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz
        for n in range(1,Ncellz+1):
            Coeff[GridMap[m-1,n-1],GridMap[m-2,n-1]]=qwestcoeff+0.5*Convcoeff  # Coeff for m-1,n
            b[GridMap[m-1,n-1],0]=-Estcoeff*T[GridMap[m-1,n-1]]-Convcoeff*Tamb-Qgencoeff
            if n==1:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qwestcoeff-qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n]]=qnorthsouthcoeff  # Coeff for m,n+1
            elif n==Ncellz:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qwestcoeff-qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-2]]=qnorthsouthcoeff  # Coeff for m,n-1
            else:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qwestcoeff-2*qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-2]]=qnorthsouthcoeff  # Coeff for m,n-1
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n]]=qnorthsouthcoeff  # Coeff for m,n+1
    
        # z=0 boundary
        n=1
        for m in range(1,Ncellr+1):
            qeastcoeff=kr*2*np.pi*(m+0.5)*delz
            qwestcoeff=kr*2*np.pi*(m-0.5)*delz
            qnorthsouthcoeff=kz*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)/delz
            Estcoeff=rho*cp*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz/delt
            Convcoeff=hz0*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)
            Qgencoeff=qgen*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz
            Coeff[GridMap[m-1,n-1],GridMap[m-1,n]]=qnorthsouthcoeff+0.5*Convcoeff  # Coeff for m,n+1
            b[GridMap[m-1,n-1],0]=-Estcoeff*T[GridMap[m-1,n-1]]-Convcoeff*Tamb-Qgencoeff
            if m==1:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qeastcoeff-qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m,n-1]]=qeastcoeff  # Coeff for m+1,n
            elif m==Ncellr:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qwestcoeff-qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m-2,n-1]]=qwestcoeff  # Coeff for m-1,n
            else:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qwestcoeff-qeastcoeff-qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m,n-1]]=qeastcoeff  # Coeff for m+1,n
                Coeff[GridMap[m-1,n-1],GridMap[m-2,n-1]]=qwestcoeff  # Coeff for m-1,n
    
        # z=H boundary
        n=Ncellz
        for m in range(1,Ncellr+1):
            qeastcoeff=kr*2*np.pi*(m+0.5)*delz
            qwestcoeff=kr*2*np.pi*(m-0.5)*delz
            qnorthsouthcoeff=kz*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)/delz
            Estcoeff=rho*cp*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz/delt
            Convcoeff=hzH*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)
            Qgencoeff=qgen*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz
            Coeff[GridMap[m-1,n-1],GridMap[m-1,n-2]]=qnorthsouthcoeff+0.5*Convcoeff  # Coeff for m,n-1
            b[GridMap[m-1,n-1],0]=-Estcoeff*T[GridMap[m-1,n-1]]-Convcoeff*Tamb-Qgencoeff
            if m==1:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qeastcoeff-qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m,n-1]]=qeastcoeff  # Coeff for m+1,n
            elif m==Ncellr:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qwestcoeff-qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m-2,n-1]]=qwestcoeff  # Coeff for m-1,n
            else:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qeastcoeff-qwestcoeff-qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m,n-1]]=qeastcoeff  # Coeff for m+1,n
                Coeff[GridMap[m-1,n-1],GridMap[m-2,n-1]]=qwestcoeff  # Coeff for m-1,n

        #   print(Coeff)
        #   print(b)    
        T=np.linalg.solve(Coeff,b)
        T_all[:,i]=T.reshape((Ncellr*Ncellz,))
    print(T_all[:,i])
    return T_all

# This functions solves with Arrhenius generation term if you want to have Qgen due to Arrhenius only, give qgen=0
def Arrheniusgen(GridMap,kr,kz,hR,hz0,hzH,rho,cp,R,H,Ncellr,Ncellz,delr,delz,delt,Tinit,Tamb,tsample,qgen,Q0,E,Picard,Newton,Nmax,tol):
    Runicons=8.314 # Universal gas constant in J/molK
    Coeff=np.zeros((Ncellr*Ncellz,Ncellr*Ncellz))
    b=np.zeros((Ncellr*Ncellz,1))
    T_all=np.empty((Ncellr*Ncellz,tsample.shape[0]))
    T=Tinit*np.ones([Ncellr*Ncellz,1])
    T_all[:,0]=T.reshape(Ncellr*Ncellz,)
    for i in range(1,tsample.shape[0]):
        t=tsample[i]

# Filling coefficient and forcing matrix for all things other than generation          
# Coefficient and forcing matrix for interior nodes
        for m in range(2,Ncellr):
            qwestcoeff=kr*2*np.pi*(m-0.5)*delz
            qeastcoeff=kr*2*np.pi*(m+0.5)*delz
            qnorthsouthcoeff=kz*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)/delz
            Estcoeff=rho*cp*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz/delt
            Qgencoeff=qgen*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz
            for n in range(2,Ncellz):
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qeastcoeff-qwestcoeff-2*qnorthsouthcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m-2,n-1]]=qwestcoeff  # Coeff for m-1,n
                Coeff[GridMap[m-1,n-1],GridMap[m,n-1]]=qeastcoeff  # Coeff for m+1,n
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-2]]=qnorthsouthcoeff  # Coeff for m,n-1
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n]]=qnorthsouthcoeff  # Coeff for m,n+1
                b[GridMap[m-1,n-1],0]=-Estcoeff*T[GridMap[m-1,n-1]]-Qgencoeff

# Coefficient and forcing matrix for boundary nodes
# r=0 boundary
        m=1
        qeastcoeff=kr*2*np.pi*(m+0.5)*delz
        qnorthsouthcoeff=kz*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)/delz
        Estcoeff=rho*cp*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz/delt
        Qgencoeff=qgen*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz
        for n in range(1,Ncellz+1):
            Coeff[GridMap[m-1,n-1],GridMap[m,n-1]]=qeastcoeff  # Coeff for m+1,n
            b[GridMap[m-1,n-1],0]=-Estcoeff*T[GridMap[m-1,n-1]]-Qgencoeff
            if n==1:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qeastcoeff-qnorthsouthcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n]]=qnorthsouthcoeff  # Coeff for m,n+1
            elif n==Ncellz:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qeastcoeff-qnorthsouthcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-2]]=qnorthsouthcoeff  # Coeff for m,n-1
            else:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qeastcoeff-2*qnorthsouthcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-2]]=qnorthsouthcoeff  # Coeff for m,n-1
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n]]=qnorthsouthcoeff  # Coeff for m,n+1

        # r=R boundary
        m=Ncellr
        qwestcoeff=kr*2*np.pi*(m-0.5)*delz
        qnorthsouthcoeff=kz*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)/delz
        Estcoeff=rho*cp*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz/delt
        Convcoeff=hR*2*np.pi*R*delz
        Qgencoeff=qgen*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz
        for n in range(1,Ncellz+1):
            Coeff[GridMap[m-1,n-1],GridMap[m-2,n-1]]=qwestcoeff+0.5*Convcoeff  # Coeff for m-1,n
            b[GridMap[m-1,n-1],0]=-Estcoeff*T[GridMap[m-1,n-1]]-Convcoeff*Tamb-Qgencoeff
            if n==1:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qwestcoeff-qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n]]=qnorthsouthcoeff  # Coeff for m,n+1
            elif n==Ncellz:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qwestcoeff-qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-2]]=qnorthsouthcoeff  # Coeff for m,n-1
            else:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qwestcoeff-2*qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-2]]=qnorthsouthcoeff  # Coeff for m,n-1
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n]]=qnorthsouthcoeff  # Coeff for m,n+1
    
        # z=0 boundary
        n=1
        for m in range(1,Ncellr+1):
            qeastcoeff=kr*2*np.pi*(m+0.5)*delz
            qwestcoeff=kr*2*np.pi*(m-0.5)*delz
            qnorthsouthcoeff=kz*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)/delz
            Estcoeff=rho*cp*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz/delt
            Convcoeff=hz0*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)
            Qgencoeff=qgen*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz
            Coeff[GridMap[m-1,n-1],GridMap[m-1,n]]=qnorthsouthcoeff+0.5*Convcoeff  # Coeff for m,n+1
            b[GridMap[m-1,n-1],0]=-Estcoeff*T[GridMap[m-1,n-1]]-Convcoeff*Tamb-Qgencoeff
            if m==1:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qeastcoeff-qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m,n-1]]=qeastcoeff  # Coeff for m+1,n
            elif m==Ncellr:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qwestcoeff-qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m-2,n-1]]=qwestcoeff  # Coeff for m-1,n
            else:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qwestcoeff-qeastcoeff-qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m,n-1]]=qeastcoeff  # Coeff for m+1,n
                Coeff[GridMap[m-1,n-1],GridMap[m-2,n-1]]=qwestcoeff  # Coeff for m-1,n
    
        # z=H boundary
        n=Ncellz
        for m in range(1,Ncellr+1):
            qeastcoeff=kr*2*np.pi*(m+0.5)*delz
            qwestcoeff=kr*2*np.pi*(m-0.5)*delz
            qnorthsouthcoeff=kz*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)/delz
            Estcoeff=rho*cp*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz/delt
            Convcoeff=hzH*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)
            Qgencoeff=qgen*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz
            Coeff[GridMap[m-1,n-1],GridMap[m-1,n-2]]=qnorthsouthcoeff+0.5*Convcoeff  # Coeff for m,n-1
            b[GridMap[m-1,n-1],0]=-Estcoeff*T[GridMap[m-1,n-1]]-Convcoeff*Tamb-Qgencoeff
            if m==1:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qeastcoeff-qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m,n-1]]=qeastcoeff  # Coeff for m+1,n
            elif m==Ncellr:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qwestcoeff-qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m-2,n-1]]=qwestcoeff  # Coeff for m-1,n
            else:
                Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qeastcoeff-qwestcoeff-qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
                Coeff[GridMap[m-1,n-1],GridMap[m,n-1]]=qeastcoeff  # Coeff for m+1,n
                Coeff[GridMap[m-1,n-1],GridMap[m-2,n-1]]=qwestcoeff  # Coeff for m-1,n

        #   print(Coeff)
        #   print(b)
        
        # Here comes the treatment of non-linear Arrhenius generation

        # Using PICARD's method (Zeroth order)
        if Picard:
            for k in range(0,Nmax+1):
                if k==0:
                    Tk=T
                for m in range(1,Ncellr+1):
                    for n in range(1,Ncellz+1):
                        bk[GridMap[m-1,n-1],0]=b[GridMap[m-1,n-1],0]-Q0*np.exp(-E/Runicons/Tk[GridMap[m-1,n-1],0])
                Tknew=np.linalg.solve(Coeff,bk)
                Maxerr=np.abs(np.max(Tknew-Tk))
                if Maxerr>tol:
                    Tk=Tknew
                    if k==Nmax:
                        print("Maximum error %6.6f is higher than the tolerance after Nmax iterations!" %(Maxerr))
                else:
                    break
            T=Tknew
            T_all[:,i]=Tknew.reshape(Ncellr*Ncellz,)

        # Using Newton's method (First order)
        if Newton:
            for k in range(0,Nmax):
                if k==0:
                    Tk=T
                for m in range(1,Ncellr+1):
                    for n in range(1,Ncellz+1):
                        bk[GridMap[m-1,n-1],0]=b[GridMap[m-1,n-1],0]-Q0*np.exp(-E/Runicons/Tk[GridMap[m-1,n-1],0])+Q0*np.exp(-E/Runicons/Tk[GridMap[m-1,n-1],0])*(E/Runicons/(Tk[GridMap[m-1,n-1],0]**2))*Tk[GridMap[m-1,n-1],0]
                        Coeffk[GridMap[m-1,n-1],GridMap[m-1,n-1]]=Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]+Q0*np.exp(-E/Runicons/Tk[GridMap[m-1,n-1],0])*(E/Runicons/(Tk[GridMap[m-1,n-1],0]**2))
                Tknew=np.linalg.solve(Coeffk,bk)
                Maxerr=np.abs(np.max(Tknew-Tk))
                if Maxerr>tol:
                    Tk=Tknew
                    if k==Nmax:
                        print("Maximum error %6.6f is higher than the tolerance after Nmax iterations!" %(Maxerr))
                else:
                    break
            T=Tknew
            T_all[:,i]=Tknew.reshape(Ncellr*Ncellz,)
            print(t)    
    return T_all