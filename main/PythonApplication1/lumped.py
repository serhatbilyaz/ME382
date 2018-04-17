#This code solves the cylinder heat transfer problem with lumped assumption

import numpy as np

def calc_uniformgen(rho,cp,R,H,hR,hz0,hzH,Tamb,Tinit,delt,tsample,qgen):
    Runicon=8.314
    Vol=H*np.pi*(R**2)
    Aside=2*np.pi*R*H
    Atop=np.pi*(R**2)
    Abottom=Atop

    Tk=Tinit
    Tall=np.zeros(tsample.shape[0])
    Tall[0]=Tinit
    Estcoeff=rho*Vol*cp/delt
    Convcoeff=hz0*Abottom+hR*Aside+hzH*Atop
    Coeff=-Estcoeff-Convcoeff
    for i in range(1,tsample.shape[0]):
        t=tsample(i)
        b=-Convcoeff*Tamb-Estcoeff*Tall(i-1)
        C=qgen
        b=b-C*Vol
        Tknew=np.linalg.solve(Coeff,b)
        if Tknew>1e5:
            break
        Tall[i]=Tknew
    return Tall


def calc_Arrheniusgen(rho,cp,R,H,hR,hz0,hzH,Tamb,Tinit,delt,tsample,Q0,Ea,Picard,Nmax,tol):
    Runicon=8.314
    Vol=H*np.pi*(R**2)
    Aside=2*np.pi*R*H
    Atop=np.pi*(R**2)
    Abottom=Atop

    Tk=Tinit
    Tall=np.zeros(tsample.shape[0])
    Tall[0]=Tinit
    Estcoeff=rho*Vol*cp/delt
    Convcoeff=hz0*Abottom+hR*Aside+hzH*Atop
    Coeff=-Estcoeff-Convcoeff
    for i in range(1,tsample.shape[0]):
        t=tsample(i)
        for k in range(1,Nmax):
            b=-Convcoeff*Tamb-Estcoeff*Tall(i-1)
            C=Q0*np.exp(-Ea/Runicon/Tk)
            b=b-C*Vol
            Tknew=np.linalg.solve(Coeff,b)
            Res=Tknew-Tk
            if Res>tol:
                Tk=Tknew
            else:
                break
            if Tknew>1e5:
                break
        Tall[i]=Tknew
    return Tall









