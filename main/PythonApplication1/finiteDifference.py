#This script contains all subroutines that solve the problem with finite
#difference approach
import numpy as np
#import scipy
from numpy import pi
#from scipy import special
from scipy import linalg
from scipy.sparse.linalg import gmres
#from scipy.sparse.linalg import bicg
#from scipy.sparse.linalg import bicgstab
from scipy.sparse.linalg import inv
from scipy.sparse.linalg import LinearOperator, spilu

def uniformgen(GridMap,kr,kz,hR,hz0,hzH,rho,cp,R,H,Ncellr,Ncellz,delr,delz,delt,Tinit,Tamb,tsample,qgen,A_antoine,B_antoine,C_antoine,Mgas,P0, GM, tolGMRES):
    Runicons = 8.314 # Universal gas constant in J/molK
    Vol_T = (np.pi * R ** 2) * H
    mtot = rho * Vol_T
    Pvap_all = np.empty((tsample.shape[0],))
    mG_all = np.empty((tsample.shape[0],))
    Tavg_all = np.empty((tsample.shape[0],))
    Tavg_all[0] = Tinit    
     
    Coeff = np.zeros((Ncellr * Ncellz,Ncellr * Ncellz))
    b = np.zeros((Ncellr * Ncellz,1))
    T_all = np.empty((Ncellr * Ncellz,tsample.shape[0]))
    T = Tinit * np.ones([Ncellr * Ncellz,1])
    T_all[:,0] = T.reshape(Ncellr * Ncellz,)

    # Building the coefficient matrix out of the time loop
    # Coefficient and forcing matrix for interior nodes
    for m in range(2,Ncellr):
        qwestcoeff = kr * 2 * np.pi * (m - 0.5) * delz
        qeastcoeff = kr * 2 * np.pi * (m + 0.5) * delz
        qnorthsouthcoeff = kz * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) / delz
        Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
        Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
        for n in range(2,Ncellz):
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qeastcoeff - qwestcoeff - 2 * qnorthsouthcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 2,n - 1]] = qwestcoeff  # Coeff for m-1,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m,n - 1]] = qeastcoeff  # Coeff for m+1,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 2]] = qnorthsouthcoeff  # Coeff for m,n-1
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n]] = qnorthsouthcoeff  # Coeff for m,n+1
            

# Coefficient and forcing matrix for boundary nodes
# r=0 boundary
    m = 1
    qeastcoeff = kr * 2 * np.pi * (m + 0.5) * delz
    qnorthsouthcoeff = kz * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) / delz
    Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
    Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
    for n in range(1,Ncellz + 1):
        Coeff[GridMap[m - 1,n - 1],GridMap[m,n - 1]] = qeastcoeff  # Coeff for m+1,n
        
        if n == 1:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qeastcoeff - qnorthsouthcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n]] = qnorthsouthcoeff  # Coeff for m,n+1
        elif n == Ncellz:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qeastcoeff - qnorthsouthcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 2]] = qnorthsouthcoeff  # Coeff for m,n-1
        else:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qeastcoeff - 2 * qnorthsouthcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 2]] = qnorthsouthcoeff  # Coeff for m,n-1
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n]] = qnorthsouthcoeff  # Coeff for m,n+1

    # r=R boundary
    m = Ncellr
    qwestcoeff = kr * 2 * np.pi * (m - 0.5) * delz
    qnorthsouthcoeff = kz * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) / delz
    Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
    Convcoeff = hR * 2 * np.pi * R * delz
    Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
    for n in range(1,Ncellz + 1):
        Coeff[GridMap[m - 1,n - 1],GridMap[m - 2,n - 1]] = qwestcoeff + 0.5 * Convcoeff  # Coeff for m-1,n
        
        if n == 1:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qwestcoeff - qnorthsouthcoeff - 1.5 * Convcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n]] = qnorthsouthcoeff  # Coeff for m,n+1
        elif n == Ncellz:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qwestcoeff - qnorthsouthcoeff - 1.5 * Convcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 2]] = qnorthsouthcoeff  # Coeff for m,n-1
        else:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qwestcoeff - 2 * qnorthsouthcoeff - 1.5 * Convcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 2]] = qnorthsouthcoeff  # Coeff for m,n-1
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n]] = qnorthsouthcoeff  # Coeff for m,n+1
    
    # z=0 boundary
    n = 1
    for m in range(1,Ncellr + 1):
        qeastcoeff = kr * 2 * np.pi * (m + 0.5) * delz
        qwestcoeff = kr * 2 * np.pi * (m - 0.5) * delz
        qnorthsouthcoeff = kz * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) / delz
        Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
        Convcoeff = hz0 * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2)
        Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
        Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n]] = qnorthsouthcoeff + 0.5 * Convcoeff  # Coeff for m,n+1
        
        if m == 1:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qeastcoeff - qnorthsouthcoeff - 1.5 * Convcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m,n - 1]] = qeastcoeff  # Coeff for m+1,n
        elif m == Ncellr:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qwestcoeff - qnorthsouthcoeff - 1.5 * Convcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 2,n - 1]] = qwestcoeff  # Coeff for m-1,n
        else:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qwestcoeff - qeastcoeff - qnorthsouthcoeff - 1.5 * Convcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m,n - 1]] = qeastcoeff  # Coeff for m+1,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 2,n - 1]] = qwestcoeff  # Coeff for m-1,n
    
    # z=H boundary
    n = Ncellz
    for m in range(1,Ncellr + 1):
        qeastcoeff = kr * 2 * np.pi * (m + 0.5) * delz
        qwestcoeff = kr * 2 * np.pi * (m - 0.5) * delz
        qnorthsouthcoeff = kz * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) / delz
        Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
        Convcoeff = hzH * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2)
        Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
        Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 2]] = qnorthsouthcoeff + 0.5 * Convcoeff  # Coeff for m,n-1
        
        if m == 1:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qeastcoeff - qnorthsouthcoeff - 1.5 * Convcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m,n - 1]] = qeastcoeff  # Coeff for m+1,n
        elif m == Ncellr:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qwestcoeff - qnorthsouthcoeff - 1.5 * Convcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 2,n - 1]] = qwestcoeff  # Coeff for m-1,n
        else:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qeastcoeff - qwestcoeff - qnorthsouthcoeff - 1.5 * Convcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m,n - 1]] = qeastcoeff  # Coeff for m+1,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 2,n - 1]] = qwestcoeff  # Coeff for m-1,n
    if GM==0:
        #invCoeff=inv(Coeff)
        invCoeff=np.linalg.inv(Coeff)


    for i in range(1,tsample.shape[0]):
        t = tsample[i]
        print("t_FD=%4.4f" % (t))  

        for m in range(2,Ncellr):
            Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
            Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
            for n in range(2,Ncellz):
                b[GridMap[m - 1,n - 1],0] = -Estcoeff * T[GridMap[m - 1,n - 1]] - Qgencoeff
        # r=0 boundary
        m = 1
        Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
        Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
        for n in range(1,Ncellz + 1):
            b[GridMap[m - 1,n - 1],0] = -Estcoeff * T[GridMap[m - 1,n - 1]] - Qgencoeff
        # r=R boundary
        m = Ncellr
        Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
        Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
        Convcoeff = hR * 2 * np.pi * R * delz
        for n in range(1,Ncellz + 1):
            b[GridMap[m - 1,n - 1],0] = -Estcoeff * T[GridMap[m - 1,n - 1]] - Convcoeff * Tamb - Qgencoeff
        # z=0 boundary
        n = 1
        for m in range(1,Ncellr + 1):
            Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
            Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
            Convcoeff = hzH * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2)
            b[GridMap[m - 1,n - 1],0] = -Estcoeff * T[GridMap[m - 1,n - 1]] - Convcoeff * Tamb - Qgencoeff
        # z=H boundary
        n = Ncellz
        for m in range(1,Ncellr + 1):
            Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
            Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
            Convcoeff = hzH * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2)
            b[GridMap[m - 1,n - 1],0] = -Estcoeff * T[GridMap[m - 1,n - 1]] - Convcoeff * Tamb - Qgencoeff
        
        #   print(Coeff)
        #print("Condition number of Coeff is %8.8f" %(np.linalg.cond(Coeff)))
        #   print(b)
        print(np.shape(Coeff))
        if GM == 0:
            #T = np.linalg.solve(Coeff,b)
            T=np.matmul(invCoeff,b)
        else:
            ilu = spilu(Coeff)
            Mx = lambda x: ilu.solve(x)
            M = LinearOperator((Ncellr * Ncellz,Ncellr * Ncellz),Mx)
            Tgmres = gmres(Coeff, b, x0=None, tol=tolGMRES,M=M)
            #Tgmres=bicg(Coeff, b, x0=None, tol=tol)
            #Tgmres=bicgstab(Coeff, b, x0=None, tol=tol)
            T = np.array(Tgmres[0])
            gmresinfo = Tgmres[1]

        Tavg = np.mean(T)
        Tavg_all[i] = Tavg
        T_all[:,i] = T.reshape((Ncellr * Ncellz,))
    #print(T_all[:,i])
    return (T_all,Tavg_all)

# This functions solves with Arrhenius generation term if you want to have Qgen
# due to Arrhenius only, give qgen=0
def Arrheniusgen(GridMap,kr,kz,hR,hz0,hzH,rho,cp,R,H,Ncellr,Ncellz,delr,delz,delt,Tinit,Tamb,tsample,qgen,Q0,E,Picard,Newton,Nmax,tolPic,tolGMRES,A_antoine,B_antoine,C_antoine,Mgas,P0,GM):
    Runicons = 8.314 # Universal gas constant in J/molK
    Vol_T = (np.pi * R ** 2) * H
    mtot = rho * Vol_T

    Coeff = np.zeros((Ncellr * Ncellz,Ncellr * Ncellz))
    Coeffk = np.zeros((Ncellr * Ncellz,Ncellr * Ncellz))
    Hf = np.zeros((Ncellr * Ncellz,Ncellr * Ncellz))
    b = np.zeros((Ncellr * Ncellz,1))
    bk = np.zeros((Ncellr * Ncellz,1))
    bknew = np.zeros((Ncellr * Ncellz,1))
    step = np.zeros((Ncellr * Ncellz,1))
    T_all = np.empty((Ncellr * Ncellz,tsample.shape[0])) * np.nan
    #Pvap_all = np.empty((tsample.shape[0],))
    #mG_all = np.empty((tsample.shape[0],))
    #VolL_all = np.empty((tsample.shape[0],))
    Tavg_all = np.empty((tsample.shape[0],)) * np.nan
    Tavg_all[0] = Tinit
    T = Tinit * np.ones([Ncellr * Ncellz,1])
    T_all[:,0] = T.reshape(Ncellr * Ncellz,)
# Filling coefficient and forcing matrix for all things other than generation
# Coefficient and forcing matrix for interior nodes
    for m in range(2,Ncellr):
        qwestcoeff = kr * 2 * np.pi * (m - 0.5) * delz
        qeastcoeff = kr * 2 * np.pi * (m + 0.5) * delz
        qnorthsouthcoeff = kz * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) / delz
        Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
        Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
        for n in range(2,Ncellz):
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qeastcoeff - qwestcoeff - 2 * qnorthsouthcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 2,n - 1]] = qwestcoeff  # Coeff for m-1,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m,n - 1]] = qeastcoeff  # Coeff for m+1,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 2]] = qnorthsouthcoeff  # Coeff for m,n-1
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n]] = qnorthsouthcoeff  # Coeff for m,n+1
            

# Coefficient and forcing matrix for boundary nodes
# r=0 boundary
    m = 1
    qeastcoeff = kr * 2 * np.pi * (m + 0.5) * delz
    qnorthsouthcoeff = kz * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) / delz
    Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
    Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
    for n in range(1,Ncellz + 1):
        Coeff[GridMap[m - 1,n - 1],GridMap[m,n - 1]] = qeastcoeff  # Coeff for m+1,n
        if n == 1:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qeastcoeff - qnorthsouthcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n]] = qnorthsouthcoeff  # Coeff for m,n+1
        elif n == Ncellz:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qeastcoeff - qnorthsouthcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 2]] = qnorthsouthcoeff  # Coeff for m,n-1
        else:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qeastcoeff - 2 * qnorthsouthcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 2]] = qnorthsouthcoeff  # Coeff for m,n-1
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n]] = qnorthsouthcoeff  # Coeff for m,n+1

    # r=R boundary
    m = Ncellr
    qwestcoeff = kr * 2 * np.pi * (m - 0.5) * delz
    qnorthsouthcoeff = kz * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) / delz
    Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
    Convcoeff = hR * 2 * np.pi * R * delz
    Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
    for n in range(1,Ncellz + 1):
        Coeff[GridMap[m - 1,n - 1],GridMap[m - 2,n - 1]] = qwestcoeff + 0.5 * Convcoeff  # Coeff for m-1,n
        if n == 1:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qwestcoeff - qnorthsouthcoeff - 1.5 * Convcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n]] = qnorthsouthcoeff  # Coeff for m,n+1
        elif n == Ncellz:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qwestcoeff - qnorthsouthcoeff - 1.5 * Convcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 2]] = qnorthsouthcoeff  # Coeff for m,n-1
        else:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qwestcoeff - 2 * qnorthsouthcoeff - 1.5 * Convcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 2]] = qnorthsouthcoeff  # Coeff for m,n-1
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n]] = qnorthsouthcoeff  # Coeff for m,n+1
    
    # z=0 boundary
    n = 1
    for m in range(1,Ncellr + 1):
        qeastcoeff = kr * 2 * np.pi * (m + 0.5) * delz
        qwestcoeff = kr * 2 * np.pi * (m - 0.5) * delz
        qnorthsouthcoeff = kz * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) / delz
        Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
        Convcoeff = hz0 * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2)
        Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
        Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n]] = qnorthsouthcoeff + 0.5 * Convcoeff  # Coeff for m,n+1
        if m == 1:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qeastcoeff - qnorthsouthcoeff - 1.5 * Convcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m,n - 1]] = qeastcoeff  # Coeff for m+1,n
        elif m == Ncellr:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qwestcoeff - qnorthsouthcoeff - 1.5 * Convcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 2,n - 1]] = qwestcoeff  # Coeff for m-1,n
        else:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qwestcoeff - qeastcoeff - qnorthsouthcoeff - 1.5 * Convcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m,n - 1]] = qeastcoeff  # Coeff for m+1,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 2,n - 1]] = qwestcoeff  # Coeff for m-1,n
    
    # z=H boundary
    n = Ncellz
    for m in range(1,Ncellr + 1):
        qeastcoeff = kr * 2 * np.pi * (m + 0.5) * delz
        qwestcoeff = kr * 2 * np.pi * (m - 0.5) * delz
        qnorthsouthcoeff = kz * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) / delz
        Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
        Convcoeff = hzH * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2)
        Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
        Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 2]] = qnorthsouthcoeff + 0.5 * Convcoeff  # Coeff for m,n-1
        if m == 1:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qeastcoeff - qnorthsouthcoeff - 1.5 * Convcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m,n - 1]] = qeastcoeff  # Coeff for m+1,n
        elif m == Ncellr:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qwestcoeff - qnorthsouthcoeff - 1.5 * Convcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 2,n - 1]] = qwestcoeff  # Coeff for m-1,n
        else:
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 1,n - 1]] = -Estcoeff - qeastcoeff - qwestcoeff - qnorthsouthcoeff - 1.5 * Convcoeff  # Coeff for m,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m,n - 1]] = qeastcoeff  # Coeff for m+1,n
            Coeff[GridMap[m - 1,n - 1],GridMap[m - 2,n - 1]] = qwestcoeff  # Coeff for m-1,n    
        #   print(Coeff)
        #   print(b)
        #print("Condition number of Coeff is %8.8f" %(np.linalg.cond(Coeff)))
    print(np.shape(Coeff))
    if GM==0:
        #invCoeff=inv(Coeff)
        invCoeff=np.linalg.inv(Coeff)


    for i in range(1,tsample.shape[0]):
        t = tsample[i]
        print("t_FD=%4.4f" % (t)) 

    # Build b vector without taking care of the Arrhenius generation
        for m in range(2,Ncellr):
            Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
            Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
            for n in range(2,Ncellz):
                b[GridMap[m - 1,n - 1],0] = -Estcoeff * T[GridMap[m - 1,n - 1]] - Qgencoeff
        # r=0 boundary
        m = 1
        Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
        Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
        for n in range(1,Ncellz + 1):
            b[GridMap[m - 1,n - 1],0] = -Estcoeff * T[GridMap[m - 1,n - 1]] - Qgencoeff
        # r=R boundary
        m = Ncellr
        Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
        Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
        Convcoeff = hR * 2 * np.pi * R * delz
        for n in range(1,Ncellz + 1):
            b[GridMap[m - 1,n - 1],0] = -Estcoeff * T[GridMap[m - 1,n - 1]] - Convcoeff * Tamb - Qgencoeff
        # z=0 boundary
        n = 1
        for m in range(1,Ncellr + 1):
            Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
            Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
            Convcoeff = hzH * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2)
            b[GridMap[m - 1,n - 1],0] = -Estcoeff * T[GridMap[m - 1,n - 1]] - Convcoeff * Tamb - Qgencoeff
        # z=H boundary
        n = Ncellz
        for m in range(1,Ncellr + 1):
            Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
            Qgencoeff = qgen * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
            Convcoeff = hzH * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2)
            b[GridMap[m - 1,n - 1],0] = -Estcoeff * T[GridMap[m - 1,n - 1]] - Convcoeff * Tamb - Qgencoeff

        # Here comes the treatment of non-linear Arrhenius generation

        # Using PICARD's method (Zeroth order)
        if Picard:
            for k in range(0,Nmax + 1):
                if k == 0:
                    Tk = T
                for m in range(1,Ncellr + 1):
                    Estcoeff = rho * cp * np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz / delt
                    Vol = np.pi * ((m + 0.5) ** 2 - (m - 0.5) ** 2) * (delr ** 2) * delz
                    for n in range(1,Ncellz + 1):
                        C = Q0 * np.exp(-E / Runicons / Tk[GridMap[m - 1,n - 1],0])
                        bk[GridMap[m - 1,n - 1],0] = b[GridMap[m - 1,n - 1],0] - C * Vol

                #print("Condition number of Coeffk is %8.8f"
                #%(np.linalg.cond(Coeffk)))
                #print(Coeff)
                #print(bk)
                #TknewLU=np.linalg.solve(Coeff,bk)
                if GM == 0:
                    #Tknew = np.linalg.solve(Coeff,bk)
                    Tknew=np.matmul(invCoeff,bk)
                else:
                    ilu = spilu(Coeff)
                    Mx = lambda x: ilu.solve(x)
                    M = LinearOperator((Ncellr * Ncellz,Ncellr * Ncellz),Mx)
                    Tgmres = gmres(Coeff, b, x0=None, tol=tolGMRES,M=M)
                    Tknew = np.reshape(np.array(Tgmres[0]),(Ncellr * Ncellz,1))
                    gmresinfo = Tgmres[1]
                
                #step=Tknew-Tk
                #for linesearchstep in range(1,2):
                #for m in range(1,Ncellr+1):
                #   for n in range(1,Ncellz+1):
                #     C=Q0*np.exp(-E/Runicons/Tknew[GridMap[m-1,n-1],0])
                #     bknew[GridMap[m-1,n-1],0]=b[GridMap[m-1,n-1],0]-C*Vol
                        
                #f1=np.matmul(Coeff,Tk)-bk
                #f2=np.matmul(Coeff,Tknew)-bknew
                #fdif=f2-f1
                #print("Maximum of f2-f1 is %6.6f"%(np.max(fdif)))

                Maxerr = np.max(np.abs(Tknew - Tk))
                #print("Maxerr is %6.6f at %4.0f th iteration" %(Maxerr,k))
                step = (Tknew - Tk)
                if Maxerr > tolPic:
                    Tk = Tk + step
                    if k == Nmax:
                        print("Maximum error %6.6f is higher than the tolerance after Nmax iterations!" % (Maxerr))
                    if Maxerr > 10 ** 6:
                        print("Picard iteration diverged at t=%4.4f" % (t))
                        return (T_all,Tavg_all)
                else:
                    break
        
        T = Tknew
        #print(T)

        Tavg = np.mean(T)
        if np.isnan(Tavg):
            print("NaN Tavg is found")
        Tavg_all[i] = Tavg
        T_all[:,i] = Tknew.reshape(Ncellr * Ncellz,)
                
    return (T_all,Tavg_all)