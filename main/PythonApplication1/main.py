
# This is a code to calculate 2D axisymmetric heat transfer in cylindrical coordinates (r,z)

import numpy as np
#import scipy
#import matplotlib
#from matplotlib import pyplot
from numpy import pi
#from scipy import special
import analytical_nogen
import lineplot
import contourplot
#import pandas as pd

# Inputs
R=1.2 # Radius of the cylinder (m)
H=1 # Height of the cylinder (m)

Ncellr=10 # Number of cells in r direction
Ncellz=10 # Number of cells in z direction

delt=10 # Time step in seconds
tfinal=1000 # Final time in seconds

Tinit=400 # Uniform initial temperature 
Tamb=300 # Ambient temperature

kr=40 # Heat conductivity in r direction (W/m/K)
kz=40 # Heat conductivity in z direction (W/m/K)

hR=20 # Convection coefficient at r=R (W/m^2/K)
hzH=20 # Convection coefficient at z=H (W/m^2/K)
hz0=20 # Convection coefficient at z=0 (W/m^2/K)

cp=490 # Specific heat of the cylinder (J/kgK)
rho=8000 # Density of the cylinder (kg/m^3)
delr=R/Ncellr
delz=H/Ncellz

# Grid numbering
GridMap=np.zeros((Ncellr,Ncellz),dtype=np.int)
for i in range(0,Ncellz):
    GridMap[i,:]=(i*Ncellr)+np.arange(0,Ncellr,1)
print(GridMap)

# Locations of grid points
zLoc=np.zeros((Ncellr,Ncellz))
rLoc=np.zeros((Ncellr,Ncellz))
for m in range(0,Ncellr):
    for n in range(0,Ncellz):
        rLoc[m,n]=delr*(m+0.5)
        zLoc[m,n]=delz*(n+0.5)
print(zLoc)
print(rLoc)
# Initializing T
T=Tinit*np.ones([Ncellr*Ncellz,1])
Tlocal=Tinit*np.ones([Ncellr,Ncellz])
#print(T)

# Solution
Coeff=np.zeros((Ncellr*Ncellz,Ncellr*Ncellz))
b=np.zeros((Ncellr*Ncellz,1))
tsample=np.linspace(0,tfinal,tfinal/delt)
for t in tsample:

# Coefficient and forcing matrix for interior nodes
    for m in range(2,Ncellr):
        qwestcoeff=kr*2*np.pi*(m-0.5)*delz
        qeastcoeff=kr*2*np.pi*(m+0.5)*delz
        qnorthsouthcoeff=kz*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)/delz
        Estcoeff=rho*cp*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz/delt
        for n in range(2,Ncellz):
            Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qeastcoeff-qwestcoeff-2*qnorthsouthcoeff  # Coeff for m,n
            Coeff[GridMap[m-1,n-1],GridMap[m-2,n-1]]=qwestcoeff  # Coeff for m-1,n
            Coeff[GridMap[m-1,n-1],GridMap[m,n-1]]=qeastcoeff  # Coeff for m+1,n
            Coeff[GridMap[m-1,n-1],GridMap[m-1,n-2]]=qnorthsouthcoeff  # Coeff for m,n-1
            Coeff[GridMap[m-1,n-1],GridMap[m-1,n]]=qnorthsouthcoeff  # Coeff for m,n+1
            b[GridMap[m-1,n-1],0]=-Estcoeff*T[GridMap[m-1,n-1]]

# Coefficient and forcing matrix for boundary nodes
# r=0 boundary
    m=1
    qeastcoeff=kr*2*np.pi*(m+0.5)*delz
    qnorthsouthcoeff=kz*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)/delz
    Estcoeff=rho*cp*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz/delt
    for n in range(1,Ncellz+1):
        Coeff[GridMap[m-1,n-1],GridMap[m,n-1]]=qeastcoeff  # Coeff for m+1,n
        b[GridMap[m-1,n-1],0]=-Estcoeff*T[GridMap[m-1,n-1]]
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
    for n in range(1,Ncellz+1):
        Coeff[GridMap[m-1,n-1],GridMap[m-2,n-1]]=qwestcoeff+0.5*Convcoeff  # Coeff for m-1,n
        b[GridMap[m-1,n-1],0]=-Estcoeff*T[GridMap[m-1,n-1]]-Convcoeff*Tamb
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

        Coeff[GridMap[m-1,n-1],GridMap[m-1,n]]=qnorthsouthcoeff+0.5*Convcoeff  # Coeff for m,n+1
        b[GridMap[m-1,n-1],0]=-Estcoeff*T[GridMap[m-1,n-1]]-Convcoeff*Tamb
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

        Coeff[GridMap[m-1,n-1],GridMap[m-1,n-2]]=qnorthsouthcoeff+0.5*Convcoeff  # Coeff for m,n-1
        b[GridMap[m-1,n-1],0]=-Estcoeff*T[GridMap[m-1,n-1]]-Convcoeff*Tamb
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
    for m in range(1,Ncellr+1):
        for n in range(1,Ncellz+1):
            Tlocal[n-1,m-1]=T[GridMap[m-1,n-1],0]
    print(t)
    #print(np.flip(Tlocal,0))


#Analytical solution (Uncomment if you want to solve the entire domain)

Theta=np.zeros((Ncellr,Ncellz))
for m in range(0,Ncellr):
    for n in range(0,Ncellz):        
        Theta[n,m]=analytical_nogen.calc_nogen(H,R,zLoc[m,n],rLoc[m,n],tfinal,rho,cp,hzH,hR,kz,kr)
print("Theta is",Theta)
Temp_analytical=Theta*(Tinit-Tamb)+Tamb
print("Analytical solution is:")
print(np.flip(Temp_analytical,0))

print("Numerical solution is:")
print(np.flip(Tlocal,0))


# Obtain contour subplots 
contourplot.compareanaly(Tlocal,Temp_analytical,rLoc,zLoc,tfinal)


# Plotting along r at certain z

z=H/2; # Plot location in z
n=(z/delz)-0.5 #Temp index in terms of z
n=int(np.floor(n)) # Find the closest integer to be an index
rplot_alongr=rLoc[:,n] #r locations for plotting
zplot_alongr=zLoc[0,n] #Evaluate real r (According to rounded index)
#print(rplot_alongr[3])
Temp_FD_alongr=Tlocal[n,:] 
Temp_analy_alongr=np.zeros((Ncellr))
for m in range(0,Ncellr):
    Temp_analy_alongr[m]=analytical_nogen.calc_nogen(H,R,zplot_alongr,rplot_alongr[m],tfinal,rho,cp,hzH,hR,kz,kr)
Temp_analy_alongr=Temp_analy_alongr*(Tinit-Tamb)+Tamb
#Temp_analy_alongr=Temp_analytical[n,:] # Use if you want to extract it from the entire domain solution (see above)

lineplot.alongr_plotT(Temp_FD_alongr,Temp_analy_alongr,rplot_alongr,zplot_alongr,tfinal,R)

# Plotting along z at certain r
r=R/2 # Plot location in r
m=(r/delr)-0.5 #Temp index in terms of r
m=int(np.floor(m)) # Find the closest integer to be an index
rplot_alongz=rLoc[m,0] # Evaluate real r (According to rounded index)
zplot_alongz=zLoc[m,:] # z locations for plotting

Temp_FD_alongz=Tlocal[:,m]
Temp_analy_alongz=np.zeros((Ncellz))
for n in range(0,Ncellz):
    Temp_analy_alongz[n]=analytical_nogen.calc_nogen(H,R,zplot_alongz[n],rplot_alongz,tfinal,rho,cp,hzH,hR,kz,kr)
Temp_analy_alongz=Temp_analy_alongz*(Tinit-Tamb)+Tamb

#Temp_analy_alongz=Temp_analytical[:,m] # Use if you want to extract it from the entire domain solution (see above)

lineplot.alongz_plotT(Temp_FD_alongz,Temp_analy_alongz,rplot_alongz,zplot_alongz,tfinal,H)