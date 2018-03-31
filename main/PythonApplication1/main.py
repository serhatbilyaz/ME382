
# This is a code to calculate 2D axisymmetric heat transfer in cylindrical coordinates (r,z)

import numpy as np
import matplotlib
from matplotlib import pyplot
from numpy import pi
#import pandas as pd

# Inputs
R=1 # Radius of the cylinder (m)
H=1 # Height of the cylinder (m)

Ncellr=10 # Number of cells in r direction
Ncellz=10 # Number of cells in z direction

delt=0.1 # Time step
tfinal=10 # Final time

Tinit=300 # Uniform initial temperature 
Tamb=300 # Ambient temperature

kr=0.25 # Heat conductivity in r direction (W/m/K)
kz=0.25 # Heat conductivity in z direction (W/m/K)

hR=10 # Convection coefficient at r=R (W/m^2/K)
hzH=10 # Convection coefficient at z=H (W/m^2/K)
hz0=10 # Convection coefficient at z=0 (W/m^2/K)

cp=1000 # Specific heat of the cylinder (J/kgK)
rho=1000 # Density of the cylinder (kg/m^3)
delr=R/Ncellr
delz=H/Ncellz

# Grid numbering
GridMap=np.zeros((Ncellr,Ncellz),dtype=np.int)
for i in range(0,Ncellz):
    GridMap[i,:]=(i*Ncellr)+np.arange(0,Ncellr,1)
print(GridMap)

# Initializing T
T=Tinit*np.ones([Ncellr*Ncellz,1])
#print(T)

# Solution
Coeff=np.zeros((Ncellr*Ncellz,Ncellr*Ncellz))
b=np.zeros((Ncellr*Ncellz,1))

# Coefficient and forcing matrix for interior nodes
for m in range(2,Ncellr-1):
    qwestcoeff=kr*2*np.pi*(m-0.5)*delz
    qeastcoeff=kr*2*np.pi*(m+0.5)*delz
    qnorthsouthcoeff=kz*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)/delz
    Estcoeff=rho*cp*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz/delt
    for n in range(2,Ncellz-1):
        Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=Estcoeff-qeastcoeff-qwestcoeff-2*qnorthsouthcoeff  # Coeff for m,n
        Coeff[GridMap[m-1,n-1],GridMap[m-2,n-1]]=qwestcoeff  # Coeff for m-1,n
        Coeff[GridMap[m-1,n-1],GridMap[m,n-1]]=qeastcoeff  # Coeff for m+1,n
        Coeff[GridMap[m-1,n-1],GridMap[m-1,n-2]]=qnorthsouthcoeff  # Coeff for m,n-1
        Coeff[GridMap[m-1,n-1],GridMap[m-1,n]]=qnorthsouthcoeff  # Coeff for m,n+1
        b[GridMap[m-1,n-1],0]=Estcoeff*T[GridMap[m-1,n-1]]

# Coefficient and forcing matrix for boundary nodes
# r=0 boundary
m=1
qeastcoeff=kr*2*np.pi*(m+0.5)*delz
qnorthsouthcoeff=kz*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)/delz
Estcoeff=rho*cp*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz/delt
for n in range(1,Ncellz):
    Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=Estcoeff-qeastcoeff-2*qnorthsouthcoeff  # Coeff for m,n
    Coeff[GridMap[m-1,n-1],GridMap[m,n-1]]=qeastcoeff  # Coeff for m+1,n
    Coeff[GridMap[m-1,n-1],GridMap[m-1,n-2]]=qnorthsouthcoeff  # Coeff for m,n-1
    Coeff[GridMap[m-1,n-1],GridMap[m-1,n]]=qnorthsouthcoeff  # Coeff for m,n+1
    b[GridMap[m-1,n-1],0]=Estcoeff*T[GridMap[m-1,n-1]]
# r=R boundary
m=Ncellr
qwestcoeff=kr*2*np.pi*(m-0.5)*delz
qnorthsouthcoeff=kz*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)/delz
Estcoeff=rho*cp*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz/delt
Convcoeff=hR*2*np.pi*R*delz
for n in range(1,Ncellz):
    Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qwestcoeff-2*qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
    Coeff[GridMap[m-1,n-1],GridMap[m-2,n-1]]=qwestcoeff+0.5*Convcoeff  # Coeff for m-1,n
    Coeff[GridMap[m-1,n-1],GridMap[m-1,n-2]]=qnorthsouthcoeff  # Coeff for m,n-1
    Coeff[GridMap[m-1,n-1],GridMap[m-1,n]]=qnorthsouthcoeff  # Coeff for m,n+1
    b[GridMap[m-1,n-1],0]=-Estcoeff*T[GridMap[m-1,n-1]]-Convcoeff*Tamb
# z=0 boundary
n=1
for m in range(1,Ncellr):
    qeastcoeff=kr*2*np.pi*(m+0.5)*delz
    qwestcoeff=kr*2*np.pi*(m-0.5)*delz
    qnorthsouthcoeff=kz*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)/delz
    Estcoeff=rho*cp*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz/delt
    Convcoeff=hz0*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)

    Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qwestcoeff-qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
    Coeff[GridMap[m-1,n-1],GridMap[m,n-1]]=qeastcoeff  # Coeff for m+1,n
    Coeff[GridMap[m-1,n-1],GridMap[m-2,n-1]]=qwestcoeff  # Coeff for m-1,n
    Coeff[GridMap[m-1,n-1],GridMap[m-1,n]]=qnorthsouthcoeff+0.5*Convcoeff  # Coeff for m,n+1
    b[GridMap[m-1,n-1],0]=-Estcoeff*T[GridMap[m-1,n-1]]-Convcoeff*Tamb
# z=H boundary
n=Ncellz
for m in range(1,Ncellr):
    qeastcoeff=kr*2*np.pi*(m+0.5)*delz
    qwestcoeff=kr*2*np.pi*(m-0.5)*delz
    qnorthsouthcoeff=kz*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)/delz
    Estcoeff=rho*cp*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)*delz/delt
    Convcoeff=hzH*np.pi*((m+0.5)**2-(m-0.5)**2)*(delr**2)

    Coeff[GridMap[m-1,n-1],GridMap[m-1,n-1]]=-Estcoeff-qwestcoeff-qnorthsouthcoeff-1.5*Convcoeff  # Coeff for m,n
    Coeff[GridMap[m-1,n-1],GridMap[m,n-1]]=qeastcoeff  # Coeff for m+1,n
    Coeff[GridMap[m-1,n-1],GridMap[m-2,n-1]]=qwestcoeff  # Coeff for m-1,n
    Coeff[GridMap[m-1,n-1],GridMap[m-1,n-2]]=qnorthsouthcoeff+0.5*Convcoeff  # Coeff for m,n-1
    b[GridMap[m-1,n-1],0]=-Estcoeff*T[GridMap[m-1,n-1]]-Convcoeff*Tamb

#   print(Coeff)
#   print(b)