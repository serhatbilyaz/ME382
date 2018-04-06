
# This is a code to calculate 2D axisymmetric heat transfer in cylindrical coordinates (r,z)

import numpy as np
#import scipy
#import matplotlib
#from matplotlib import pyplot
from numpy import pi
#from scipy import special
import finiteDifference
import analytical_nogen
import lineplot
import contourplot
import createanimation
#import pandas as pd


# Inputs
R=1.2 # Radius of the cylinder (m)
H=1 # Height of the cylinder (m)

Ncellr=12 # Number of cells in r direction
Ncellz=5 # Number of cells in z direction

delt=10 # Time step in seconds
tfinal=1000 # Final time in seconds

Tinit=400 # Uniform initial temperature 
Tamb=300 # Ambient temperature

qgen=0 # Volumetric heat generation (W/m^3)

kr=40 # Heat conductivity in r direction (W/m/K)
kz=40 # Heat conductivity in z direction (W/m/K)

hR=20 # Convection coefficient at r=R (W/m^2/K)
hzH=20 # Convection coefficient at z=H (W/m^2/K)
hz0=20 # Convection coefficient at z=0 (W/m^2/K)

cp=490 # Specific heat of the cylinder (J/kgK)
rho=8000 # Density of the cylinder (kg/m^3)
delr=R/Ncellr
delz=H/Ncellz

# Grid numbering (Row increase => r increase, Column increase => z increase)
GridMap=np.zeros((Ncellz,Ncellr),dtype=np.int)
for i in range(0,Ncellz):
    GridMap[i,:]=(i*Ncellr)+np.arange(0,Ncellr,1)
print(GridMap)

#FD solution uses transposed GridMap (Column increase => r increase, Row increase => z increase)
GridMapFD=GridMap
GridMapFD=np.transpose(GridMapFD)
print(GridMapFD)

# Locations of grid points (Consistent with original GridMap, r increases in row direction, z increases in column direction)
zLoc=np.zeros((Ncellz,Ncellr))
rLoc=np.zeros((Ncellz,Ncellr))
for m in range(0,Ncellz):
    for n in range(0,Ncellr):
        rLoc[m,n]=delr*(n+0.5)
        zLoc[m,n]=delz*(m+0.5)
print(zLoc)
print(rLoc)

# Sample the time domain
tsample=np.linspace(0,tfinal,tfinal/delt+1)
print(tsample)

# Initialize T_FD_all and T_analy_all
T_FD_all=np.empty((Ncellr*Ncellz,tsample.shape[0]))*np.nan
T_analy_all=np.empty((Ncellr*Ncellz,tsample.shape[0]))*np.nan

# Finite Difference Solution
T_FD_all=finiteDifference.uniformgen(GridMapFD,kr,kz,hR,hz0,hzH,rho,cp,R,H,Ncellr,Ncellz,delr,delz,delt,Tinit,Tamb,tsample,qgen)

# Analytical Solution
#T_analy_all=analytical_nogen.calc_whole(H,R,zLoc,rLoc,tsample,rho,cp,hzH,hR,kz,kr,Tinit,Tamb,GridMap,Ncellr,Ncellz)


tplot=tfinal
# Obtain contour subplots 
contourplot.compareanaly(T_FD_all,T_analy_all,rLoc,zLoc,tsample,delt,Ncellr,Ncellz,tplot)
# Plotting along r at certain z
zplot=H/2 # Plot location
lineplot.alongr_plotT(T_FD_all,T_analy_all,rLoc,zLoc,tsample,R,delz,delt,Ncellr,Ncellz,zplot,tplot)
# Plotting along z at certain r
rplot=R/2 # Plot location
lineplot.alongz_plotT(T_FD_all,T_analy_all,rLoc,zLoc,tsample,R,delz,delt,Ncellr,Ncellz,rplot,tplot)

# Create animation 
#createanimation.compareanaly(T_FD_all,T_analy_all,rLoc,zLoc,tsample,Ncellr,Ncellz)
