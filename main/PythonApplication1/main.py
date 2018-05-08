
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
#import contourplot
#import createanimation
import pandas as pd
import csv
import time

Name='AG2'
print('%s'%(Name))
# Inputs
R=0.009 # Radius of the cylinder (m)
H=0.065 # Height of the cylinder (m)

Ncellr=60 # Number of cells in r direction
Ncellz=60 # Number of cells in z direction

delt=0.5 # Time step in seconds
tfinal=50*60 # Final time in seconds

Tinit=273+25 # Uniform initial temperature 
Tamb=273+155 # Ambient temperature

qgen=0*10**6 # Volumetric heat generation (W/m^3)

kr=5 # Heat conductivity in r direction (W/m/K)
kz=5 # Heat conductivity in z direction (W/m/K)

hR=20 # Convection coefficient at r=R (W/m^2/K)
hzH=20 # Convection coefficient at z=H (W/m^2/K)
hz0=20 # Convection coefficient at z=0 (W/m^2/K)

cp=1040 # Specific heat of the cylinder (J/kgK)
rho=2843 # Density of the cylinder (kg/m^3)

# Inputs for Arrhenius gen solution
Ea=140000 # Activation energy (J/mol)
Q0=1*10**22  # Arrhenius generation term constant
# Selection of nonlinear iterative solution method and its attributes
Picard=1# Use Picard's method
Newton=0  # Use Newton's method (Actually, don't use it now, it is not stable)
Nmax=10000  # Max number of iterations
tol=1*10**-4 # Accepted tolerance, also used for GMRES
GM=1    #if GM=0, LU solver, if GM=1 gmres solver

#Inputs for vent gas release
P0=100000 # Background pressure in Pa
Mgas=0.1181311 # Molecular weight of the vent gas (kg/mol)
#Mmix=28 # Molecular weight of the mixture (e.g. air in the background)
#Antoine equation constants
A_antoine=4.77616
B_antoine=1721.904 
C_antoine=-37.959


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

# Read experimental data
with open('expdatacsv.csv', newline='') as csvfile:
    Expdata=pd.read_csv(csvfile)
Expdataarray=np.array(Expdata)

# Initialize T_FD_all and T_analy_all
T_FD_all=np.empty((Ncellr*Ncellz,tsample.shape[0]))*np.nan
T_analy_all=np.empty((Ncellr*Ncellz,tsample.shape[0]))*np.nan

# Finite Difference Solution
tic=time.time()
#FD_all_data=finiteDifference.uniformgen(GridMapFD,kr,kz,hR,hz0,hzH,rho,cp,R,H,Ncellr,Ncellz,delr,delz,delt,Tinit,Tamb,tsample,qgen,A_antoine,B_antoine,C_antoine,Mgas,P0, GM, tol)
FD_all_data=finiteDifference.Arrheniusgen(GridMapFD,kr,kz,hR,hz0,hzH,rho,cp,R,H,Ncellr,Ncellz,delr,delz,delt,Tinit,Tamb,tsample,qgen,Q0,Ea,Picard,Newton,Nmax,tol,A_antoine,B_antoine,C_antoine,Mgas,P0,GM)
time_FD=time.time()-tic
T_FD_all=FD_all_data[0]
Tavg_FD_all=FD_all_data[1]

# Analytical Solution
tic=time.time()
Analy_data=analytical_nogen.calc_whole(H,R,zLoc,rLoc,tsample,rho,cp,hzH,hR,kz,kr,Tinit,Tamb,GridMap,Ncellr,Ncellz)
time_Analy=time.time()-tic
T_analy_all=Analy_data[0]
Tavg_analy_all=Analy_data[1]
Analy_err_all=Analy_data[2]
print(T_FD_all[:,100])

# writing finite dif temp to csv file 
writing_FD = open('%s_FD_all.csv'%Name, 'w')
with writing_FD:
    writer = csv.writer(writing_FD)
    writer.writerows(T_FD_all)

# writing finite dif temp to csv file 
writing_analy = open('%s_analy_all.csv'%Name, 'w')
with writing_analy:
    writer = csv.writer(writing_analy)
    writer.writerows(T_analy_all)

# writing FDavg temp to csv file 
#writing_avgFD = open('%s_FDavg.csv'%Name, 'w')
#with writing_avgFD:
#    writer = csv.writer(writing_avgFD)
#    writer.writerows(Tavg_FD_all)
np.savetxt('%s_FDavg.csv'%Name,Tavg_FD_all)

# writing FDavg temp to csv file 
#writing_avgAnaly = open('%s_Analyavg.csv'%Name, 'w')
#with writing_avgAnaly:
#    writer = csv.writer(writing_avgAnaly)
#    writer.writerows(Tavg_analy_all)
np.savetxt('%s_Analyavg.csv'%Name,Tavg_analy_all)


parameters = [Ncellr, Ncellz, delt, tfinal, time_FD, time_Analy]
#print(parameters)

# writing 'Ncellr, Ncellz, delt, tfinal' to csv file 
np.savetxt('%s_parameters.csv'%Name, parameters)

print("files written to csv")

tplot=2000
# Obtain contour subplots 
#contourplot.compareanaly(T_FD_all,T_analy_all,rLoc,zLoc,tsample,delt,Ncellr,Ncellz,tplot)

# Obtain lineplots
# Plotting along r at certain z
zplot=H/2 # Plot location
lineplot.alongr_plotT(T_FD_all,T_analy_all,rLoc,zLoc,tsample,R,delz,delt,Ncellr,Ncellz,zplot,tplot)
# Plotting along z at certain r
rplot=R/2 # Plot location
lineplot.alongz_plotT(T_FD_all,T_analy_all,rLoc,zLoc,tsample,H,delr,delt,Ncellr,Ncellz,rplot,tplot)
# Plotting over t on r=R and at certain z
lineplot.surfaceTovertime_onR(T_FD_all,T_analy_all,GridMap,GridMapFD,zLoc,tsample,R,delz,Ncellr,zplot)
# Plotting vent gas amount over time
#lineplot.mG_overtime(mG_FD_all,tsample)
# Plotting Pvapor over time
#lineplot.Pvap_overtime(Pvap_FD_all,tsample)
# Plotting Tavg over time
lineplot.Tavg_overtime(Tavg_FD_all,Tavg_analy_all,Expdataarray,tsample)
# Plotting analytical truncation errors
lineplot.Analy_err_overtime(Analy_err_all,tsample)
# Plotting VolL over time
#lineplot.VL_overtime(VolL_FD_all,tsample)
# Create animation 
#createanimation.compareanaly(T_FD_all,T_analy_all,rLoc,zLoc,tsample,Ncellr,Ncellz)
