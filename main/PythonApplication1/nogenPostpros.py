import numpy as np
import matplotlib
import pandas as pd
from matplotlib import pyplot

H=1
R=1
z=H/2
r=R/2
t=100   #pick time at which error is calcualted 
namepara='AG7_parameters.csv'   #change names of files here
nameFDall='AG7_FD_all.csv'
nameAnalyall='AG7_FD_all.csv'

#reading parameter file
param_data=pd.read_csv(namepara,header=None)    
param=np.array(param_data)

Ncellr=10 # Ncellr=param[0] <- it didnt like when i pull the value out of the array, it gives an error in line 29
Ncellz=10 # Ncellz=param[1]
delt=param[2]
tfinal=param[3]
delr=R/Ncellr
delz=H/Ncellz
print(delt)
print(tfinal)
tsample=np.linspace(0,tfinal,tfinal/delt+1)

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


#reading analytical file
analy_data=pd.read_csv(nameAnalyall,header=None)
T_analy_all=np.array(analy_data)

#reading FD file
FD_data=pd.read_csv(nameFDall,header=None)
T_FD_all=np.array(FD_data)

#finding time, r and z locations based on the lineplot file 
i=t     #something might be wring here.. it didnt like when i had an experssion i=t/delt, i should be the time location of intreset 
print(t)

#finding r location 
m=(z/delz)-0.5 #Temp index in terms of z
m=int(np.floor(m)) # Find the closest integer to be an index
rplot_alongr=rLoc[m,:] #r locations for plotting
zplot_alongr=zLoc[m,0] #Evaluate real z (According to rounded index)

#finding z location
n=(r/delr)-0.5 #Temp index in terms of r
n=int(np.floor(n)) # Find the closest integer to be an index
rplot_alongz=rLoc[0,n] # Evaluate real r (According to rounded index)
zplot_alongz=zLoc[:,n].reshape(Ncellz,) # z locations for plotting

# selecting temp array at time i and finding temp along r
Temp_FD=T_FD_all[:,i]
Temp_FD=Temp_FD.reshape((Ncellz,Ncellr))
Temp_analy=T_analy_all[:,i]
Temp_analy=Temp_analy.reshape((Ncellz,Ncellr))

#error along r
Temp_FD_alongr=Temp_FD[m,:]
Temp_analy_alongr=Temp_analy[m,:]
Err_alongr=np.abs(Temp_analy_alongr-Temp_FD_alongr)*100/Temp_analy_alongr

#error along z
Temp_FD_alongz=Temp_FD[:,n]
Temp_analy_alongz=Temp_analy[:,n]
Err_alongz=np.abs(Temp_analy_alongz-Temp_FD_alongz)*100/Temp_analy_alongz




