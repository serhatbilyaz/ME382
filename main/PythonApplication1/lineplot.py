# This script contains lineplot subroutines

import numpy as np
import matplotlib
from matplotlib import pyplot

def alongr_plotT(T_FD_all,T_analy_all,rLoc,zLoc,tsample,R,delz,delt,Ncellr,Ncellz,z,t):
    print("Creating lineplot along r ...")
    i=t/delt
    i=int(np.around(i)) # Time index
    if i==tsample.shape[0]:
        i=i-1
    t=tsample[i]
    
    m=(z/delz)-0.5 #Temp index in terms of z
    m=int(np.floor(m)) # Find the closest integer to be an index
    rplot_alongr=rLoc[m,:] #r locations for plotting
    zplot_alongr=zLoc[m,0] #Evaluate real z (According to rounded index)
    #print(rplot_alongr[3])
    Temp_FD_alongr=T_FD_all[:,i]
    Temp_FD_alongr=Temp_FD_alongr.reshape((Ncellz,Ncellr))[m,:]
    Temp_analy_alongr=T_analy_all[:,i]
    Temp_analy_alongr=Temp_analy_alongr.reshape((Ncellz,Ncellr))[m,:]
    pyplot.figure()
    pyplot.plot(rplot_alongr,Temp_FD_alongr-273,'--',linewidth=3, label="Finite Difference")
    pyplot.plot(rplot_alongr,Temp_analy_alongr-273,linewidth=3, label="Analytical")

    bounds=np.array([0, R])
    pyplot.xlim(bounds) 
    pyplot.xlabel('r (m)', fontsize=10)
    pyplot.ylabel('T (C)', fontsize=10)
    pyplot.title('Temperature distribution at z=%2.4f m at t=%4.2f s' %(zplot_alongr,t))
    pyplot.legend(loc='best')
    pyplot.savefig('T_alongr_z%2.4f_t%4.1f.pdf'%(zplot_alongr,t), bbox_inches='tight')
    return

def alongz_plotT(T_FD_all,T_analy_all,rLoc,zLoc,tsample,H,delr,delt,Ncellr,Ncellz,r,t):
    print("Creating lineplot along z ...")
    i=t/delt
    i=int(np.around(i)) # Time index
    if i==tsample.shape[0]:
        i=i-1
    t=tsample[i]
    #print(t)
    #print(i)
    n=(r/delr)-0.5 #Temp index in terms of r
    n=int(np.floor(n)) # Find the closest integer to be an index
    rplot_alongz=rLoc[0,n] # Evaluate real r (According to rounded index)
    zplot_alongz=zLoc[:,n].reshape(Ncellz,) # z locations for plotting
    
    Temp_FD_alongz=T_FD_all[:,i]
    Temp_FD_alongz=Temp_FD_alongz.reshape((Ncellz,Ncellr))[:,n]
    Temp_analy_alongz=T_analy_all[:,i]
    Temp_analy_alongz=Temp_analy_alongz.reshape((Ncellz,Ncellr))[:,n]
    pyplot.figure()
    pyplot.plot(zplot_alongz,Temp_FD_alongz,'--',linewidth=3, label="Finite Difference")
    pyplot.plot(zplot_alongz,Temp_analy_alongz,linewidth=3, label="Analytical")

    bounds=np.array([0, H])
    pyplot.xlim(bounds) 
    pyplot.xlabel('z (m)', fontsize=10)
    pyplot.ylabel('T (C)', fontsize=10)
    pyplot.title('Temperature distribution at r=%2.4f m at t=%4.2f s' %(rplot_alongz,t))
    pyplot.legend(loc='best')
    pyplot.savefig('T_alongz_r%2.4f_t%4.1f.pdf'%(rplot_alongz,t), bbox_inches='tight')
    return

def surfaceTovertime_onR(T_FD_all,T_analy_all,GridMap,GridMapFD,zLoc,tsample,R,delz,Ncellr,z):
    print("Creating lineplot of surface temperature at r=R over time ...")
    m=(z/delz)-0.5 #Temp index in terms of z
    m=int(np.floor(m)) # Find the closest integer to be an index
    zplot=zLoc[m,0] #Evaluate real z (According to rounded index)

    Temp_FD_overt=T_FD_all[GridMapFD[Ncellr-1,m],:]
    Temp_analy_overt=T_analy_all[GridMap[m,Ncellr-1],:]

    pyplot.figure()
    pyplot.plot(tsample/60,Temp_FD_overt-273,'--',linewidth=3, label="Finite Difference")
    pyplot.plot(tsample/60,Temp_analy_overt-273,linewidth=3, label="Analytical")

    bounds=np.array([0, 200])
    pyplot.ylim(bounds) 
    pyplot.xlabel('t (min)', fontsize=10)
    pyplot.ylabel('T (C)', fontsize=10)
    pyplot.title('Temperature at r=R and z=%2.2f s' %(zplot))
    pyplot.xticks(np.arange(0, tsample[tsample.shape[0]-1]/60, step=5))
    pyplot.legend(loc='best')
    pyplot.savefig('Tovert_onR_z%2.2f.pdf'%(zplot), bbox_inches='tight')
    return

def mG_overtime(mG_FD_all,tsample):
    print("Creating lineplot of vent gas amount over time ...")
    
    pyplot.figure()
    pyplot.plot(tsample,mG_FD_all,'-',linewidth=3, label="Finite Difference")
#    pyplot.plot(tsample,Temp_analy_overt,linewidth=3, label="Analytical")

    #bounds=np.array([0, H])
    #pyplot.xlim(bounds) 
    pyplot.xlabel('t (s)', fontsize=10)
    pyplot.ylabel('m_G (kg)', fontsize=10)
    pyplot.title('Accumulation of vent gas')
    pyplot.legend(loc='best')
    pyplot.savefig('mG_overt.pdf', bbox_inches='tight')
    return

def VL_overtime(VL_FD_all,tsample):
    print("Creating lineplot of liquid volume over time ...")
    
    pyplot.figure()
    pyplot.plot(tsample,VL_FD_all,'-',linewidth=3, label="Finite Difference")
#    pyplot.plot(tsample,Temp_analy_overt,linewidth=3, label="Analytical")

    #bounds=np.array([0, H])
    #pyplot.xlim(bounds) 
    pyplot.xlabel('t (s)', fontsize=10)
    pyplot.ylabel('V_L (m^3)', fontsize=10)
    pyplot.title('Change of liquid volume')
    pyplot.legend(loc='best')
    pyplot.savefig('VL_overt.pdf', bbox_inches='tight')
    return

def Pvap_overtime(Pvap_FD_all,tsample):
    print("Creating lineplot of Pvap over time ...")
    
    pyplot.figure()
    pyplot.plot(tsample,Pvap_FD_all,'-',linewidth=3, label="Finite Difference")
#    pyplot.plot(tsample,Temp_analy_overt,linewidth=3, label="Analytical")

    #bounds=np.array([0, H])
    #pyplot.xlim(bounds) 
    pyplot.xlabel('t (s)', fontsize=10)
    pyplot.ylabel('Pvap (Pa)', fontsize=10)
    pyplot.title('Change of vapor pressure')
    pyplot.legend(loc='best')
    pyplot.savefig('Pvap_overt.pdf', bbox_inches='tight')
    return

def Tavg_overtime(Tavg_FD_all,Tavg_analy_all,tsample):
    print("Creating lineplot of Tavg over time ...")
    
    pyplot.figure()
    pyplot.plot(tsample/60,Tavg_FD_all-273,'--',linewidth=3, label="Finite Difference")
    pyplot.plot(tsample/60,Tavg_analy_all-273,linewidth=3, label="Analytical (No Qgen)")

    bounds=np.array([0, 200])
    pyplot.ylim(bounds) 
    pyplot.xlabel('t (min)', fontsize=10)
    pyplot.ylabel('Temperature (C)', fontsize=10)
    pyplot.title('Change of average temperature')
    pyplot.legend(loc='best')
    pyplot.xticks(np.arange(0, tsample[tsample.shape[0]-1]/60, step=5))
    pyplot.savefig('Tavg_overt.pdf', bbox_inches='tight')
    return