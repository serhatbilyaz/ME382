# This script contains lineplot subroutines

import numpy as np
import matplotlib
from matplotlib import pyplot

def alongr_plotT(T_FD_all,T_analy_all,rLoc,zLoc,tsample,R,delz,delt,Ncellr,Ncellz,z,t):
    i=t/delt
    i=int(np.around(i)) # Time index
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
    pyplot.plot(rplot_alongr,Temp_FD_alongr,'--',linewidth=3, label="Finite Difference")
    pyplot.plot(rplot_alongr,Temp_analy_alongr,linewidth=3, label="Analytical")

    bounds=np.array([0, R])
    pyplot.xlim(bounds) 
    pyplot.xlabel('r (m)', fontsize=10)
    pyplot.ylabel('T (K)', fontsize=10)
    pyplot.title('Temperature distribution at z=%2.2f m at t=%4.2f s' %(zplot_alongr,t))
    pyplot.legend(loc='best')
    pyplot.savefig('T_alongr_z%2.2f_t%4.1f.pdf'%(zplot_alongr,t), bbox_inches='tight')
    return

def alongz_plotT(T_FD_all,T_analy_all,rLoc,zLoc,tsample,H,delr,delt,Ncellr,Ncellz,r,t):
    i=t/delt
    i=int(np.around(i)) # Time index
    t=tsample[i]
    print(t)
    print(i)
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
    pyplot.ylabel('T (K)', fontsize=10)
    pyplot.title('Temperature distribution at r=%2.2f m at t=%4.2f s' %(rplot_alongz,t))
    pyplot.legend(loc='best')
    pyplot.savefig('T_alongz_r%2.2f_t%4.1f.pdf'%(rplot_alongz,t), bbox_inches='tight')
    return