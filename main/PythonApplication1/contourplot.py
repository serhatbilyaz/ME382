# This script plots the temperature contour plots

import numpy as np
import matplotlib
from matplotlib import pyplot


def compareanaly(T_FD,T_analy,rLoc,zLoc,t):
    T_FD=np.rot90(T_FD,3)
    T_analy=np.rot90(T_analy,3)
    fig, (ax1,ax2)=pyplot.subplots(figsize=(12,4),ncols=2)
    FD=ax1.contourf(rLoc,zLoc,T_FD,np.arange(390,401,0.2))
    cbar=fig.colorbar(FD,ax=ax1)
    cbar.ax.set_ylabel('Temperature (K)')
    ax1.set_xlabel("r (m)")
    ax1.set_ylabel("z (m)")
    ax1.set_title("Finite Difference at t=%4.2f s" %(t))
    Analy=ax2.contourf(rLoc,zLoc,T_analy,np.arange(390,401,0.2))
    cbar2=fig.colorbar(Analy,ax=ax2)
    cbar2.ax.set_ylabel('Temperature (K)')
    ax2.set_xlabel("r (m)")
    ax2.set_ylabel("z (m)")
    ax2.set_title("Analytical at t=%4.2f s" %(t))
    
    pyplot.savefig('Tcontour_compare_t%4.1f.pdf'%(t), bbox_inches='tight')
