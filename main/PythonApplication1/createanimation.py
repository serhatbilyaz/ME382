#This script creates animation

import numpy as np
import matplotlib
from matplotlib import pyplot
from matplotlib import animation


def compareanaly(T_FD_all,T_analy_all,rLoc,zLoc,tsample,Ncellr,Ncellz):
    print("Creating contourplot animation ...")
    #T_FD=np.rot90(T_FD,3)
    #T_analy=np.rot90(T_analy,3)

    #First set up the figure, axis and plot element we want to animate
    fig, (ax1,ax2)=pyplot.subplots(figsize=(12,4),ncols=2)
    ax1.set_xlabel("r (m)")
    ax1.set_ylabel("z (m)")
    ax2.set_xlabel("r (m)")
    ax2.set_ylabel("z (m)")
    
    FD=ax1.contourf(rLoc,zLoc,np.empty((Ncellz,Ncellr)),np.arange(25,200,1))
    Analy=ax2.contourf(rLoc,zLoc,np.empty((Ncellz,Ncellr)),np.arange(25,200,1))

    def init():
        FD=ax1.contourf(rLoc,zLoc,np.empty((Ncellz,Ncellr)),np.arange(25,200,1))
        Analy=ax2.contourf(rLoc,zLoc,np.empty((Ncellz,Ncellr)),np.arange(25,200,1))
        cbar=fig.colorbar(FD,ax=ax1)
        cbar.ax.set_ylabel('Temperature (C)')
        cbar2=fig.colorbar(Analy,ax=ax2)
        cbar2.ax.set_ylabel('Temperature (C)')
        ax1.set_title("Finite Difference at t=0 s")
        ax2.set_title("Analytical at t=0 s")
        return (FD,Analy)

    def animate(i):
        T_FD=T_FD_all[:,i*100]
        T_FD=T_FD.reshape((Ncellz,Ncellr))
        T_analy=T_analy_all[:,i*100]
        T_analy=T_analy.reshape((Ncellz,Ncellr))
        
        t=tsample[i*100]
        print("t=%2.2f"%t)
        FD=ax1.contourf(rLoc,zLoc,T_FD-273,np.arange(25,160,1))
        Analy=ax2.contourf(rLoc,zLoc,T_analy-273,np.arange(25,160,1))
        ax1.set_title("Finite Difference at t=%4.2f s" %(t))
        ax2.set_title("Analytical at t=%4.2f s" %(t))

        
        return (FD,Analy)

    anim=animation.FuncAnimation(fig,animate,init_func=init,frames=int((tsample.shape[0]-1)/100),interval=1000,repeat_delay=1000)
    anim.save('compare_contours.mp4',fps=8,extra_args=['-vcodec', 'libx264'])
    #pyplot.show()