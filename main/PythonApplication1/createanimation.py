#This script creates animation

import numpy as np
import matplotlib
from matplotlib import pyplot
from matplotlib import animation


def compareanaly(T_FD_all,T_analy_all,rLoc,zLoc,tsample,Ncellr,Ncellz):
    #T_FD=np.rot90(T_FD,3)
    #T_analy=np.rot90(T_analy,3)

    #First set up the figure, axis and plot element we want to animate
    fig, (ax1,ax2)=pyplot.subplots(figsize=(12,4),ncols=2)
    ax1.set_xlabel("r (m)")
    ax1.set_ylabel("z (m)")
    ax2.set_xlabel("r (m)")
    ax2.set_ylabel("z (m)")
    
    FD=ax1.contourf(rLoc,zLoc,np.empty((Ncellr,Ncellz)),np.arange(390,401,0.2))
    Analy=ax2.contourf(rLoc,zLoc,np.empty((Ncellr,Ncellz)),np.arange(390,401,0.2))

    def init():
        FD=ax1.contourf(rLoc,zLoc,np.empty((Ncellr,Ncellz)),np.arange(390,401,0.2))
        Analy=ax2.contourf(rLoc,zLoc,np.empty((Ncellr,Ncellz)),np.arange(390,401,0.2))
        cbar=fig.colorbar(FD,ax=ax1)
        cbar.ax.set_ylabel('Temperature (K)')
        cbar2=fig.colorbar(Analy,ax=ax2)
        cbar2.ax.set_ylabel('Temperature (K)')
        ax1.set_title("Finite Difference at t=0 s")
        ax2.set_title("Analytical at t=0 s")
        return (FD,Analy)

    def animate(i):
        T_FD=T_FD_all[:,i].reshape((Ncellr,Ncellz))
        T_analy=T_analy_all[:,i].reshape((Ncellr,Ncellz))
        t=tsample[i]
        FD=ax1.contourf(rLoc,zLoc,T_FD,np.arange(390,401,0.2))
        Analy=ax2.contourf(rLoc,zLoc,T_analy,np.arange(390,401,0.2))
        ax1.set_title("Finite Difference at t=%4.2f s" %(t))
        ax2.set_title("Analytical at t=%4.2f s" %(t))

        
        return (FD,Analy)

    anim=animation.FuncAnimation(fig,animate,init_func=init,frames=tsample.shape[0],interval=1000,repeat_delay=1000)
    anim.save('compare_contours.mp4',fps=8,extra_args=['-vcodec', 'libx264'])
    #pyplot.show()