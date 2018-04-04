#This script creates animation

import numpy as np
import matplotlib
from matplotlib import pyplot
from matplotlib import animation
import types

def compareanaly(T_FD_all,T_analy_all,rLoc,zLoc,tsample,Ncellr,Ncellz):
    #T_FD=np.rot90(T_FD,3)
    #T_analy=np.rot90(T_analy,3)
    def setvisible(self,vis):
        for c in self.collections: c.set_visible(vis)
    def setanimated(self,ani):
        for c in self.collections: c.set_animated(ani)
    #First set up the figure, axis and plot element we want to animate

    fig, (ax1,ax2)=pyplot.subplots(figsize=(12,4),ncols=2)
    ax1.set_xlabel("r (m)")
    ax1.set_ylabel("z (m)")
    
    ax2.set_xlabel("r (m)")
    ax2.set_ylabel("z (m)")

    #cbar=fig.colorbar(FD,ax=ax1)
    #cbar.ax.set_ylabel('Temperature (K)')
    #cbar=fig.colorbar([],ax=ax1)
    #cbar.ax.set_ylabel('Temperature (K)')
    #cbar2=fig.colorbar([],ax=ax2)
    #cbar2.ax.set_ylabel('Temperature (K)')
    
    FD=ax1.contourf(rLoc,zLoc,np.empty((Ncellr,Ncellz)),np.arange(390,401,0.2))
    Analy=ax2.contourf(rLoc,zLoc,np.empty((Ncellr,Ncellz)),np.arange(390,401,0.2))

    def init():
        FD=ax1.contourf(rLoc,zLoc,np.empty((Ncellr,Ncellz)),np.arange(390,401,0.2))
        Analy=ax2.contourf(rLoc,zLoc,np.empty((Ncellr,Ncellz)),np.arange(390,401,0.2))
        cbar=fig.colorbar(FD,ax=ax1)
        cbar.ax.set_ylabel('Temperature (K)')
        cbar2=fig.colorbar(Analy,ax=ax2)
        cbar2.ax.set_ylabel('Temperature (K)')
        #FD.set_visible = types.MethodType(setvisible,FD)
        FD.set_animated = types.MethodType(setanimated,FD)
        #Analy.set_visible = types.MethodType(setvisible,Analy)
        Analy.set_animated = types.MethodType(setanimated,Analy)
        FD.axes = pyplot.gca()
        Analy.axes = pyplot.gca()
        FD.figure = FD
        Analy.figure = Analy
        #FD.set_data(rLoc,zLoc,np.empty((Ncellr,Ncellz)))
        #cbar.set_data([])
        #Analy.set_data(rLoc,zLoc,np.empty((Ncellr,Ncellz)))
        #cbar2.set_data([])
        ax1.set_title("Finite Difference at t=0 s")
        ax2.set_title("Analytical at t=0 s")
        return (FD,Analy)

    def animate(i):
        T_FD=T_FD_all[:,i].reshape((Ncellr,Ncellz))
        T_analy=T_analy_all[:,i].reshape((Ncellr,Ncellz))
        t=tsample[i]
        FD=ax1.contourf(rLoc,zLoc,T_FD,np.arange(390,401,0.2))
        Analy=ax2.contourf(rLoc,zLoc,T_analy,np.arange(390,401,0.2))
        
        FD.axes = pyplot.gca()
        Analy.axes = pyplot.gca()
        FD.figure = FD
        Analy.figure = Analy
        #FD.set_data(rLoc,zLoc,T_FD)
        #cbar.set_data(FD)
        #Analy.set_data(rLoc,zLoc,T_analy)
        #cbar2.set_data(Analy)
        ax1.set_title("Finite Difference at t=%4.2f s" %(t))
        ax2.set_title("Analytical at t=%4.2f s" %(t))
        #FD.set_visible = types.MethodType(setvisible,FD)
        FD.set_animated = types.MethodType(setanimated,FD)
        #Analy.set_visible = types.MethodType(setvisible,Analy)
        Analy.set_animated = types.MethodType(setanimated,Analy)
        
        return (FD,Analy)
    #FD=ax1.contourf(rLoc,zLoc,T_FD,np.arange(390,401,0.2))
    #cbar=fig.colorbar(FD,ax=ax1)

        ## Bug fix for Quad Contour set not having the attributes 'set_visible' and 'set_animated'
        # ** uncomment the following 2 lines to make the code work:**
        #FD.set_visible = types.MethodType(setvisible,im)
        #FD.set_animated = types.MethodType(setanimated,FD)
        #Analy.set_animated = types.MethodType(setanimated,Analy)
    #FD.axes = plt.gca()
    #FD.figure = fig
    
    #Analy=ax2.contourf(rLoc,zLoc,T_analy,np.arange(390,401,0.2))
    #cbar2=fig.colorbar(Analy,ax=ax2)

    anim=animation.FuncAnimation(fig,animate,init_func=init,frames=tsample.shape[0],interval=1000)
    anim.save('compare_contours.mp4',fps=30,extra_args=['-vcodec', 'libx264'])
    #pyplot.show()