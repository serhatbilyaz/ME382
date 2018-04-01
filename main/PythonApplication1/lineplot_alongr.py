# This script plots the radial distribution of temperature for a given z

import numpy as np
import matplotlib
from matplotlib import pyplot

def plotT(T,z,R): 
    #z=?
    #Temp=? (from T)
    pyplot.figure()
    pyplot.plot(r,Temp,'--',linewidth=3, label="z=")
    #pyplot.plot(ExpX_Uy,ExpUy,linewidth=3, label="Exp")

    bounds=np.array([0, R])
    pyplot.xlim(bounds) 
    pyplot.xlabel('r (m)', fontsize=10)
    pyplot.ylabel('$T$ (K)', fontsize=10)
    pyplot.legend(loc='upper left')
    pyplot.savefig('T_alongr.pdf', bbox_inches='tight')
    return 