# This script plots the radial distribution of temperature for a given z

import numpy as np
import matplotlib
from matplotlib import pyplot

def alongr_plotT(Temp_FD,Temp_analy,r,z,t,R): 
    pyplot.figure()
    pyplot.plot(r,Temp_FD,'--',linewidth=3, label="Finite Difference")
    pyplot.plot(r,Temp_analy,linewidth=3, label="Analytical")

    bounds=np.array([0, R])
    pyplot.xlim(bounds) 
    pyplot.xlabel('r (m)', fontsize=10)
    pyplot.ylabel('T (K)', fontsize=10)
    pyplot.title('Temperature distribution at z=%2.2f m at t=%4.2f s' %(z,t))
    pyplot.legend(loc='best')
    pyplot.savefig('T_alongr_z%2.2f_t%4.1f.pdf'%(z,t), bbox_inches='tight')
    return

def alongz_plotT(Temp_FD,Temp_analy,r,z,t,H): 
    pyplot.figure()
    pyplot.plot(z,Temp_FD,'--',linewidth=3, label="Finite Difference")
    pyplot.plot(z,Temp_analy,linewidth=3, label="Analytical")

    bounds=np.array([0, H])
    pyplot.xlim(bounds) 
    pyplot.xlabel('z (m)', fontsize=10)
    pyplot.ylabel('T (K)', fontsize=10)
    pyplot.title('Temperature distribution at r=%2.2f m at t=%4.2f s' %(r,t))
    pyplot.legend(loc='best')
    pyplot.savefig('T_alongz_r%2.2f_t%4.1f.pdf'%(r,t), bbox_inches='tight')
    return