# This script reads the solution of different cases and compare them

import numpy as np
import matplotlib
import pandas as pd
from matplotlib import pyplot


NoGen=0
ArrGen=0
GridIndependence=0
LinearSolver=0
CompareExp=0
CompareAnaly=0

if ArrGen:
#    if GridIndependence:
    if LinearSolver:
        # Read parameters
        AG2_param_data=pd.read_csv('AG2/AG2_parameters.csv',header=None)
        AG2_param=np.array(AG2_param_data)
        
        delt=AG2_param[2]
        tfinal=AG2_param[3]
        tsample=np.linspace(0,tfinal,tfinal/delt+1)
        AG1_FDavg_data=pd.read_csv('AG1/AG1_FDavg.csv',header=None)
        AG1_FDavg=np.array(AG1_FDavg_data)

        AG2_FDavg_data=pd.read_csv('AG2/AG2_FDavg.csv',header=None)
        AG2_FDavg=np.array(AG2_FDavg_data)

        AG3_FDavg_data=pd.read_csv('AG3/AG3_FDavg.csv',header=None)
        AG3_FDavg=np.array(AG3_FDavg_data)

        AG4_FDavg_data=pd.read_csv('AG4/AG4_FDavg.csv',header=None)
        AG4_FDavg=np.array(AG4_FDavg_data)

        print("Creating lineplot of Tavg over time for ArrGen LinearSolver...")
        pyplot.figure()
        pyplot.plot(tsample/60,AG1_FDavg-273,'-',linewidth=3, label="LU solver")
        pyplot.plot(tsample/60,AG2_FDavg-273,'--',linewidth=3, label="GMRES tol 1e-2")
        pyplot.plot(tsample/60,AG3_FDavg-273,'-.',linewidth=3, label="GMRES tol 1e-3")
        pyplot.plot(tsample/60,AG4_FDavg-273,'-',linewidth=3, label="GMRES tol 1e-4")

#        pyplot.plot(Expdataarray[:,0],Expdataarray[:,1],'-.',linewidth=3,label="Experimental")

        bounds=np.array([0, 200])
        pyplot.ylim(bounds) 
        pyplot.xlabel('t (min)', fontsize=10)
        pyplot.ylabel('Temperature (C)', fontsize=10)
        pyplot.title('Change of average temperature')
        pyplot.legend(loc='best')
        pyplot.xticks(np.arange(0, tsample[tsample.shape[0]-1]/60, step=5))
        pyplot.savefig('ArrGen_LinSol_Tavg.pdf', bbox_inches='tight')
        #AG2_T_FD_all_data=pd.read_csv('AG2_FD_all.csv')
        #AG2_T_FD_all=np.array(AG2_T_FD_all_data)

 #   if CompareExp:


#print(AG2_FDavg)

#if NoGen:
#    if GridIndependence:

#    if CompareAnaly: