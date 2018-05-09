# This script reads the solution of different cases and compare them

import numpy as np
import matplotlib
import pandas as pd
from matplotlib import pyplot


NoGen=0
ArrGen=1
GridIndependence=1
LinearSolver=0
CompareExp=0
CompareAnaly=0

if ArrGen:
    if GridIndependence:
        #AG5_param_data=pd.read_csv('AG5/AG5_parameters.csv',header=None)
        #AG5_param=np.array(AG5_param_data)

        #AG6_param_data=pd.read_csv('AG6/AG6_parameters.csv',header=None)
        #AG6_param=np.array(AG6_param_data)

        AG7_param_data=pd.read_csv('AG7/AG7_parameters.csv',header=None)
        AG7_param=np.array(AG7_param_data)

        AG8_param_data=pd.read_csv('AG8/AG8_parameters.csv',header=None)
        AG8_param=np.array(AG8_param_data)

        AG9_param_data=pd.read_csv('AG9/AG9_parameters.csv',header=None)
        AG9_param=np.array(AG9_param_data)

        #deltAG5=AG5_param[2]
        #tfinalAG5=AG5_param[3]
        #print(deltAG5)
        #print(tfinalAG5)
        #tsampleAG5=np.linspace(0,tfinalAG5,tfinalAG5/deltAG5+1)

        #deltAG6=AG6_param[2]
        #tfinalAG6=AG6_param[3]
        #print(deltAG6)
        #print(tfinalAG6)
        #tsampleAG6=np.linspace(0,tfinalAG6,tfinalAG6/deltAG6+1)

        deltAG7=AG7_param[2]
        tfinalAG7=AG7_param[3]
        print(deltAG7)
        print(tfinalAG7)
        tsampleAG7=np.linspace(0,tfinalAG7,tfinalAG7/deltAG7+1)

        deltAG8=AG8_param[2]
        tfinalAG8=AG8_param[3]
        print(deltAG8)
        print(tfinalAG8)
        tsampleAG8=np.linspace(0,tfinalAG8,tfinalAG8/deltAG8+1)

        deltAG9=AG9_param[2]
        tfinalAG9=AG9_param[3]
        print(deltAG9)
        print(tfinalAG9)
        tsampleAG9=np.linspace(0,tfinalAG9,tfinalAG9/deltAG9+1)

        #AG5_FDavg_data=pd.read_csv('AG5/AG5_FDavg.csv',header=None)
        #AG5_FDavg=np.array(AG5_FDavg_data)

        #AG6_FDavg_data=pd.read_csv('AG6/AG6_FDavg.csv',header=None)
        #AG6_FDavg=np.array(AG6_FDavg_data)

        AG7_FDavg_data=pd.read_csv('AG7/AG7_FDavg.csv',header=None)
        AG7_FDavg=np.array(AG7_FDavg_data)

        AG8_FDavg_data=pd.read_csv('AG8/AG8_FDavg.csv',header=None)
        AG8_FDavg=np.array(AG8_FDavg_data)

        AG9_FDavg_data=pd.read_csv('AG9/AG9_FDavg.csv',header=None)
        AG9_FDavg=np.array(AG9_FDavg_data)
        print(np.shape(tsampleAG7))
        print(np.shape(tsampleAG9))
        print(np.shape(AG9_FDavg))
        TinterpAG9=np.interp(tsampleAG7,tsampleAG9,np.reshape(AG9_FDavg,np.shape(tsampleAG9)))
        TinterpAG8=np.interp(tsampleAG7,tsampleAG8,np.reshape(AG8_FDavg,np.shape(tsampleAG8)))
        #print(np.shape(TinterpAG9))
        #print(np.shape(AG7_FDavg))
        
        #TinterpAG7=np.interp(tsampleAG5,tsampleAG7,np.reshape(AG7_FDavg,np.shape(tsampleAG7)))
        #TinterpAG6=np.interp(tsampleAG5,tsampleAG6,np.reshape(AG6_FDavg,np.shape(tsampleAG6)))

        Err9_7=np.abs(TinterpAG9-np.reshape(AG7_FDavg,np.shape(tsampleAG7)))*100/np.reshape(AG7_FDavg,np.shape(tsampleAG7))
        Err8_7=np.abs(TinterpAG8-np.reshape(AG7_FDavg,np.shape(tsampleAG7)))*100/np.reshape(AG7_FDavg,np.shape(tsampleAG7))

        print("Creating lineplot of Tavg over time for ArrGen GridIndependence...")
        pyplot.figure()
        #pyplot.plot(tsampleAG5/60,AG5_FDavg-273,'-',linewidth=3, label="N=60 dt=0.2")
        #pyplot.plot(tsampleAG6/60,AG6_FDavg-273,'--',linewidth=3, label="N=50 dt=0.3")
        pyplot.plot(tsampleAG7/60,AG7_FDavg-273,'-.',linewidth=3, label="N=40 dt=0.5")
        pyplot.plot(tsampleAG8/60,AG8_FDavg-273,':',linewidth=3, label="N=30 dt=0.75")
        pyplot.plot(tsampleAG9/60,AG9_FDavg-273,'-',linewidth=3, label="N=20 dt=1")
#        pyplot.plot(Expdataarray[:,0],Expdataarray[:,1],'-.',linewidth=3,label="Experimental")
        bounds=np.array([0, 200])
        pyplot.ylim(bounds) 
        pyplot.xlabel('t (min)', fontsize=10)
        pyplot.ylabel('Temperature (C)', fontsize=10)
        pyplot.title('Change of average temperature')
        pyplot.legend(loc='best')
        pyplot.xticks(np.arange(0, tsampleAG9[tsampleAG9.shape[0]-1]/60, step=5))
        pyplot.savefig('ArrGen_GridIndep_Tavg.pdf', bbox_inches='tight')

        print("Creating lineplot of Tavg over time for ArrGen GridIndependence...")
        pyplot.figure()
        pyplot.plot(tsampleAG7/60,Err8_7,'--',linewidth=3, label="Err8_7")
        pyplot.plot(tsampleAG7/60,Err9_7,':',linewidth=3, label="Err9_7")

#        pyplot.plot(Expdataarray[:,0],Expdataarray[:,1],'-.',linewidth=3,label="Experimental")
        #bounds=np.array([0, 10])
        #pyplot.ylim(bounds) 
        pyplot.xlabel('t (min)', fontsize=10)
        pyplot.ylabel('% Error', fontsize=10)
        pyplot.title('Percent error in average temperature')
        pyplot.legend(loc='best')
        pyplot.xticks(np.arange(0, tsampleAG9[tsampleAG9.shape[0]-1]/60, step=5))
        pyplot.savefig('ArrGen_GridIndepErr_Tavg.pdf', bbox_inches='tight')



    if LinearSolver:
        # Read parameters
        AG2_param_data=pd.read_csv('AG2/AG2_parameters.csv',header=None)
        AG2_param=np.array(AG2_param_data)
        
        delt=AG2_param[2]
        tfinal=AG2_param[3]
        print(delt)
        print(tfinal)
        tsample=np.linspace(0,tfinal,tfinal/delt+1)

        AG1_FDavg_data=pd.read_csv('AG1/AG1_FDavg.csv',header=None)
        AG1_FDavg=np.array(AG1_FDavg_data)

        AG2_FDavg_data=pd.read_csv('AG2/AG2_FDavg.csv',header=None)
        AG2_FDavg=np.array(AG2_FDavg_data)

        AG3_FDavg_data=pd.read_csv('AG3/AG3_FDavg.csv',header=None)
        AG3_FDavg=np.array(AG3_FDavg_data)

        AG4_FDavg_data=pd.read_csv('AG4/AG4_FDavg.csv',header=None)
        AG4_FDavg=np.array(AG4_FDavg_data)
	#print(AG4_FDavg)
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