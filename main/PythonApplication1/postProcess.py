# This script reads the solution of different cases and compare them

import numpy as np
import matplotlib
import pandas as pd
from matplotlib import pyplot
import lineplot
import contourplot
import createanimation

NoGen=0
ArrGen=0
GridIndependence=0
LinearSolver=0
CompareExp=0
CompareAnaly=0
Tavgcomparison=0
rdist=0
zdist=0
Contourplot=1
Animation=0
Analytruncation=0
printtimes=0
R=0.009 # Radius of the cylinder (m)
H=0.065 # Height of the cylinder (m)






if ArrGen:
    if GridIndependence:
        AG5_param_data=pd.read_csv('AG5/AG5_parameters.csv',header=None)
        AG5_param=np.array(AG5_param_data)

        #AG6_param_data=pd.read_csv('AG6/AG6_parameters.csv',header=None)
        #AG6_param=np.array(AG6_param_data)

        AG7_param_data=pd.read_csv('AG7/AG7_parameters.csv',header=None)
        AG7_param=np.array(AG7_param_data)

        AG1_param_data=pd.read_csv('AG1/AG1_parameters.csv',header=None)
        AG1_param=np.array(AG1_param_data)

        AG9_param_data=pd.read_csv('AG9/AG9_parameters.csv',header=None)
        AG9_param=np.array(AG9_param_data)

        deltAG5=AG5_param[2]
        tfinalAG5=AG5_param[3]
        print(deltAG5)
        print(tfinalAG5)
        tsampleAG5=np.linspace(0,tfinalAG5,tfinalAG5/deltAG5+1)

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

        deltAG1=AG1_param[2]
        tfinalAG1=AG1_param[3]
        print(deltAG1)
        print(tfinalAG1)
        tsampleAG1=np.linspace(0,tfinalAG1,tfinalAG1/deltAG1+1)

        deltAG9=AG9_param[2]
        tfinalAG9=AG9_param[3]
        print(deltAG9)
        print(tfinalAG9)
        tsampleAG9=np.linspace(0,tfinalAG9,tfinalAG9/deltAG9+1)

        AG5_FDavg_data=pd.read_csv('AG5/AG5_FDavg.csv',header=None)
        AG5_FDavg=np.array(AG5_FDavg_data)

        #AG6_FDavg_data=pd.read_csv('AG6/AG6_FDavg.csv',header=None)
        #AG6_FDavg=np.array(AG6_FDavg_data)

        AG7_FDavg_data=pd.read_csv('AG7/AG7_FDavg.csv',header=None)
        AG7_FDavg=np.array(AG7_FDavg_data)

        AG1_FDavg_data=pd.read_csv('AG1/AG1_FDavg.csv',header=None)
        AG1_FDavg=np.array(AG1_FDavg_data)

        AG9_FDavg_data=pd.read_csv('AG9/AG9_FDavg.csv',header=None)
        AG9_FDavg=np.array(AG9_FDavg_data)
        #print(np.shape(tsampleAG7))
        #print(np.shape(tsampleAG9))
        #print(np.shape(AG9_FDavg))
        TinterpAG9=np.interp(tsampleAG5,tsampleAG9,np.reshape(AG9_FDavg,np.shape(tsampleAG9)))
        TinterpAG1=np.interp(tsampleAG5,tsampleAG1,np.reshape(AG1_FDavg,np.shape(tsampleAG1)))
        TinterpAG7=np.interp(tsampleAG5,tsampleAG7,np.reshape(AG7_FDavg,np.shape(tsampleAG7)))
        #print(np.shape(TinterpAG9))
        #print(np.shape(AG7_FDavg))
        
        #TinterpAG7=np.interp(tsampleAG5,tsampleAG7,np.reshape(AG7_FDavg,np.shape(tsampleAG7)))
        #TinterpAG6=np.interp(tsampleAG5,tsampleAG6,np.reshape(AG6_FDavg,np.shape(tsampleAG6)))

        Err9_5=np.abs(TinterpAG9-np.reshape(AG5_FDavg,np.shape(tsampleAG5)))*100/np.reshape(AG5_FDavg,np.shape(tsampleAG5))
        Err1_5=np.abs(TinterpAG1-np.reshape(AG5_FDavg,np.shape(tsampleAG5)))*100/np.reshape(AG5_FDavg,np.shape(tsampleAG5))
        Err7_5=np.abs(TinterpAG7-np.reshape(AG5_FDavg,np.shape(tsampleAG5)))*100/np.reshape(AG5_FDavg,np.shape(tsampleAG5))

        print("Creating lineplot of Tavg over time for ArrGen GridIndependence...")
        pyplot.figure()
        pyplot.plot(tsampleAG5/60,AG5_FDavg-273,'-',linewidth=3, label="N=60 dt=0.2")
        #pyplot.plot(tsampleAG6/60,AG6_FDavg-273,'--',linewidth=3, label="N=50 dt=0.3")
        pyplot.plot(tsampleAG7/60,AG7_FDavg-273,'-.',linewidth=3, label="N=40 dt=0.5")
        pyplot.plot(tsampleAG1/60,AG1_FDavg-273,':',linewidth=3, label="N=30 dt=0.75")
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

        print("Creating Error lineplot of Tavg over time for ArrGen GridIndependence...")
        pyplot.figure()
        pyplot.plot(tsampleAG5/60,Err7_5,':',linewidth=3, label="N=40 dt=0.5")
        pyplot.plot(tsampleAG5/60,Err1_5,'--',linewidth=3, label="N=30 dt=0.75")
        pyplot.plot(tsampleAG5/60,Err9_5,':',linewidth=3, label="N=20 dt=1")

#        pyplot.plot(Expdataarray[:,0],Expdataarray[:,1],'-.',linewidth=3,label="Experimental")
        #bounds=np.array([0, 10])
        #pyplot.ylim(bounds) 
        pyplot.xlabel('t (min)', fontsize=10)
        pyplot.ylabel('% Error', fontsize=10)
        pyplot.title('Percent error in average temperature compared to N=60 dt=0.2')
        pyplot.legend(loc='best')
        pyplot.xticks(np.arange(0, tsampleAG9[tsampleAG9.shape[0]-1]/60, step=5))
        pyplot.savefig('ArrGen_GridIndepErr_Tavg.pdf', bbox_inches='tight')



    if LinearSolver:
        # Read parameters
        AG2_param_data=pd.read_csv('AG2_P_GMRES_true/AG2_parameters.csv',header=None)
        AG2_param=np.array(AG2_param_data)
        
        delt=AG2_param[2]
        tfinal=AG2_param[3]
        print(delt)
        print(tfinal)
        tsample=np.linspace(0,tfinal,tfinal/delt+1)

        AG1_FDavg_data=pd.read_csv('AG1/AG1_FDavg.csv',header=None)
        AG1_FDavg=np.array(AG1_FDavg_data)

        AG2_FDavg_data=pd.read_csv('AG2_P_GMRES_true/AG2_FDavg.csv',header=None)
        AG2_FDavg=np.array(AG2_FDavg_data)

        AG3_FDavg_data=pd.read_csv('AG3_P_GMRES_true/AG3_FDavg.csv',header=None)
        AG3_FDavg=np.array(AG3_FDavg_data)

        AG4_FDavg_data=pd.read_csv('AG4_P_GMRES_true/AG4_FDavg.csv',header=None)
        AG4_FDavg=np.array(AG4_FDavg_data)

        Err2=np.abs(AG2_FDavg-AG1_FDavg)*100/AG1_FDavg
        Err3=np.abs(AG3_FDavg-AG1_FDavg)*100/AG1_FDavg
        Err4=np.abs(AG4_FDavg-AG1_FDavg)*100/AG1_FDavg

	#print(AG4_FDavg)
        print("Creating lineplot of Tavg over time for ArrGen LinearSolver...")
        pyplot.figure()
        pyplot.plot(tsample/60,AG1_FDavg-273,'-',linewidth=3, label="LU solver")
        pyplot.plot(tsample/60,AG2_FDavg-273,'--',linewidth=3, label="GMRES tol 1e-4")
        pyplot.plot(tsample/60,AG3_FDavg-273,'-.',linewidth=3, label="GMRES tol 1e-5")
        pyplot.plot(tsample/60,AG4_FDavg-273,'-',linewidth=3, label="GMRES tol 1e-6")

#        pyplot.plot(Expdataarray[:,0],Expdataarray[:,1],'-.',linewidth=3,label="Experimental")

        bounds=np.array([0, 200])
        pyplot.ylim(bounds) 
        pyplot.xlabel('t (min)', fontsize=10)
        pyplot.ylabel('Temperature (C)', fontsize=10)
        pyplot.title('Change of average temperature')
        pyplot.legend(loc='best')
        pyplot.xticks(np.arange(0, tsample[tsample.shape[0]-1]/60, step=5))
        pyplot.savefig('ArrGen_LinSol_Tavg.pdf', bbox_inches='tight')

        print("Creating lineplot of ErrTavg over time for ArrGen LinearSolver...")
        pyplot.figure()
        pyplot.plot(tsample/60,Err2,'--',linewidth=3, label="GMRES tol 1e-4")
        pyplot.plot(tsample/60,Err3,'-.',linewidth=3, label="GMRES tol 1e-5")
        pyplot.plot(tsample/60,Err4,':',linewidth=3, label="GMRES tol 1e-6")

#        pyplot.plot(Expdataarray[:,0],Expdataarray[:,1],'-.',linewidth=3,label="Experimental")

        #bounds=np.array([0, 200])
        #pyplot.ylim(bounds) 
        pyplot.xlabel('t (min)', fontsize=10)
        pyplot.ylabel('% Error ', fontsize=10)
        pyplot.title('Percent error compared to LU solution')
        pyplot.legend(loc='best')
        pyplot.xticks(np.arange(0, tsample[tsample.shape[0]-1]/60, step=5))
        pyplot.savefig('ArrGen_ErrLinSol_Tavg.pdf', bbox_inches='tight')


    if CompareExp:
        # Read experimental data
        Expdata=pd.read_csv('expdatacsv.csv',header=None)
        Expdataarray=np.array(Expdata)
        ExpTime=Expdataarray[:,0]
        ExpTemp=Expdataarray[:,1]

        AG5_param_data=pd.read_csv('AG5/AG5_parameters.csv',header=None)
        AG5_param=np.array(AG5_param_data)
        NcellrAG5=int(AG5_param[0])
        NcellzAG5=int(AG5_param[1])
        deltAG5=AG5_param[2]
        tfinalAG5=AG5_param[3]
        print(deltAG5)
        print(tfinalAG5)
        print(NcellrAG5)
        print(NcellzAG5)
        tsampleAG5=np.linspace(0,tfinalAG5,tfinalAG5/deltAG5+1)

        # Grid numbering (Row increase => r increase, Column increase => z increase)
        GridMap=np.zeros((NcellzAG5,NcellrAG5),dtype=np.int)
        for i in range(0,NcellzAG5):
            GridMap[i,:]=(i*NcellrAG5)+np.arange(0,NcellrAG5,1)
        GridMapFD=GridMap
        GridMapFD=np.transpose(GridMapFD)

        
        AG5_FDavg_data=pd.read_csv('AG5/AG5_FDavg.csv',header=None)
        AG5_FDavg=np.reshape(np.array(AG5_FDavg_data),np.shape(tsampleAG5))

        #Taking all data and extract surface temperature from all data
        AG5_T_FD_all_data=pd.read_csv('AG5/AG5_FD_all.csv',header=None)
        AG5_T_FD_all=np.array(AG5_T_FD_all_data)
        
        
        m=int(NcellzAG5/2)
        AG5_T_surf=AG5_T_FD_all[GridMapFD[NcellrAG5-1,m],:]
        #Temp_analy_overt=T_analy_all[GridMap[m,Ncellr-1],:]

        TinterpsurfAG5=np.interp(ExpTime*60,tsampleAG5,np.reshape(AG5_FDavg,np.shape(tsampleAG5)))
        TinterpAG5=np.interp(ExpTime*60,tsampleAG5,AG5_FDavg)
        ErrExp_avg=np.abs(TinterpAG5-273.15-ExpTemp)*100/(ExpTemp+273)
        ErrExp_surf=np.abs(TinterpsurfAG5-273.15-ExpTemp)*100/(ExpTemp+273)

        print("Creating lineplot of Tavg over time for ArrGen CompareExp...")
        pyplot.figure()
        pyplot.plot(ExpTime,ExpTemp,'-',linewidth=3, label="Experimental")
        pyplot.plot(tsampleAG5/60,AG5_FDavg-273,'-.',linewidth=3, label="N=60 dt=0.2 avg")
        pyplot.plot(tsampleAG5/60,AG5_T_surf-273,'--',linewidth=3, label="N=60 dt=0.2 r=R")
        bounds=np.array([0, 200])
        pyplot.ylim(bounds) 
        pyplot.xlabel('t (min)', fontsize=10)
        pyplot.ylabel('Temperature (C)', fontsize=10)
        pyplot.title('Change of average temperature')
        pyplot.legend(loc='best')
        pyplot.xticks(np.arange(0, tsampleAG5[tsampleAG5.shape[0]-1]/60, step=5))
        pyplot.savefig('ArrGen_CompareExp_Tavg.pdf', bbox_inches='tight')

        print("Creating Error lineplot of Tavg over time for ArrGen CompareExp...")
        pyplot.figure()
        pyplot.plot(ExpTime,ErrExp_avg,'--',linewidth=3, label="Tavg")
        pyplot.plot(ExpTime,ErrExp_surf,'-',linewidth=3, label="Tsurf")
        #pyplot.plot(tsampleAG5/60,Err9_7,':',linewidth=3, label="Err9_7")

#        pyplot.plot(Expdataarray[:,0],Expdataarray[:,1],'-.',linewidth=3,label="Experimental")
        #bounds=np.array([0, 10])
        #pyplot.ylim(bounds) 
        pyplot.xlabel('t (min)', fontsize=10)
        pyplot.ylabel('% Error', fontsize=10)
        pyplot.title('Percent error of average temperature with respect to the experimental data')
        pyplot.legend(loc='best')
        pyplot.xticks(np.arange(0, tsampleAG5[tsampleAG5.shape[0]-1]/60, step=5))
        pyplot.savefig('ArrGen_CompareExpErr_Tavg.pdf', bbox_inches='tight')
    if rdist:
        tplot=2000
        zplot=H/2

        AG5_param_data=pd.read_csv('AG5/AG5_parameters.csv',header=None)
        AG5_param=np.array(AG5_param_data)
        NcellrAG5=int(AG5_param[0])
        NcellzAG5=int(AG5_param[1])
        deltAG5=AG5_param[2]
        tfinalAG5=AG5_param[3]
        print(deltAG5)
        print(tfinalAG5)
        print(NcellrAG5)
        print(NcellzAG5)
        tsampleAG5=np.linspace(0,tfinalAG5,tfinalAG5/deltAG5+1)
        delr=R/NcellrAG5
        delz=H/NcellzAG5
        # Grid numbering (Row increase => r increase, Column increase => z increase)
        GridMap=np.zeros((NcellzAG5,NcellrAG5),dtype=np.int)
        for i in range(0,NcellzAG5):
            GridMap[i,:]=(i*NcellrAG5)+np.arange(0,NcellrAG5,1)
        GridMapFD=GridMap
        GridMapFD=np.transpose(GridMapFD)

        # Locations of grid points (Consistent with original GridMap, r increases in row direction, z increases in column direction)
        zLoc=np.zeros((NcellzAG5,NcellrAG5))
        rLoc=np.zeros((NcellzAG5,NcellrAG5))
        for m in range(0,NcellzAG5):
            for n in range(0,NcellrAG5):
                rLoc[m,n]=delr*(n+0.5)
                zLoc[m,n]=delz*(m+0.5)
        AG5_T_FD_all_data=pd.read_csv('AG5/AG5_FD_all.csv',header=None)
        AG5_T_FD_all=np.array(AG5_T_FD_all_data)

        #AG5_T_Analy_all_data=pd.read_csv('AG5/AG5_Analy_all.csv',header=None)
        #AG5_T_Analy_all=np.array(AG5_T_Analy_all_data)
        AG5_T_Analy_all=np.empty((NcellrAG5*NcellzAG5,tsampleAG5.shape[0]))*np.nan
        lineplot.alongr_plotT(AG5_T_FD_all,AG5_T_Analy_all,rLoc,zLoc,tsampleAG5,R,delz,deltAG5,NcellrAG5,NcellzAG5,zplot,tplot)
    if zdist:
        tplot=2000
        rplot=R/2

        AG5_param_data=pd.read_csv('AG5/AG5_parameters.csv',header=None)
        AG5_param=np.array(AG5_param_data)
        NcellrAG5=int(AG5_param[0])
        NcellzAG5=int(AG5_param[1])
        deltAG5=AG5_param[2]
        tfinalAG5=AG5_param[3]
        print(deltAG5)
        print(tfinalAG5)
        print(NcellrAG5)
        print(NcellzAG5)
        tsampleAG5=np.linspace(0,tfinalAG5,tfinalAG5/deltAG5+1)
        delr=R/NcellrAG5
        delz=H/NcellzAG5
        # Grid numbering (Row increase => r increase, Column increase => z increase)
        GridMap=np.zeros((NcellzAG5,NcellrAG5),dtype=np.int)
        for i in range(0,NcellzAG5):
            GridMap[i,:]=(i*NcellrAG5)+np.arange(0,NcellrAG5,1)
        GridMapFD=GridMap
        GridMapFD=np.transpose(GridMapFD)

        # Locations of grid points (Consistent with original GridMap, r increases in row direction, z increases in column direction)
        zLoc=np.zeros((NcellzAG5,NcellrAG5))
        rLoc=np.zeros((NcellzAG5,NcellrAG5))
        for m in range(0,NcellzAG5):
            for n in range(0,NcellrAG5):
                rLoc[m,n]=delr*(n+0.5)
                zLoc[m,n]=delz*(m+0.5)
        AG5_T_FD_all_data=pd.read_csv('AG5/AG5_FD_all.csv',header=None)
        AG5_T_FD_all=np.array(AG5_T_FD_all_data)

        #AG5_T_Analy_all_data=pd.read_csv('AG5/AG5_Analy_all.csv',header=None)
        #AG5_T_Analy_all=np.array(AG5_T_Analy_all_data)
        AG5_T_Analy_all=np.empty((NcellrAG5*NcellzAG5,tsampleAG5.shape[0]))*np.nan
        lineplot.alongz_plotT(AG5_T_FD_all,AG5_T_Analy_all,rLoc,zLoc,tsampleAG5,H,delr,deltAG5,NcellrAG5,NcellzAG5,rplot,tplot)


if NoGen:
    if GridIndependence:

        NG5_param_data=pd.read_csv('NG5_P_GMRES/NG5_parameters.csv',header=None)
        NG5_param=np.array(NG5_param_data)
        deltNG5=NG5_param[2]
        tfinalNG5=NG5_param[3]
        tsampleNG5=np.linspace(0,tfinalNG5,tfinalNG5/deltNG5+1)

        #NG6_param_data=pd.read_csv('NG6/NG6_parameters.csv',header=None)
        #NG6_param=np.array(NG6_param_data)
        #deltNG6=NG6_param[2]
        #tfinalNG6=NG6_param[3]
        #tsampleNG6=np.linspace(0,tfinalNG6,tfinalNG6/deltNG6+1)

        NG7_param_data=pd.read_csv('NG7_P_GMRES/NG7_parameters.csv',header=None)
        NG7_param=np.array(NG7_param_data)
        deltNG7=NG7_param[2]
        tfinalNG7=NG7_param[3]
        tsampleNG7=np.linspace(0,tfinalNG7,tfinalNG7/deltNG7+1)

        #NG8_param_data=pd.read_csv('NG8/NG8_parameters.csv',header=None)
        #NG8_param=np.array(NG8_param_data)
        #deltNG8=NG8_param[2]
        #tfinalNG8=NG8_param[3]
        #tsampleNG8=np.linspace(0,tfinalNG8,tfinalNG8/deltNG8+1)

        NG9_param_data=pd.read_csv('NG9_P_GMRES/NG9_parameters.csv',header=None)
        NG9_param=np.array(NG9_param_data)
        deltNG9=NG9_param[2]
        tfinalNG9=NG9_param[3]
        tsampleNG9=np.linspace(0,tfinalNG9,tfinalNG9/deltNG9+1)

        #NG10_param_data=pd.read_csv('NG10_P_GMRES/NG10_parameters.csv',header=None)
        #NG10_param=np.array(NG10_param_data)
        #deltNG10=NG10_param[2]
        #tfinalNG10=NG10_param[3]
        #tsampleNG10=np.linspace(0,tfinalNG10,tfinalNG10/deltNG10+1)
        #print('Just before reading NG10')
        NG5_T_FD_all_data=pd.read_csv('NG5_P_GMRES/NG5_FD_all.csv',header=None)
        NG5_T_FD_all=np.array(NG5_T_FD_all_data)
        #NG6_T_FD_all_data=pd.read_csv('NG6/NG6_FD_all.csv',header=None)
        #NG6_T_FD_all=np.array(NG6_T_FD_all_data)
        NG7_T_FD_all_data=pd.read_csv('NG7_P_GMRES/NG7_FD_all.csv',header=None)
        NG7_T_FD_all=np.array(NG7_T_FD_all_data)
        #NG8_T_FD_all_data=pd.read_csv('NG8/NG8_FD_all.csv',header=None)
        #NG8_T_FD_all=np.array(NG8_T_FD_all_data)
        NG9_T_FD_all_data=pd.read_csv('NG9_P_GMRES/NG9_FD_all.csv',header=None)
        NG9_T_FD_all=np.array(NG9_T_FD_all_data)
        #print('Just before reading NG10')
        #NG10_T_FD_all_data=pd.read_csv('NG10_P_GMRES/NG10_FD_all.csv',header=None)
        #NG10_T_FD_all=np.array(NG10_T_FD_all_data)
        #print('Just after reading NG10')

        NG5_T_Analy_all_data=pd.read_csv('NG5/NG5_Analy_all.csv',header=None)
        NG5_T_Analy_all=np.array(NG5_T_Analy_all_data)
        #NG6_T_Analy_all_data=pd.read_csv('NG6/NG6_Analy_all.csv',header=None)
        #NG6_T_Analy_all=np.array(NG6_T_Analy_all_data)
        NG7_T_Analy_all_data=pd.read_csv('NG7/NG7_Analy_all.csv',header=None)
        NG7_T_Analy_all=np.array(NG7_T_Analy_all_data)
        #NG8_T_Analy_all_data=pd.read_csv('NG8/NG8_Analy_all.csv',header=None)
        #NG8_T_Analy_all=np.array(NG8_T_Analy_all_data)
        NG9_T_Analy_all_data=pd.read_csv('NG9/NG9_Analy_all.csv',header=None)
        NG9_T_Analy_all=np.array(NG9_T_Analy_all_data)
        #NG10_T_Analy_all_data=pd.read_csv('NG10_P_GMRES/NG10_Analy_all.csv',header=None)
        #NG10_T_Analy_all=np.array(NG10_T_Analy_all_data)
        #print('Just after reading NG10')
        Diff5=np.abs(NG5_T_FD_all-NG5_T_Analy_all)/NG5_T_Analy_all
        #Diff6=np.abs(NG6_T_FD_all-NG6_T_Analy_all)/NG6_T_Analy_all
        Diff7=np.abs(NG7_T_FD_all-NG7_T_Analy_all)/NG7_T_Analy_all
        #Diff8=np.abs(NG8_T_FD_all-NG8_T_Analy_all)/NG8_T_Analy_all
        Diff9=np.abs(NG9_T_FD_all-NG9_T_Analy_all)/NG9_T_Analy_all
        #Diff10=np.abs(NG10_T_FD_all-NG10_T_Analy_all)
 
        Maxerr5=np.reshape(Diff5.max(axis=0),np.shape(tsampleNG5))
        #Maxerr6=np.reshape(Diff6.max(axis=0),np.shape(tsampleNG6))
        Maxerr7=np.reshape(Diff7.max(axis=0),np.shape(tsampleNG7))        
        #Maxerr8=np.reshape(Diff8.max(axis=0),np.shape(tsampleNG8))
        Maxerr9=np.reshape(Diff9.max(axis=0),np.shape(tsampleNG9))
        #Maxerr10=np.reshape(Diff10.max(axis=0),np.shape(tsampleNG10))
        #print(Diff7[:,1000])
        
        print("Creating MaxError lineplot over time for NoGen GridIndependence...")
        pyplot.figure()
        #pyplot.plot(tsampleNG10/60,Maxerr10,'-',linewidth=3, label="N=80 dt=0.1")
        pyplot.plot(tsampleNG5/60,100*Maxerr5,'-',linewidth=3, label="N=60 dt=0.2")
        #pyplot.plot(tsampleNG6/60,Maxerr6,'--',linewidth=3, label="N=50 dt=0.3")
        pyplot.plot(tsampleNG7/60,100*Maxerr7,'-.',linewidth=3, label="N=40 dt=0.5")
        #pyplot.plot(tsampleNG8/60,Maxerr8,':',linewidth=3, label="N=30 dt=0.75")
        pyplot.plot(tsampleNG9/60,100*Maxerr9,'-',linewidth=3, label="N=20 dt=1")

#        pyplot.plot(Expdataarray[:,0],Expdataarray[:,1],'-.',linewidth=3,label="Experimental")
        #bounds=np.array([0, 10])
        #pyplot.ylim(bounds) 
        pyplot.xlabel('t (min)', fontsize=10)
        pyplot.ylabel('% Error', fontsize=10)
        pyplot.title('Max % Error in the domain over time')
        pyplot.legend(loc='best')
        pyplot.xticks(np.arange(0, tsampleNG7[tsampleNG7.shape[0]-1]/60, step=5))
        pyplot.savefig('NoGen_GridIndepErr.pdf', bbox_inches='tight')

    if LinearSolver:
               # Read parameters
        NG2_param_data=pd.read_csv('NG2/NG2_parameters.csv',header=None)
        NG2_param=np.array(NG2_param_data)
        
        delt=NG2_param[2]
        tfinal=NG2_param[3]
        print(delt)
        print(tfinal)
        tsample=np.linspace(0,tfinal,tfinal/delt+1)

        NG8_FDavg_data=pd.read_csv('NG8/NG8_FDavg.csv',header=None)
        NG8_FDavg=np.array(NG8_FDavg_data)

        NG2_FDavg_data=pd.read_csv('NG2/NG2_FDavg.csv',header=None)
        NG2_FDavg=np.array(NG2_FDavg_data)

        NG3_FDavg_data=pd.read_csv('NG3/NG3_FDavg.csv',header=None)
        NG3_FDavg=np.array(NG3_FDavg_data)

        NG4_FDavg_data=pd.read_csv('NG4/NG4_FDavg.csv',header=None)
        NG4_FDavg=np.array(NG4_FDavg_data)

        NG8_T_FD_all_data=pd.read_csv('NG8/NG8_FD_all.csv',header=None)
        NG8_T_FD_all=np.array(NG8_T_FD_all_data)
        NG2_T_FD_all_data=pd.read_csv('NG2/NG2_FD_all.csv',header=None)
        NG2_T_FD_all=np.array(NG2_T_FD_all_data)
        NG3_T_FD_all_data=pd.read_csv('NG3/NG3_FD_all.csv',header=None)
        NG3_T_FD_all=np.array(NG3_T_FD_all_data)
        NG4_T_FD_all_data=pd.read_csv('NG4/NG4_FD_all.csv',header=None)
        NG4_T_FD_all=np.array(NG4_T_FD_all_data)

        Diff2=np.abs(NG2_T_FD_all-NG8_T_FD_all)/NG8_T_FD_all
        Diff3=np.abs(NG3_T_FD_all-NG8_T_FD_all)/NG8_T_FD_all
        Diff4=np.abs(NG4_T_FD_all-NG8_T_FD_all)/NG8_T_FD_all

        Maxerr2=np.reshape(Diff2.max(axis=0),np.shape(NG2_FDavg_data))
        Maxerr3=np.reshape(Diff3.max(axis=0),np.shape(NG2_FDavg_data))
        Maxerr4=np.reshape(Diff4.max(axis=0),np.shape(NG2_FDavg_data)) 

	#print(NG4_FDavg)
        print("Creating lineplot of Tavg over time for NoGen LinearSolver...")
        pyplot.figure()
        pyplot.plot(tsample/60,NG8_FDavg-273,'--',linewidth=3, label="LU solver")
        pyplot.plot(tsample/60,NG2_FDavg-273,':',linewidth=3, label="GMRES tol 1e-4")
        pyplot.plot(tsample/60,NG3_FDavg-273,'-.',linewidth=3, label="GMRES tol 1e-5")
        pyplot.plot(tsample/60,NG4_FDavg-273,':',linewidth=3, label="GMRES tol 1e-6")

#        pyplot.plot(Expdataarray[:,0],Expdataarray[:,1],'-.',linewidth=3,label="Experimental")

        bounds=np.array([0, 200])
        pyplot.ylim(bounds) 
        pyplot.xlabel('t (min)', fontsize=10)
        pyplot.ylabel('Temperature (C)', fontsize=10)
        pyplot.title('Change of average temperature')
        pyplot.legend(loc='best')
        pyplot.xticks(np.arange(0, tsample[tsample.shape[0]-1]/60, step=5))
        pyplot.savefig('NoGen_LinSol_Tavg.pdf', bbox_inches='tight')


        print("Creating MaxError lineplot over time for NoGen LinearSolver...")
        pyplot.figure()
        pyplot.plot(tsample/60,100*Maxerr2,'-',linewidth=3, label="GMRES 1e-4")
        pyplot.plot(tsample/60,100*Maxerr3,'--',linewidth=3, label="GMRES 1e-5")
        pyplot.plot(tsample/60,100*Maxerr4,'-.',linewidth=3, label="GMRES 1e-6")

#        pyplot.plot(Expdataarray[:,0],Expdataarray[:,1],'-.',linewidth=3,label="Experimental")
        #bounds=np.array([0, 10])
        #pyplot.ylim(bounds) 
        pyplot.xlabel('t (min)', fontsize=10)
        pyplot.ylabel('% Error', fontsize=10)
        pyplot.title('Max Error in the domain relative to LU solver')
        pyplot.legend(loc='best')
        pyplot.xticks(np.arange(0, tsample[tsample.shape[0]-1]/60, step=5))
        pyplot.savefig('NoGen_LinSolErr.pdf', bbox_inches='tight')

    if Tavgcomparison:
        # Read parameters
        NG5_param_data=pd.read_csv('NG5_P_GMRES/NG5_parameters.csv',header=None)
        NG5_param=np.array(NG5_param_data)
        
        delt=NG5_param[2]
        tfinal=NG5_param[3]
        print(delt)
        print(tfinal)
        tsample=np.linspace(0,tfinal,tfinal/delt+1)

        NG5_Analyavg_data=pd.read_csv('NG5/NG5_Analyavg.csv',header=None)
        NG5_Analyavg=np.array(NG5_Analyavg_data)

        NG5_FDavg_data=pd.read_csv('NG5_P_GMRES/NG5_FDavg.csv',header=None)
        NG5_FDavg=np.array(NG5_FDavg_data)

        print("Creating lineplot of Tavg over time for best case of NoGen...")
        pyplot.figure()
        pyplot.plot(tsample/60,NG5_Analyavg-273,'-',linewidth=3, label="Analytical")
        pyplot.plot(tsample/60,NG5_FDavg-273,'--',linewidth=3, label="Finite Difference")


#        pyplot.plot(Expdataarray[:,0],Expdataarray[:,1],'-.',linewidth=3,label="Experimental")

        bounds=np.array([0, 200])
        pyplot.ylim(bounds) 
        pyplot.xlabel('t (min)', fontsize=10)
        pyplot.ylabel('Temperature (C)', fontsize=10)
        pyplot.title('Change of average temperature')
        pyplot.legend(loc='best')
        pyplot.xticks(np.arange(0, tsample[tsample.shape[0]-1]/60, step=5))
        pyplot.savefig('NoGen_best_Tavg.pdf', bbox_inches='tight')
    if rdist:
        tplot=2000
        zplot=H/2

        NG5_param_data=pd.read_csv('NG5_P_GMRES/NG5_parameters.csv',header=None)
        NG5_param=np.array(NG5_param_data)
        NcellrNG5=int(NG5_param[0])
        NcellzNG5=int(NG5_param[1])
        deltNG5=NG5_param[2]
        tfinalNG5=NG5_param[3]
        print(deltNG5)
        print(tfinalNG5)
        print(NcellrNG5)
        print(NcellzNG5)
        tsampleNG5=np.linspace(0,tfinalNG5,tfinalNG5/deltNG5+1)
        delr=R/NcellrNG5
        delz=H/NcellzNG5
        # Grid numbering (Row increase => r increase, Column increase => z increase)
        GridMap=np.zeros((NcellzNG5,NcellrNG5),dtype=np.int)
        for i in range(0,NcellzNG5):
            GridMap[i,:]=(i*NcellrNG5)+np.arange(0,NcellrNG5,1)
        GridMapFD=GridMap
        GridMapFD=np.transpose(GridMapFD)

        # Locations of grid points (Consistent with original GridMap, r increases in row direction, z increases in column direction)
        zLoc=np.zeros((NcellzNG5,NcellrNG5))
        rLoc=np.zeros((NcellzNG5,NcellrNG5))
        for m in range(0,NcellzNG5):
            for n in range(0,NcellrNG5):
                rLoc[m,n]=delr*(n+0.5)
                zLoc[m,n]=delz*(m+0.5)
        NG5_T_FD_all_data=pd.read_csv('NG5_P_GMRES/NG5_FD_all.csv',header=None)
        NG5_T_FD_all=np.array(NG5_T_FD_all_data)

        NG5_T_Analy_all_data=pd.read_csv('NG5/NG5_Analy_all.csv',header=None)
        NG5_T_Analy_all=np.array(NG5_T_Analy_all_data)

        lineplot.alongr_plotT(NG5_T_FD_all,NG5_T_Analy_all,rLoc,zLoc,tsampleNG5,R,delz,deltNG5,NcellrNG5,NcellzNG5,zplot,tplot)
    if zdist:
        tplot=2000
        rplot=R/2

        NG5_param_data=pd.read_csv('NG5_P_GMRES/NG5_parameters.csv',header=None)
        NG5_param=np.array(NG5_param_data)
        NcellrNG5=int(NG5_param[0])
        NcellzNG5=int(NG5_param[1])
        deltNG5=NG5_param[2]
        tfinalNG5=NG5_param[3]
        print(deltNG5)
        print(tfinalNG5)
        print(NcellrNG5)
        print(NcellzNG5)
        tsampleNG5=np.linspace(0,tfinalNG5,tfinalNG5/deltNG5+1)
        delr=R/NcellrNG5
        delz=H/NcellzNG5
        # Grid numbering (Row increase => r increase, Column increase => z increase)
        GridMap=np.zeros((NcellzNG5,NcellrNG5),dtype=np.int)
        for i in range(0,NcellzNG5):
            GridMap[i,:]=(i*NcellrNG5)+np.arange(0,NcellrNG5,1)
        GridMapFD=GridMap
        GridMapFD=np.transpose(GridMapFD)

        # Locations of grid points (Consistent with original GridMap, r increases in row direction, z increases in column direction)
        zLoc=np.zeros((NcellzNG5,NcellrNG5))
        rLoc=np.zeros((NcellzNG5,NcellrNG5))
        for m in range(0,NcellzNG5):
            for n in range(0,NcellrNG5):
                rLoc[m,n]=delr*(n+0.5)
                zLoc[m,n]=delz*(m+0.5)
        NG5_T_FD_all_data=pd.read_csv('NG5_P_GMRES/NG5_FD_all.csv',header=None)
        NG5_T_FD_all=np.array(NG5_T_FD_all_data)

        NG5_T_Analy_all_data=pd.read_csv('NG5/NG5_Analy_all.csv',header=None)
        NG5_T_Analy_all=np.array(NG5_T_Analy_all_data)

        lineplot.alongz_plotT(NG5_T_FD_all,NG5_T_Analy_all,rLoc,zLoc,tsampleNG5,H,delr,deltNG5,NcellrNG5,NcellzNG5,rplot,tplot)
    if Analytruncation:
        NG10_param_data=pd.read_csv('NG10_P_GMRES/NG10_parameters.csv',header=None)
        NG10_param=np.array(NG10_param_data)
        deltNG10=NG10_param[2]
        tfinalNG10=NG10_param[3]
        tsampleNG10=np.linspace(0,tfinalNG10,tfinalNG10/deltNG10+1)

        NG10_Analyerr_data=pd.read_csv('NG10_P_GMRES/NG10_analyErr_all.csv',header=None)
        NG10_Analyerr=np.array(NG10_Analyerr_data)
        NG10_Analyerr_69=np.reshape(NG10_Analyerr[:,0],np.shape(tsampleNG10))
        NG10_Analyerr_79=np.reshape(NG10_Analyerr[:,1],np.shape(tsampleNG10))
        NG10_Analyerr_89=np.reshape(NG10_Analyerr[:,2],np.shape(tsampleNG10))

        print("Creating Truncation error for lineplot over time for NoGen Analytruncation...")
        pyplot.figure()
        pyplot.plot(tsampleNG10/60,NG10_Analyerr_69,'-',linewidth=3, label="6 term")
        pyplot.plot(tsampleNG10/60,NG10_Analyerr_79,'-.',linewidth=3, label="7 term")
        pyplot.plot(tsampleNG10/60,NG10_Analyerr_89,'--',linewidth=3, label="8 term")
        
#        pyplot.plot(Expdataarray[:,0],Expdataarray[:,1],'-.',linewidth=3,label="Experimental")
        bounds=np.array([0, 0.2])
        pyplot.xlim(bounds) 
        pyplot.xlabel('t (min)', fontsize=10)
        pyplot.ylabel('% Max Error in Temp', fontsize=10)
        pyplot.title('Max % Error in the domain over time')
        pyplot.legend(loc='best')
        #pyplot.xticks(np.arange(0, tsampleNG10[tsampleNG10.shape[0]-1]/60, step=5))
        pyplot.savefig('NoGen_Analytruncation.pdf', bbox_inches='tight')
       

if Contourplot:
    tplot=2200
    NG5_param_data=pd.read_csv('NG5/NG5_parameters.csv',header=None)
    NG5_param=np.array(NG5_param_data)
    NcellrNG5=int(NG5_param[0])
    NcellzNG5=int(NG5_param[1])
    deltNG5=NG5_param[2]
    tfinalNG5=NG5_param[3]
    print(deltNG5)
    print(tfinalNG5)
    print(NcellrNG5)
    print(NcellzNG5)
    tsampleNG5=np.linspace(0,tfinalNG5,tfinalNG5/deltNG5+1)
    delr=R/NcellrNG5
    delz=H/NcellzNG5
    # Grid numbering (Row increase => r increase, Column increase => z increase)
    GridMap=np.zeros((NcellzNG5,NcellrNG5),dtype=np.int)
    for i in range(0,NcellzNG5):
        GridMap[i,:]=(i*NcellrNG5)+np.arange(0,NcellrNG5,1)
        GridMapFD=GridMap
        GridMapFD=np.transpose(GridMapFD)

        # Locations of grid points (Consistent with original GridMap, r increases in row direction, z increases in column direction)
    zLoc=np.zeros((NcellzNG5,NcellrNG5))
    rLoc=np.zeros((NcellzNG5,NcellrNG5))
    for m in range(0,NcellzNG5):
        for n in range(0,NcellrNG5):
            rLoc[m,n]=delr*(n+0.5)
            zLoc[m,n]=delz*(m+0.5)
    NG5_T_FD_all_data=pd.read_csv('NG5/NG5_FD_all.csv',header=None)
    NG5_T_FD_all=np.array(NG5_T_FD_all_data)

    AG5_T_FD_all_data=pd.read_csv('AG5/AG5_FD_all.csv',header=None)
    AG5_T_FD_all=np.array(AG5_T_FD_all_data)

    contourplot.compareanaly(NG5_T_FD_all,AG5_T_FD_all,rLoc,zLoc,tsampleNG5,deltNG5,NcellrNG5,NcellzNG5,tplot)
if Animation:

    NG5_param_data=pd.read_csv('NG5/NG5_parameters.csv',header=None)
    NG5_param=np.array(NG5_param_data)
    NcellrNG5=int(NG5_param[0])
    NcellzNG5=int(NG5_param[1])
    deltNG5=NG5_param[2]
    tfinalNG5=NG5_param[3]
    #print(deltNG5)
    #print(tfinalNG5)
    #print(NcellrNG5)
    #print(NcellzNG5)
    tsampleNG5=np.linspace(0,tfinalNG5,tfinalNG5/deltNG5+1)
    print(int((tsampleNG5.shape[0]-1)/100))
    delr=R/NcellrNG5
    delz=H/NcellzNG5
    # Grid numbering (Row increase => r increase, Column increase => z increase)
    GridMap=np.zeros((NcellzNG5,NcellrNG5),dtype=np.int)
    for i in range(0,NcellzNG5):
        GridMap[i,:]=(i*NcellrNG5)+np.arange(0,NcellrNG5,1)
        GridMapFD=GridMap
        GridMapFD=np.transpose(GridMapFD)

        # Locations of grid points (Consistent with original GridMap, r increases in row direction, z increases in column direction)
    zLoc=np.zeros((NcellzNG5,NcellrNG5))
    rLoc=np.zeros((NcellzNG5,NcellrNG5))
    for m in range(0,NcellzNG5):
        for n in range(0,NcellrNG5):
            rLoc[m,n]=delr*(n+0.5)
            zLoc[m,n]=delz*(m+0.5)
    NG5_T_FD_all_data=pd.read_csv('NG5/NG5_FD_all.csv',header=None)
    NG5_T_FD_all=np.array(NG5_T_FD_all_data)

    AG5_T_FD_all_data=pd.read_csv('AG5/AG5_FD_all.csv',header=None)
    AG5_T_FD_all=np.array(AG5_T_FD_all_data)

    createanimation.compareanaly(NG5_T_FD_all,AG5_T_FD_all,rLoc,zLoc,tsampleNG5,NcellrNG5,NcellzNG5)

if printtimes:
    NG5_param_data=pd.read_csv('NG5/NG5_parameters.csv',header=None)
    NG5_param=np.array(NG5_param_data)
    NG5_tFD=NG5_param[4]
    NG5_tAnaly=NG5_param[5]