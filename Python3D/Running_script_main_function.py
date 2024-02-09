
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 09:07:30 2022

@author: Andrea Piccolo
Slab Avalanching small project

1) Function that reads the vtk file and produce picture
    a) Read vtk, plot phase, topography 

"""

import sys,os
import numpy as np
from time import perf_counter 
import getopt
import argparse
from Save_Data_Base import _write_h5_database
    

"""
Import useful function 
"""
from Read_VTK_files_LAMEM import * #Read_VTK_files_LAMEM

from Read_VTK_files_LAMEM import  _file_list

from Parser_File import * 
from Slab_detachment import * 

def _run_script_visualization(ptsave,Folder,Test_Name,l_path,Data_Base_path,Group):
    t_INIZIO = perf_counter()

    
    dic_val= {
    "OP"       : "_Oceanic_Lit_ [ ]",
    "CC1"       : "_Continental_Crust_ [ ]",
    "CC2"      : "_Continental_Crust_2 [ ]",
    "Sed"      : "_Sediments_ [ ]",
    "Lit"      : "_Lithosphere_all_ [ ]",
    "T"        : "temperature [C]",
    "velocity" : "velocity [cm/yr]",
    "disp"     : "tot_displ [km]",
    "vis"      : "visc_creep [Pa*s]",
    "nu"       : "visc_total [Pa*s]",
    "eps"      : "j2_strain_rate [1/s]",
    "tau"      : "j2_dev_stress [MPa]",
    "Rho"      : "density [kg/m^3]",
    "Psi"      : "plast_dissip [W/m^3]"}
    
    phase_dictionary={
    "0" : ["Air","white"],
    "1" : ["Upper Crust", "silver"],
    "2" : ["Lower Crust", "dimgray"],
    "3" : ["Continental Lithosperhic mantle","lightgoldenrodyellow"],
    "4" : ["Continental Lithosperhic mantle2","palegoldenrod"],
    "5" : ["Upper Mantle", "aliceblue"],
    "6" : ["Slab","sandybrown"],
    "7" : ["Oceanic Crust","palegreen"],
    "8" : ["Weak Zone Mantle","rosybrown"],
    "9" : ["Upper Crust2","linen"],
    "10" : ["Lower  Crust 2","gainsboro"],
    "11" : ["Pelagic Sediment","lavender"],
    "12" : ["Prism","burlywood"],
    "13" : ["Passive Margin","darkkhaki"],
    "14" : ["Mantle 410 km","palegreen"],
    "15" : ["Slab 410 km","peachpuff"],
    "16" : ["Eclogite","palevioletred"],
    "17" : ["Lower  Mantle","darkseagreen"],
    "18" : ["Lower Slab","tan"]
    }
    
    fname=("SDet")
    dyn             ='%s.pvtr'      %fname
    surf            ='%s_surf.pvts' %fname
    phase           ='%s_phase.pvtr'%fname
    ptrac ='%s_passive_tracers.pvtu' %fname
    ###############################################
    fname2 = fname+".pvd"
    fname=os.path.join(Folder,Test_Name,fname2)
    time, Flist, n_tim_step =_file_list(fname)    # Retrive the file list and the time step information
    
    #nstp = 41
    #time = time[0:nstp]
    #Flist = Flist[0:nstp]
    
    ##############################################
    ts0 = Flist[0]
    fn=ts0[1:ts0.find('/')]
    Filename_0 = [os.path.join(Folder,Test_Name,fn,dyn),os.path.join(Folder,Test_Name,fn,phase)]
    Filename_0ptr = os.path.join(Folder,Test_Name,fn,ptrac)
    folder_ = os.path.join(Folder,Test_Name)
    IG       = Initial_Geometry(os.path.join(folder_,'Test_Data_Base.mat'))
    C = Coordinate_System(Filename_0,ptsave,(-700.0,700.0),(-500,500),(-600.0,50.0))
    Phase_DB = Phase_Data_Base(folder_)
    Initial_Condition = Initial_condition(Phase_DB.Phase_6_,Phase_DB.Phase_5_,IG)
    Initial_Condition.tc = Initial_Condition.tc/3.5
    ###############################################
    FSurf = Free_S_Slab_break_off(C,len(time)) # Create the instance of free surface class 
    DYN   = VAL(C,dic_val)  # Create the instance of the .pvd file 
    Ph    = Phase_det(C,phase_dictionary) # Create the instance of the Phase field
    Slab  = SLAB(C,IG,len(time))        # Create the instance of the Slab.
    PT    = Passive_Tracers(Filename_0ptr) # Create the instance of Passive tracers 
    ################################################
    ###############################################################################
    # Read Information related to the initial condition
    ###############################################################################
    ipic= 0 
    for istp in Flist:#[19:]:
    ########################### Files and time stepping information############
        t1_start = perf_counter()
        fn=istp[1:istp.find('/')]
        t_cur=time[ipic]        
        if ipic>0:
            dt = time[ipic]-time[ipic-1]
        else:
            dt = 1.0 
        print('Timestep n', str(ipic))
            
        Filename_dyn=os.path.join(Folder,Test_Name,fn,dyn)
        Filename_ph=os.path.join(Folder,Test_Name,fn,phase)
        Filename_s=os.path.join(Folder,Test_Name,fn,surf)
        Filename_ptr = os.path.join(Folder,Test_Name,fn,ptrac)
        ##############Passive Tracers##############################################
        t1 = perf_counter()
        PT._update_PTracer(Filename_ptr)
        t2 = perf_counter()
        Values_time = t2-t1
        print("PT value took","{:02}".format(Values_time))
        
        ###### Retrieve the field that is needed ##################################
        t1 = perf_counter()
        DYN._update_Val(Filename_dyn,C,Initial_Condition)
        t2 = perf_counter()
        Values_time = t2-t1
        print("Dynamic value took","{:02}".format(Values_time))
        ###########################################################################
        t1 = perf_counter()
        Ph._update_phase(Filename_ph,C)  
        Ph._interpolate_dyn_phase(DYN,C)
        t2 = perf_counter()    
        average_plot_time = t2-t1
        print("Phase_update","{:02}".format(average_plot_time))

        ###########################################################################
        t1= perf_counter()
        FSurf._Update_(Filename_s,C,ipic)
        FSurf._update_extra_variables(DYN,C,dt,ipic)
        FSurf.ASCI_FILE_ALT(ipic,t_cur,Test_Name,ptsave,C)
        t2 = perf_counter()
        print("Free surface ","{:02}".format(t2-t1))
        #======Select the chosen ===========# 
        list_phase = [1,2,12,13,16,17]
        levels     = [-15.0,-14.0]
        if ipic == 0:
            BSPT = Basement_Passive_Tracer(PT,C,FSurf,list_phase,levels,len(time),ipic)
            #BSPT._plot_passive_tracers(FSurf,C.x,C.y,time,ipic,ptsave,'T')
        else:
            BSPT._update_PTDB(PT,ipic,time)
            #BSPT._plot_passive_tracers(FSurf,C.x,C.y,time,ipic,ptsave,'dzdt')
            #BSPT._plot_passive_tracers(FSurf,C.x,C.y,time,ipic,ptsave,'dPdt')
            #BSPT._plot_passive_tracers(FSurf,C.x,C.y,time,ipic,ptsave,'dTdt')

            

        ###########################################################################
        t1 = perf_counter()
        Slab. _update_C(C,FSurf,Ph,IG,ipic,t_cur,dt)
        Slab._plot_average_C(t_cur,C.xp,C.zp,ptsave,ipic,IG.Slab,Initial_Condition,time)  
        #FSurf._plot_maps_FS(t_cur,C.y,C.x,ptsave,ipic,Slab)
        #FSurf._plot_1D_plots_Free_surface(ipic,ptsave,Slab,t_cur,Initial_Condition.D0/1e3)
        t2 = perf_counter()
        print("Slab routine ","{:02}".format(t2-t1))
        ###########################################################################
        t1 = perf_counter()   
        plot_map_time = t2-t1
        ###########################################################################
        #  Free Surface 
        #######################################################################
        FS_time = t2-t1
        t1_end = perf_counter()
        tstep_time_pr = t1_end-t1_start
        minutes, seconds = divmod(t1_end-t1_start, 60)
        ipic +=1 
        minutes, seconds = divmod(tstep_time_pr, 60)
        print("===========================================")
        print("||Script took ","{:02}".format(int(minutes)),"minutes and","{:05.2f}".format(seconds),' seconds ||' )
        print("===========================================")
    _write_h5_database(Data_Base_path,Group,Test_Name,Slab,Initial_Condition,IG,C,FSurf,Phase_DB,BSPT,time)

