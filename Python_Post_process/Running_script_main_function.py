
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

    

"""
Import useful function 
"""
from Read_VTK_files_LAMEM import * #Read_VTK_files_LAMEM

from Read_VTK_files_LAMEM import  _file_list

from Parser_File import * 
from Slab_detachment import * 
from Slab_detachment import _plot_D_D0

def _run_script_visualization(ptsave,Folder,Test_Name,l_path,vIC):
    t_INIZIO = perf_counter()

    
    dic_val= {
    "OP"       : "_Oceanic_Lit_ [ ]",
    "CC1"       : "_Continental_Crust_ [ ]",
    "T"        : "temperature [C]",
    "velocity" : "velocity [cm/yr]",
    "disp"     : "tot_displ [km]",
    "vis"      : "visc_total [Pa*s]",
    "nu"       : "visc_creep [Pa*s]",
    "eps"      : "j2_strain_rate [1/s]",
    "tau"      : "j2_dev_stress [MPa]",
    "Rho"      : "density [kg/m^3]",
    "gamma"    : "plast_strain [ ]"}
    
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
    #passive_tracers ='%s_passive_tracers.pvtu' %fname
    ###############################################
    fname2 = fname+".pvd"
    fname=os.path.join(Folder,Test_Name,fname2)
    time, Flist, n_tim_step =_file_list(fname)    # Retrive the file list and the time step information
    ##############################################
    ts0 = Flist[0]
    fn=ts0[1:ts0.find('/')]
    Filename_0 = [os.path.join(Folder,Test_Name,fn,dyn),os.path.join(Folder,Test_Name,fn,phase)]
    folder_ = os.path.join(Folder,Test_Name)
    C = Coordinate_System(Filename_0,ptsave,(-600.0,600.0),(-600.0,50.0))
    Phase_DB = Phase_Data_Base(folder_)
    Initial_Condition = Initial_condition(Phase_DB.Phase_6_,Phase_DB.Phase_5_,vIC)
    Initial_Condition.tc = Initial_Condition.tc/3.5
    ###############################################
    FSurf = FS(C,len(time)) # Create the instance of free surface class 
    DYN   = VAL(C,dic_val)  # Create the instance of the .pvd file 
    Ph    = Phase(C,phase_dictionary) # Create the instance of the Phase field
    Slab  = SLAB(C,len(time))        # Create the instance of the Slab. 
    ################################################
    ###############################################################################
    # Read Information related to the initial condition
    ###############################################################################
    ipic= 0 
    for istp in Flist:
    ########################### Files and time stepping information############
        t1_start = perf_counter()
        fn=istp[1:istp.find('/')]
        t_cur=time[ipic]        
        Filename_dyn=os.path.join(Folder,Test_Name,fn,dyn)
        Filename_ph=os.path.join(Folder,Test_Name,fn,phase)
        Filename_s=os.path.join(Folder,Test_Name,fn,surf)
    #    Filename_ptr=os.path.join(Folder,Test_Name,fn,passive_tracers)
        
        ###### Retrieve the field that is needed ##################################
        t1 = perf_counter()
        DYN._update_Val(Filename_dyn,C,Initial_Condition)
        t2 = perf_counter()
        Values_time = t2-t1
        print("Dynamic value took","{:02}".format(Values_time))
        ###########################################################################
        t1 = perf_counter()
        Ph._update_phase(Filename_ph,C)  
        t2 = perf_counter()    
        average_plot_time = t2-t1
        print("Phase_update","{:02}".format(average_plot_time))

        ###########################################################################
        t1= perf_counter()
        FSurf._Update_(Filename_s,C,ipic)
        FSurf.ASCI_FILE(ipic,t_cur,Test_Name,ptsave,C.x)
        t2 = perf_counter()
        print("Free surface ","{:02}".format(t2-t1))

        ###########################################################################
        t1 = perf_counter()
        Slab. _update_C(C,DYN,ipic,time)
        t2 = perf_counter()
        print("Slab routine ","{:02}".format(t2-t1))
        ###########################################################################
        t1 = perf_counter()   
        if (ipic % 10 == 0):
            # Plot the field of interest of dynamic 
            DYN._plot_maps_V(t_cur,C.z,C.x,ptsave,ipic)
            # Plot the phase plot 
            Ph._plot_phase_field(C, DYN, ptsave, ipic, t_cur,FSurf,t_cur/Initial_Condition.td)
            # Plot the Slab averages 
            Slab._plot_average_C(t_cur,C.z,ptsave,ipic)          
            t2 = perf_counter()
        plot_map_time = t2-t1
        ###########################################################################
        #  Free Surface 
        #######################################################################
          
        FS_time = t2-t1
        t1_end = perf_counter()
        tstep_time_pr = t1_end-t1_start
        minutes, seconds = divmod(t1_end-t1_start, 60)
          
        ipic +=1 
        
        
        
        
        t_FINE = perf_counter()      
        
        minutes, seconds = divmod(t_FINE-t_INIZIO, 60)
        print("===========================================")
        
        print("||Script took ","{:02}".format(int(minutes)),"minutes and","{:05.2f}".format(seconds),' seconds ||' )
        
        print("===========================================")
    
    ind_z_t = Slab._plot_slab_time(time/Initial_Condition.td,C.z,ptsave,Test_Name,ptsave,30,Initial_Condition)
    _plot_D_D0(Slab,Initial_Condition,ptsave,time,Test_Name,10)


