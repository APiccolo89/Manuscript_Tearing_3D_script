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
import scipy.io as sio  
import h5py   
"""
Import useful function 
"""
# Read the function to read the output of LaMEM
from Read_VTK_files_LAMEM import * #Read_VTK_files_LAMEM
from Read_VTK_files_LAMEM import  _file_list
from Parser_File import *
# Read the main file for processing the simulation
from Running_script_main_function import _run_script_visualization
from Slab_detachment import *




# Folder where the tests are contained
Folder = r'/scratch/bt307806/3D_Numerical_suites' 
# Folder where the output must be saved
ptsave1 = r'/scratch/bt307806/3D_Numerical_suites/output_'

#Folder = r'/bgi/bt307806/Marcel_project/Lukas_Project/Average_Temperature_Tests' 
#ptsave1 = r'/bgi/bt307806/Marcel_project/Lukas_Project/Output'
# Check if ptsave1 exists, in case create
if not os.path.isdir(ptsave1):
    os.mkdir(ptsave1)
# List of tests {Could be as well a list automatically created}
Test_Name =['T2_AV0_v3_v10_10000']#,'T1_AV0_v10_1000','T1_AV0_v3_10000','T1_AV0_v5_10000','T1_AV0_v10_10000']
#l_path =[ os.path.join(Folder,Test_Name)]
# In case the test are contained in a server
List_Folder_Server = os.listdir(Folder)
# Loop of test
for L in Test_Name:
    # Which test is going to be post process
    print(L)
    # Where the test is located
    l_ = os.path.join(Folder,L)
    #ptsave = os.path.join(ptsave,L)
    #Folder_DB = ptsave 
    # Is this folder existing, in case create the output folder 
    if os.path.isdir(l_):
        ptsave = os.path.join(ptsave1,L)
        if not os.path.isdir(ptsave):
            os.mkdir(ptsave)
        try:
            ptsave_DB = os.path.join(ptsave1,'Data_Base')
            if not os.path.isdir(ptsave_DB):
                os.mkdir(ptsave_DB)
            _run_script_visualization(ptsave,Folder,L,l_,ptsave_DB,'Obliquity')
        except:
            print(L, 'has problem to be post processed, file corrupted, or simply empty folder go futher')
    else: 
        print(L)

