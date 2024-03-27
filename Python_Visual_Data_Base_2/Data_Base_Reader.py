#%%
# Import all the variable

import sys,os,fnmatch
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib
from matplotlib import cm
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import re
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from pylab import figure, axes, pie, title, show
import pylab
#from matplotlib import rcParams
#import json
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy.ma import masked_array
import matplotlib.cbook as cbook
from Class_Data_Base import * 
from Class_Data_Base import _merge_database
from time import perf_counter 

def make_figure_3(A,B,path_figure,figure_name):
    # figure name
    fn = os.path.join(path_figure,'%s.png'%(figure_name))
    
    # Prepare variables
    b0 = A.Det.D_x_t_det_F
    b1 = B.Det.D_x_t_det_F

    a0 = A.C.x_sp

    i10,i20 = np.where(A.time==np.nanmin(A.Det.det_vec)),np.where(A.time==np.nanmax(A.Det.det_vec))
    i11,i21 = np.where(B.time==np.nanmin(B.Det.det_vec)),np.where(B.time==np.nanmax(B.Det.det_vec))


    # Prepare figure layout 
    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(15*cm, 15*cm))  
    bx = 0.07
    by = 0.1
    sx = 0.45
    dx = 0.03
    sy = 0.8
    dy = []
    # Prepare axis of the two end member 
    ax0 = fg.add_axes([bx, by, sx, sy])
    ax1 = fg.add_axes([bx+sx+dx, by, sx, sy])
    # 

    # Print figure
    ax0.plot(a0,b0[:,i10[0][0]-10:i20[0][0]]/(A.IC.D0/1000),linewidth = 1.0,color = 'k') 
    ax1.plot(a0,b1[:,i11[0][0]-10:i21[0][0]]/(B.IC.D0/1000),linewidth = 1.0,color = 'k')
    ax0.set_ylim(0.05,1.0)
    #ax0.set_ytick(0.1,0.5,0.9)
    ax0.tick_params(axis="y",direction="in")
    ax0.tick_params(axis="x",direction="in")
    ax1.set_ylim(0.05,1.0)
    ax1.tick_params(left=True,right=True,labelleft=False) 
    ax1.tick_params(axis="y",direction="in")
    ax1.tick_params(axis="x",direction="in")

    #ax1.set_ytick(0.1,0.5,0.9)
     

    fg.savefig(fn,dpi=600)






path = r'../../output_def/Data_Base'
name_PR = 'Data_base_Slab_detachment_3D_PR_r.hdf5'
name_PR200 = 'Data_base_Slab_detachment_3D_PR_200.hdf5'
name_PR600 ='Data_base_Slab_detachment_3D_PR_600.hdf5'
name_DB    = '3D_numerical_suites.hd5f'
file_A = os.path.join(path,name_PR)
file_B = os.path.join(path,name_PR200)
file_C = os.path.join(path,name_PR600)
file_DB = os.path.join(path,name_DB)

#_merge_database(file_A,file_B,file_C,file_DB)


path_save = r'../../Data_Base_Pic_file_2'
if not os.path.isdir(path_save):
    os.mkdir(path_save)
# Merge Data Base


# [A.2] Create the database
DB = Data_Base(file_DB)


# Loop over the avaiable test and collect data, save a smaller database, and ascii file {TIP: if you want 
# to create a subclasses you do not have to put a decorator for the classes}
#DB._post_process_data(path,path_save,False,False)

TSD2 = Test(file_DB,['PR_r','TSD2'])

TSD2_V15 = Test(file_DB,['PR_r','TSD2_V15'])

path_figure = os.path.join(path_save,'figure_manuscript')
if not os.path.isdir(path_figure):
    os.mkdir(path_figure)

make_figure_3(TSD2,TSD2_V15,path_figure,'figure3')



# %%
