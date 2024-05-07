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
from figure_functions import *
from figure_functions import make_figure_3
from figure_functions import make_figure_4
from figure_functions import make_figure5
from figure_functions import _make_gif
from figure_functions import make_figure_3N




path = r'../../Data_Bases'
name_PR = 'Data_base_Slab_detachment_3D_PR_r.hdf5'
name_PR200 = 'Data_base_Slab_detachment_3D_PR_200.hdf5'
name_PR600 ='Data_base_Slab_detachment_3D_PR_600.hdf5'
name_PRNO ='Data_base_Slab_detachment_3D_PR_no.hdf5'

name_DB    = '3D_numerical_suites.hd5f'
file_A = os.path.join(path,name_PR)
file_B = os.path.join(path,name_PR200)
file_C = os.path.join(path,name_PR600)
file_D = os.path.join(path,name_PRNO)
file_DB = os.path.join(path,name_DB)


# Merge Data Base
#_merge_database(file_A,file_B,file_C,file_D,file_DB)


path_save = r'../../Data_Base_KIT_GLA'

if not os.path.isdir(path_save):
    os.mkdir(path_save)


# [A.2] Create the database
DB = Data_Base(file_DB)


# Loop over the avaiable test and collect data, save a smaller database, and ascii file {TIP: if you want 
# to create a subclasses you do not have to put a decorator for the classes}
DB._post_process_data(path,path_save,False,False,True)
path_figure = os.path.join(path_save,'figure_manuscript')
path_gif = os.path.join(path_figure,'gif')

if not os.path.isdir(path_figure):
    os.mkdir(path_figure)


if not os.path.isdir(path_gif):
    os.mkdir(path_gif)


TSD2 = Test(file_DB,['PR_r','TSD2'])

TSD2_V15 = Test(file_DB,['PR_r','TSD2_V15'])

initial_topography(TSD2,path_save)

initial_geotherm(path_save)


# Figure 3: 
"""
Figure 3: Figure depicting two end member scenarios (TSD2 and TSD2_V15)
The figure represent the evolution of D^{\dagger} along x_s (distance from left
to right along the trench) for different timestep. In the column a) we have the fast
tearing model, in the column b) we have the slow tearing model. 

Observation: [Fast tearing ] Fast tearing is characterized by the development of an instability 
within the slab at ca 200 km 

"""

make_figure_3N(TSD2,path_figure,'figure3a',[10.15,10.20,10.24])

make_figure_3N(TSD2_V15,path_figure,'figure3b',[15.41,25.09,29.69])

"""
Figure 4: 
Figure depicting the total uplift within the tearing processes in the upper row 
Figure depicting the total uplift within 1->end of tearing 

"""
make_figure_4(TSD2,TSD2_V15,path_figure,'figure4a',1)

make_figure_4(TSD2,TSD2_V15,path_figure,'figure4b',2)

make_figure_4(TSD2,TSD2_V15,path_figure,'figure4c',3)

"""
Figure 5: 

"""
make_figure5(DB,path_figure,'figure5')

"""
Figure 6: 

"""

make_figure6(DB,path_figure,'figure6')








# %%
