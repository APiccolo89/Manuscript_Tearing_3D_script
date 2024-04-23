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




# %%

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

#_merge_database(file_A,file_B,file_C,file_D,file_DB)


path_save = r'../../Data_Base_KIT_GLA'
initial_geotherm(path_save)

if not os.path.isdir(path_save):
    os.mkdir(path_save)
# Merge Data Base


# [A.2] Create the database
DB = Data_Base(file_DB)


# Loop over the avaiable test and collect data, save a smaller database, and ascii file {TIP: if you want 
# to create a subclasses you do not have to put a decorator for the classes}
#DB._post_process_data(path,path_save,False,False)
path_figure = os.path.join(path_save,'figure_manuscript')


TSD2 = Test(file_DB,['PR_r','TSD2'])
#TSD2.print_topography_timeseries(path_save,0.0)
#TSD2_V10 = Test(file_DB,['PR_r','TSD2_V10'])
#TSD2_V10.print_topography_timeseries(path_save,0.0)
#TSD2_V11 = Test(file_DB,['PR_r','TSD2_V11'])
#TSD2_V11.print_topography_timeseries(path_save,0.0)
#TSD2_V12 = Test(file_DB,['PR_r','TSD2_V12'])
#TSD2_V12.print_topography_timeseries(path_save,0.0)
#TSD2_V13 = Test(file_DB,['PR_r','TSD2_V13'])
#TSD2_V13.print_topography_timeseries(path_save,0.0)
#TSD2_V15 = Test(file_DB,['PR_r','TSD2_V15'])
#TSD2_V15.print_topography_timeseries(path_save,0.0)
TSD2_V15 = Test(file_DB,['PR_r','TSD2_V15'])

if not os.path.isdir(path_figure):
    os.mkdir(path_figure)

make_plot_dD(path_save,TSD2,TSD2_V15)
# Figure 3: 
"""
Figure 3: Figure depicting two end member scenarios (TSD2 and TSD2_V15)
The figure represent the evolution of D^{\dagger} along x_s (distance from left
to right along the trench) for different timestep. In the column a) we have the fast
tearing model, in the column b) we have the slow tearing model. 

Observation: [Fast tearing ] Fast tearing is characterized by the development of an instability 
within the slab at ca 200 km 

"""
make_figure_3(TSD2,path_figure,'figure3a',[0.1,1.0],[0.1,1.0])

make_figure_3(TSD2_V15,path_figure,'figure3g',[0.1,.2],[0.0,0.1])

"""
Figure 4: 
Figure depicting the total uplift within the tearing processes in the upper row 
Figure depicting the total uplift within 1->end of tearing 

"""
make_figure_4(TSD2,TSD2_V15,path_figure,'figure4',True)
make_figure_4(TSD2,TSD2_V15,path_figure,'figure4b',False)


"""
Figure 5: 

"""
# %%
make_figure5(DB,path_figure,'figure5')

make_figure6(DB,path_figure,'figure6')

make_figure7(DB,path_figure,'figure7')

make_figure8(DB,path_figure,'figure8')





# %%
