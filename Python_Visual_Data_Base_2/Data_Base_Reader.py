
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
import argparse 
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy.ma import masked_array
import matplotlib.cbook as cbook
from Class_Data_Base import * 
from Class_Data_Base import _merge_database
from time import perf_counter 
import figure_functions as ff 

parser = argparse.ArgumentParser()

parser.add_argument("path", help="path to the databases folder",type=str)
parser.add_argument("path_save", help="path to the choosen save folder (i.e., where the user wants to save the results)",type=str)
parser.add_argument("save_h5", help = 'Save a small database with relevant data')
parser.add_argument("print_txt", help = 'Print .txt file suitable for other software such Petrell')
parser.add_argument("merge_database",help = 'Merge database from the first step of postprocessing')

args = parser.parse_args()

path           = args.path 
path_save      = args.path_save
print(path_save)
save_h5        = args.save_h5
print_txt      = args.print_txt
merge_database = args.merge_database 







#path = r'../../Data_Bases' ;path_save = r'../../Data_Base_KIT_GLA'
#merge_database = False; save_h5 = False; print_txt = False 

name_PR    = 'Data_base_Slab_detachment_3D_PR_r.hdf5'
name_PR200 = 'Data_base_Slab_detachment_3D_PR_200.hdf5'
name_PR600 = 'Data_base_Slab_detachment_3D_PR_600.hdf5'

name_DB    = '3D_numerical_suites.hd5f'
file_A     = os.path.join(path,name_PR)
file_B     = os.path.join(path,name_PR200)
file_C     = os.path.join(path,name_PR600)
file_DB    = os.path.join(path,name_DB)

# Merge Data Base
if merge_database == True:
    _merge_database(file_A,file_B,file_C,file_DB)


if not os.path.isdir(path_save):
    os.mkdir(path_save)
    print(path_save)


# [A.2] Create the database
DB = Data_Base(file_DB)

# Loop over the avaiable test and collect data, save a smaller database, and ascii file {TIP: if you want 
# to create a subclasses you do not have to put a decorator for the classes}
DB._post_process_data(path,      # path to DB
                      path_save, # path to save the postprocessed data
                      save_h5,     # save database of postprocessed data
                      print_txt,     # print topography in .txt {compatible with Petrel}
                      True)      # Check if exist already a post processed dataset within the current 
                                 # h5 file; if not, repeat the postprocess 

path_figure = os.path.join(path_save,'FM_def')

path_gif = os.path.join(path_figure,'gif')

if not os.path.isdir(path_figure):
    os.mkdir(path_figure)


if not os.path.isdir(path_gif):
    os.mkdir(path_gif)

ff.figure9(DB,path_figure,'figure9')

TSD2 = Test(file_DB,['PR_r','TSD2'])

TSD2_V15 = Test(file_DB,['PR_r','TSD2_V15'])

TSD2_V10 = Test(file_DB,['PR_r','TSD2_V10'])

TSD2_V11 = Test(file_DB,['PR_r','TSD2_V11'])

TSD2_V12 = Test(file_DB,['PR_r','TSD2_V12'])

TSD2_V13 = Test(file_DB,['PR_r','TSD2_V13'])


# Figure 1 
ff.figure1_par(path_figure)
# Figure 2: 
"""
Figure 3_5: Figure depicting two end member scenarios (TSD2 and TSD2_V15)
The figure represent the evolution of D^{\dagger} along x_s (distance from left
to right along the trench) for different timestep. In the column a) we have the fast
tearing model, in the column b) we have the slow tearing model. 

Observation: [Fast tearing ] Fast tearing is characterized by the development of an instability 
within the slab at ca 200 km 

"""
# Figure3 
ff.figure3_5(TSD2,path_figure,'figure3',[10.15,10.20,10.24],False)
# Figure 4 PARAVIEW

# Figure 5
ff.figure3_5(TSD2_V15,path_figure,'figure5',[15.41,25.09,29.69],False)
"""
Figure 6: 
Figure depicting the total uplift within the tearing processes in the upper row 
Figure depicting the total uplift within 1->end of tearing 
"""
ff.figure6(TSD2,TSD2_V15,path_figure,'figure6',1)
"""
Figure 7
"""
ff.figure_S4_S7(TSD2,    path_figure,'Figure_7A',[])
ff.figure_S4_S7(TSD2_V15,path_figure,'Figure_7B',['[e]','[f]','[g]','[h]'])
"""
Figure 8: Figure depicting the tearing velocity as a function of the main controlling parameter 
"""
ff.figure8(DB,path_figure,'figure8')
"""
Figure 9: correlation between uplift and tearing velocity
"""
ff.figure9(DB,path_figure,'figure9')
"""
Figure 10: migration velocity of slab
"""
ff.figure_10(TSD2,TSD2_V15,[9.1,10.20,11.70],[5,25.09,38.00],path_figure,'figure10')
"""
Supplementary figure: Figure that are in the supplementary+gif_maker
"""

# Figure S2
ff.figure_S2(path_figure)
# Figure S3
ff.figure_S3(TSD2,path_figure)


ff.figure_S4_S7(TSD2_V10,path_figure,'Figure_S4',[])
ff.figure_S4_S7(TSD2_V11,path_figure,'Figure_S5',[])
ff.figure_S4_S7(TSD2_V12,path_figure,'Figure_S6',[])
ff.figure_S4_S7(TSD2_V13,path_figure,'Figure_S7',[])
ff.figure_S8(DB,path_figure,'Figure_S8')
ff.figure_S9(TSD2,path_figure,'Figure_S9',[10.18,10.20,10.24])
ff.figure_S10(DB,path_figure,'Figure_S10')
ff.figure_S11(TSD2,TSD2_V10,TSD2_V11,TSD2_V12,TSD2_V13,TSD2_V15,path_figure,'Figure_S11')
ff.figure_S12(TSD2,TSD2_V10,TSD2_V11,TSD2_V12,TSD2_V13,TSD2_V15,path_save,'Figure_S12')

# Make Gif
ff._make_gif(TSD2,path_gif)
ff._make_gif(TSD2_V10,path_gif)
ff._make_gif(TSD2_V11,path_gif)
ff._make_gif(TSD2_V12,path_gif)
ff._make_gif(TSD2_V13,path_gif)
ff._make_gif(TSD2_V15,path_gif)


