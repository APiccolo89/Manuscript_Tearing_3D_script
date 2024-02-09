#%%

# Import all the variable

import sys,os,fnmatch
import numpy as np
import matplotlib.pyplot as plt
#plt.switch_backend('agg')
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
import pandas as pd
from Class_Data_Base import * 
from time import perf_counter 



path = r'C:\Users\Andrea Piccolo\Dropbox\Bayreuth_Marcel\Conference_Abstract_Poster_Presentation\Manuscript_2\Picture_Folder_Raw\output_4\Data_Base\Data_base_Slab_detachment_3D.hdf5'

path_save = r'C:\Users\Andrea Piccolo\Dropbox\Bayreuth_Marcel\Conference_Abstract_Poster_Presentation\Manuscript_2\Picture_Folder_Raw\Data_Base_pics'
if not os.path.isdir(path_save):
    os.mkdir(path_save)

# [A.2] Create the database
DB = Data_Base(path)


# Loop over the avaiable test and collect data. 
DB._post_process_data(path,path_save)



# %%
