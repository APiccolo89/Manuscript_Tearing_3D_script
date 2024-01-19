
###
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
import pandas as pd

# Main class where to collect all the data related to a specific test. 
class Data_Base():
    def __init__(self,SDet,name_group,LM,correction): 

        self._create_field(SDet,name_group,LM,correction)
        
    def _create_field(self,SDet,name_group,LM,correction):
        
        import h5py 
        ########################################
        #  Open the h5 file 
        ########################################
        f = h5py.File(SDet, 'r')
        SG = []
        frate = []
        Sub_Group= list(f["/"+name_group].keys())
        Tests = []
        for isub in Sub_Group:
            SG.append(isub)
            test = list(f["/"+name_group+"/"+isub].keys())
            print("%s is processed and test in SG %s are currently integrated" %(name_group,isub))
            ifal = 0
            tot_test_SG = len(test)
            for itest in test:
                t_path = "/%s/%s/%s" %(name_group,isub,itest)
                print(itest)
                e = t_path+"/Initial_Data/" in f
                if e == True:
                    locals()[itest] = Test_Data()
                    self.Tests.append(itest)


                    ###############################################################
                    # Fill up the data base of the local test
                    ###############################################################
                    if(np.array(f["/%s/%s/%s/failed" %(name_group,isub,itest)])==1):
                        ifal +=1
                        # Place Holder Read initial setup 
                        locals()[itest].D0 = np.array(f[t_path+"/Initial_Data/Initial_Condition/D0"])
                        locals()[itest].L0 = np.array(f[t_path+"/Initial_Data/Initial_Condition/L0"])
                        locals()[itest].epsc = np.array(f[t_path+"/Initial_Data/Initial_Condition/epsc"])
                        locals()[itest].eta0 = np.array(f[t_path+"/Initial_Data/Initial_Condition/eta0"])
                        locals()[itest].eta0Ast = np.array(f[t_path+"/Initial_Data/Initial_Condition/eta0Ast"])
                        locals()[itest].eta0D = np.array(f[t_path+"/Initial_Data/Initial_Condition/eta0D"])
                        locals()[itest].tau0 = np.array(f[t_path+"/Initial_Data/Initial_Condition/tau0"])
                        locals()[itest].tc  = np.array(f[t_path+"/Initial_Data/Initial_Condition/tc"])
                        locals()[itest].LM  = LM

                        locals()[itest].failed   = 1
                    else: 
                            # Place Holder Read initial setup of the current test 
                    
                        locals()[itest].D0 = np.array(f[t_path+"/Initial_Data/Initial_Condition/D0"])
                        locals()[itest].L0 = np.array(f[t_path+"/Initial_Data/Initial_Condition/L0"])
                        locals()[itest].epsc = np.array(f[t_path+"/Initial_Data/Initial_Condition/epsc"])
                        locals()[itest].eta0 = np.array(f[t_path+"/Initial_Data/Initial_Condition/eta0"])
                        locals()[itest].eta0Ast = np.array(f[t_path+"/Initial_Data/Initial_Condition/eta0Ast"])
                        locals()[itest].eta0D = np.array(f[t_path+"/Initial_Data/Initial_Condition/eta0D"])
                        locals()[itest].tau0 = np.array(f[t_path+"/Initial_Data/Initial_Condition/tau0"])
                        locals()[itest].tc  = np.array(f[t_path+"/Initial_Data/Initial_Condition/tc"])
                        locals()[itest].LM  = LM
                        locals()[itest].LM  = LM
                        locals()[itest].correction  = correction
                        locals()[itest].failed   = 0
                        
                        # Extract Dvec
                                        
                        locals()[itest].D    = np.array(f[t_path+"/TimeEvolution/Slab/D"])
                        locals()[itest].time = np.array(f[t_path+"/TimeEvolution/time"])
                        locals()[itest].eta_S = 10**np.array(f[t_path+"/TimeEvolution/Slab/vis"])
                        locals()[itest].Sdetvec = np.array(f[t_path+"/Detachment/det_vec"]) 
                        locals()[itest].tau = np.array(f[t_path+"/TimeEvolution/Slab/tau"])
                        locals()[itest].eps = np.array(f[t_path+"/TimeEvolution/Slab/eps"])
                        locals()[itest].vslab = np.array(f[t_path+"/TimeEvolution/Slab/v_slab"])

                        # retrive topography 
                        Topo = np.array(f[t_path+"/TimeEvolution/Slab/eps"])
                        """
                        Compute velocity detachment 
                        """
                        locals()[itest].vz = np.array(f[t_path+"/TimeEvolution/Slab/vz"])
                
                    setattr(self, itest,locals()[itest])
                else:
                    ifal+=1 
            frate.append(ifal/tot_test_SG)
            tot_test_SG = [] 
            tests = []
        self.frate = frate 
        self.SG = SG 
        f.close()


class Test_Data():
    def __init__(self):
        self.failed = []
        self.tc = [] 
        self.tau0 = []
        self.D0  = []
        self.etaM = []
        self.L0 = []
        self.eta0 =[] 
        self.epsc = []
        self.Sdetvec= []
        self.D = [] 
        self.tau   = [] 
        self.eps = []
        self.time  = [] 
        self.LM    = []
        self.vz    = []
        self.Topo = []
        self.correction = []
