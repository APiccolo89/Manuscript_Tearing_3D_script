
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
import h5py 

# Main class where to collect all the data related to a specific test. 
class Data_Base():
    """
    [1] -> Initial draft of the class: instead of collecting all the data and storing, creating dynamic variable field
    i store the path to the variable and generate the variable whenever I need. 
    [2] -> function that reads the tests and collects their name of the test name
    [3] -> dictionary that connect path to variable that needs to be customizable. 
    [4] -> function that creates the numpy array for plotting, post processing and collecting data
    {Panda array to handle tables of data}
    """
    def __init__(self,path:str,GroupName:str): 
        
        self.path = path
        
        self.GroupName = GroupName
        
        self.Tests_Name, self.n_test= self._read_test(path,GroupName)
                
        self.detachment_velocity = np.zeros([self.n_test],dtype = float)
        
        self.Starting_tearing   = np.zeros([self.n_test],dtype = float)
        
        self.uplift            = np.zeros([self.n_test,3],dtype=float) 
        
        self.Avolume = np.zeros([self.n_test],dtype = float)
        
        self.StressLimit = np.zeros([self.n_test],dtype = float)
        
        self.Temp      = np.zeros([self.n_test],dtype = float)
        
    def _read_test(self,path:str,GroupName:str):
        
        # Read File Data base 
        # Open File
        
        f = h5py.File(path, 'r')
        # Collect the test names
        
        List= list(f["/"+GroupName].keys())
        # Close file 
        
        f.close()
        # Return the data needed
        
        return List, len(List)
    
    def _read_variable(self,keys:list,Test_name): 
        """_read_variable:
        function that has a list containing the name of the variable and the path within the
        h5 file. 

        Args:
            keys (list): list of data: a. [path,unit measure]
        
        """
        # Open Data base
        
        f = h5py.File(self.path,'r')
        
        # path 2 variable
        
        path2variable = "/%s/%s%s" %(self.GroupName,Test_name,keys[0])
        # Create array 
        buf = np.array(f[path2variable])
        
        f.close()
        return buf
        
        
        
        
        
        
        
class Test():
    def __init__(self,DB:Data_Base,Test_name:str):
        
        self.IC = IC(DB,Test_name)
        
        self.Det = Det(DB,Test_name)
        
        self.FS = FS(DB, Test_name)
        
        self.Ptr = Ptr(DB, Test_name)
        
        self.C  = C(DB, Test_name)
    
        self.time = DB._read_variable(['/time','Myr', 'Time vector'],Test_name)
        
        self.time_M = (self.time[0:1:-1]+self.time[1:1:])/2
        



# Class containing the coordinate system information 
class C():
    def __init__(self,DB:Data_Base,Test_name:str):
        
        self.dict = {'xg': ['/Coordinate_System/x','km', 'Numerical Grid x'],
                     'yg': ['/Coordinate_System/y','km', 'Numerical Grid y'],
                     'zg': ['/Coordinate_System/y','km', 'Numerical Grid z'],
                     'xp': ['/Coordinate_System/xp','km', 'Phase Grid x'],
                     'yp': ['/Coordinate_System/yp','km', 'Phase Grid y'],
                     'zp': ['/Coordinate_System/zp','km', 'Phase Grid z'],
                     'x_sp': ['/Slab_Detachment/x_s','km','X coordinate trench, phase'],
                     'y_sp': ['/Slab_Detachment/y_b', 'km','Y coordinate trench, phase'],
        }
        
        self.xg = DB._read_variable(self.dict['xg'],Test_name)
        
        self.yg = DB._read_variable(self.dict['yg'],Test_name)
        
        self.zg = DB._read_variable(self.dict['zg'],Test_name)
        
        self.xp = DB._read_variable(self.dict['xp'],Test_name)
        
        self.yp = DB._read_variable(self.dict['yp'],Test_name)
        
        self.zp = DB._read_variable(self.dict['zp'],Test_name)
        
        self.x_sp = DB._read_variable(self.dict['x_sp'],Test_name)
        
        self.y_sp = DB._read_variable(self.dict['y_sp'],Test_name)
        
# Class containing the Initial condition information  
class IC():
    def __init__(self,DB:Data_Base,Test_name:str):
        self.dict = {'L0': ['/IC/L0','km', 'Length of Slab'],
                     'D0': ['/IC/D0','km', 'Thickness of Slab'],
                     'T_av': ['/IC/T_av','C', 'Average Temperature at -100 km'],
                     'etarefS': ['/IC/eta_ref_S','Pas', 'Effective viscosity slab at reference condition'],
                     'etarefM': ['/IC/eta_ref_UM','Pas', 'Average effective viscosity of the mantle at tau0'],
                     'xiS': ['/IC/xiUS','n.d.', 'Dominance dislocation of the slab'],
                     'xiM': ['/IC/xiUM','n.d.','Dominance dislocation of the mantle'],
                     'tau0': ['/IC/tau0', 'Pa','Reference stress'],
                     'VnM' : ['/Phase_DB/Phase_5_/Dislocation/V','m3/J','Activation volume mantle'],
                     'VnS' : ['/Phase_DB/Phase_6_/Dislocation/V','m3/J','Activation volume Slab'],
                     'tau_co':['/Phase_DB/Phase_6_/Dislocation/ch','Pa','Stress limiter slab'],
        }
        self.L0 = DB._read_variable(self.dict['L0'],Test_name)
        self.D0 = DB._read_variable(self.dict['D0'],Test_name)
        self.T_av = DB._read_variable(self.dict['T_av'],Test_name)
        self.etarefS = DB._read_variable(self.dict['etarefS'],Test_name)
        self.etarefM = DB._read_variable(self.dict['etarefM'],Test_name)
        self.xiS = DB._read_variable(self.dict['xiS'],Test_name)
        self.xiM = DB._read_variable(self.dict['xiM'],Test_name)
        self.tau0 = DB._read_variable(self.dict['tau0'],Test_name)
        self.VnM = DB._read_variable(self.dict['tau0'],Test_name)
        self.VnS = DB._read_variable(self.dict['tau0'],Test_name)
        self.tau_co = DB._read_variable(self.dict['tau0'],Test_name)


# Class containing the information related to the detachment
class Det():
    def __init__(self,DB:Data_Base,Test_name:str):
        self.dict = {'D': ['/Slab_Detachment/D','km', 'Thickness of the slab with time (xs-z)'],
                     'Psi': ['/Slab_Detachment/Psi','W/m3', 'Dissipative rate energy production'],
                     'T': ['/Slab_Detachment/T','C', 'Average Temperature of the slab with time (xs-z)'],
                     'depth_vec': ['/Slab_Detachment/depth_vec','km', 'Depth of detachment '],
                     'det_vec': ['/Slab_Detachment/det_vec','Myr', 'Time of detachment '],
                     'tau_vec': ['/Slab_Detachment/tau_vec','MPa.', 'Stress at the detachment'],
                     'x_vec': ['/Slab_Detachment/x_vec','km','x position of detachment'],
                     'y_vec': ['/Slab_Detachment/y_vec', 'km','y position of detachment'],
                     'vel_tear' : ['/Slab_Detachment/average_tearing_velocity', 'cm/yrs','Velocity of detachment'],
        }
        self.D = DB._read_variable(self.dict['D'],Test_name)
        self.Psi = DB._read_variable(self.dict['Psi'],Test_name)
        self.T = DB._read_variable(self.dict['T'],Test_name)
        self.depth_vec = DB._read_variable(self.dict['depth_vec'],Test_name)
        self.det_vec = DB._read_variable(self.dict['det_vec'],Test_name)
        self.tau_vec = DB._read_variable(self.dict['tau_vec'],Test_name)
        self.x_vec = DB._read_variable(self.dict['x_vec'],Test_name)
        self.y_vec = DB._read_variable(self.dict['y_vec'],Test_name)
        self.vel_tear = DB._read_variable(self.dict['vel_tear'],Test_name)

# Class containing the information of the free surface
class FS():
    def __init__(self,DB:Data_Base,Test_name:str):
        self.dict = {'H': ['/FS/Amplitude','km', 'Amplitude'],
                     'dH': ['/FS/dH','mm/yr', 'Rate of variation of Amplitude with time'],
                     'vz_M': ['/FS/vz_M','mm/yr', 'Filtered v_z of free surface'],
                     'Thickness': ['/FS/thickness','km', 'Thickness of the lithosphere '],
                     'tau_mean': ['/FS/mean_stress','MPa', 'Mean stress'],
                     'Topo': ['/FS/Topo','km', 'Topography']
       }
        self.H = DB._read_variable(self.dict['H'],Test_name)
        self.dH = DB._read_variable(self.dict['dH'],Test_name)
        self.vz_M = DB._read_variable(self.dict['vz_M'],Test_name)
        self.Thickness = DB._read_variable(self.dict['Thickness'],Test_name)
        self.tau_mean = DB._read_variable(self.dict['tau_mean'],Test_name)
        self.Topo = DB._read_variable(self.dict['Topo'],Test_name)

 
# Class containing the Passive tracers information   
class Ptr(): 
    def __init__(self,DB:Data_Base,Test_name:str):
        self.dict = {'x': ['/PTrBas/x','km', 'x position'],
                     'y': ['/PTrBas/y','km', 'y position'],
                     'z': ['/PTrBas/z','km', 'z position'],
                     'P': ['/PTrBas/P','MPa', 'Pressure passive tracer '],
                     'T': ['/PTrBas/T','C', 'Temperature of passive tracer'],
                    }
        self.x = DB._read_variable(self.dict['x'],Test_name)
        self.y = DB._read_variable(self.dict['y'],Test_name)
        self.z = DB._read_variable(self.dict['z'],Test_name)
        self.P = DB._read_variable(self.dict['P'],Test_name)
        self.T = DB._read_variable(self.dict['T'],Test_name)























# Old code in case I need. 
#class Data_Base():
#    """
#    [1] -> Initial draft of the class: instead of collecting all the data and storing, creating dynamic variable field
#    i store the path to the variable and generate the variable whenever I need. 
#    """
#    def __init__(self,SDet,name_group,LM,correction): 
#
#        self._create_field(SDet,name_group,LM,correction)
#        
#    def _create_field(self,SDet,name_group,LM,correction):
#        
#        import h5py 
#        ########################################
#        #  Open the h5 file 
#        ########################################
#        f = h5py.File(SDet, 'r')
#        SG = []
#        frate = []
#        Sub_Group= list(f["/"+name_group].keys())
#        Tests = []
#        for isub in Sub_Group:
#            SG.append(isub)
#            test = list(f["/"+name_group+"/"+isub].keys())
#            print("%s is processed and test in SG %s are currently integrated" %(name_group,isub))
#            ifal = 0
#            tot_test_SG = len(test)
#            for itest in test:
#                t_path = "/%s/%s/%s" %(name_group,isub,itest)
#                print(itest)
#                e = t_path+"/Initial_Data/" in f
#                if e == True:
#                    locals()[itest] = Test_Data()
#                    self.Tests.append(itest)
#
#
#                    ###############################################################
#                    # Fill up the data base of the local test
#                    ###############################################################
#                    if(np.array(f["/%s/%s/%s/failed" %(name_group,isub,itest)])==1):
#                        ifal +=1
#                        # Place Holder Read initial setup 
#                        locals()[itest].D0 = np.array(f[t_path+"/Initial_Data/Initial_Condition/D0"])
#                        locals()[itest].L0 = np.array(f[t_path+"/Initial_Data/Initial_Condition/L0"])
#                        locals()[itest].epsc = np.array(f[t_path+"/Initial_Data/Initial_Condition/epsc"])
#                        locals()[itest].eta0 = np.array(f[t_path+"/Initial_Data/Initial_Condition/eta0"])
#                        locals()[itest].eta0Ast = np.array(f[t_path+"/Initial_Data/Initial_Condition/eta0Ast"])
#                        locals()[itest].eta0D = np.array(f[t_path+"/Initial_Data/Initial_Condition/eta0D"])
#                        locals()[itest].tau0 = np.array(f[t_path+"/Initial_Data/Initial_Condition/tau0"])
#                        locals()[itest].tc  = np.array(f[t_path+"/Initial_Data/Initial_Condition/tc"])
#                        locals()[itest].LM  = LM
#
#                        locals()[itest].failed   = 1
#                    else: 
#                            # Place Holder Read initial setup of the current test 
#                    
#                        locals()[itest].D0 = np.array(f[t_path+"/Initial_Data/Initial_Condition/D0"])
#                        locals()[itest].L0 = np.array(f[t_path+"/Initial_Data/Initial_Condition/L0"])
#                        locals()[itest].epsc = np.array(f[t_path+"/Initial_Data/Initial_Condition/epsc"])
#                        locals()[itest].eta0 = np.array(f[t_path+"/Initial_Data/Initial_Condition/eta0"])
#                        locals()[itest].eta0Ast = np.array(f[t_path+"/Initial_Data/Initial_Condition/eta0Ast"])
#                        locals()[itest].eta0D = np.array(f[t_path+"/Initial_Data/Initial_Condition/eta0D"])
#                        locals()[itest].tau0 = np.array(f[t_path+"/Initial_Data/Initial_Condition/tau0"])
#                        locals()[itest].tc  = np.array(f[t_path+"/Initial_Data/Initial_Condition/tc"])
#                        locals()[itest].LM  = LM
#                        locals()[itest].LM  = LM
#                        locals()[itest].correction  = correction
#                        locals()[itest].failed   = 0
#                        
#                        # Extract Dvec
#                                        
#                        locals()[itest].D    = np.array(f[t_path+"/TimeEvolution/Slab/D"])
#                        locals()[itest].time = np.array(f[t_path+"/TimeEvolution/time"])
#                        locals()[itest].eta_S = 10**np.array(f[t_path+"/TimeEvolution/Slab/vis"])
#                        locals()[itest].Sdetvec = np.array(f[t_path+"/Detachment/det_vec"]) 
#                        locals()[itest].tau = np.array(f[t_path+"/TimeEvolution/Slab/tau"])
#                        locals()[itest].eps = np.array(f[t_path+"/TimeEvolution/Slab/eps"])
#                        locals()[itest].vslab = np.array(f[t_path+"/TimeEvolution/Slab/v_slab"])
#
#                        # retrive topography 
#                        Topo = np.array(f[t_path+"/TimeEvolution/Slab/eps"])
#                        """
#                        Compute velocity detachment 
#                        """
#                        locals()[itest].vz = np.array(f[t_path+"/TimeEvolution/Slab/vz"])
#                
#                    setattr(self, itest,locals()[itest])
#                else:
#                    ifal+=1 
#            frate.append(ifal/tot_test_SG)
#            tot_test_SG = [] 
#            tests = []
#        self.frate = frate 
#        self.SG = SG 
#        f.close()
#
#
#class Test_Data():
#    def __init__(self):
#        self.failed = []
#        self.tc = [] 
#        self.tau0 = []
#        self.D0  = []
#        self.etaM = []
#        self.L0 = []
#        self.eta0 =[] 
#        self.epsc = []
#        self.Sdetvec= []
#        self.D = [] 
#        self.tau   = [] 
#        self.eps = []
#        self.time  = [] 
#        self.LM    = []
#        self.vz    = []
#        self.Topo = []
#        self.correction = []
#