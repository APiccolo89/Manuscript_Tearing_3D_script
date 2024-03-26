
###
import sys,os,fnmatch
import numpy as np
import matplotlib.pyplot as plt
#plt.switch_backend('agg')
import matplotlib
from matplotlib import cm
import scipy as scp 
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
import time
from functools import wraps
import matplotlib as mpl
from Auxilary_function import *







def timer(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print(f"Execution time of {func.__name__}: {round(end - start,2)} seconds")
        return result
    return wrapper

# Main class where to collect all the data related to a specific test. 
@timer
class Data_Base():
    """
    [1] -> Initial draft of the class: instead of collecting all the data and storing, creating dynamic variable field
    i store the path to the variable and generate the variable whenever I need. 
    [2] -> function that reads the tests and collects their name of the test name
    [3] -> dictionary that connect path to variable that needs to be customizable. 
    [4] -> function that creates the numpy array for plotting, post processing and collecting data
    {Panda array to handle tables of data}
    """
    def __init__(self,path:str): 
        
        self.path = path
        
        self.Tests_Name, self.n_test= self._read_test(path)
                
        self.detachment_velocity = np.zeros([self.n_test-1],dtype = float)
        
        self.Starting_tearing   = np.zeros([self.n_test-1],dtype = float)
        
        self.Ending_tearing    = np.zeros([self.n_test-1],dtype = float)
        
        self.uplift            = np.zeros([self.n_test-1,3],dtype=float) 
        
        self.dt             = np.zeros([self.n_test-1,3],dtype=float) 
        
        self.Avolume = np.zeros([self.n_test-1],dtype = float)
        
        self.StressLimit = np.zeros([self.n_test-1],dtype = float)
        
        self.Temp      = np.zeros([self.n_test-1],dtype = float)
        
        
    def _read_test(self,path:str):
        
        # Read File Data base 
        # Open File
        
        f = h5py.File(path, 'r')
        # Collect the test names
        
        # List group 
        
        LG = list(f["/"])
        List = []
        n_test = 0 
        for ig in LG:
            lt = list(f["/"+ig].keys())
            for it in lt:
                List.append([ig ,it])
            n_test += len(list(f["/"+ig].keys()))
            
            
        # Close file 
        
        f.close()
        # Return the data needed
        
        return List, n_test
    
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
        
        path2variable = "/%s/%s%s" %(Test_name[0],Test_name[1],keys[0])
        # Create array 
        buf = np.array(f[path2variable])
        
        f.close()
        return buf
    def _post_process_data(self,path:str,path_save:str):
        itest = 0 
        for it in range(self.n_test-1):
            test_name = self.Tests_Name[it]
            if (test_name[1] != 'TSD2_HC_3'):
                path_save_b = os.path.join(path_save,test_name[1])
                if not os.path.isdir(path_save_b):
                    os.mkdir(path_save_b)

                print(test_name)

                test = Test(self,test_name)
                # Collect initial data
                self.Avolume[itest] = test.IC.VnM*1e6

                self.StressLimit[itest] = test.IC.tau_co

                self.Temp[itest] = test.IC.T_av

                # Compute the average velocity of the slab tearing
                # Plot Slab surface

                x_sp = test.C.x_sp 

                v_test = test.Det.vel_tear[0,:]

                det_vec = test.Det.det_vec

                det_vec[test.Det.depth_vec>=-60.0]=np.nan

                corrected_depth = test.Det.depth_vec

                corrected_depth[test.Det.depth_vec>=-60.0]=np.nan

                mean_v = (1000.0*1000*100)/(np.nanmax(det_vec[(x_sp>=100) & (x_sp<=1100) ])*1e6-np.nanmin(det_vec[(x_sp>=100) & (x_sp<=1100) ])*1e6)

                self.detachment_velocity[itest] = mean_v

                # Collect starting time of detachment 

                self.Starting_tearing[itest] = np.nanmin(det_vec[(x_sp>=100) & (x_sp<=1100) ])

                self.Ending_tearing[itest]  = np.nanmax(det_vec[(x_sp>=100) & (x_sp<=1100) ])
                
                # Here I do not have a clear measure for assessing the impact of the tearing. 
                # The uplift associated to the tearing seems to starts before the actual geometrical tearing
                # So, during the first iteration I would like to plot the following picture: 
                # [A] Total dH during the geometrical tearing
                # [B] Total dH from 0.1->End of the tearing 
                # [C] Total dH with different time between 0.1->Beginning of tearing 
                # =============================================
                
                
                self.uplift[itest,0]= np.nanmean(test.FS.Uplift_det[test.FS.Uplift_det>
                                                                    np.nanmean(test.FS.Uplift_det)])
                self.uplift[itest,1]= np.nanmean(test.FS.Uplift_LT[test.FS.Uplift_det>
                                                                    np.nanmean(test.FS.Uplift_det)])
                self.uplift[itest,2] = self.uplift[itest,0]/self.uplift[itest,1]
                
                print_det_prof(test.C.x_sp,test.Det.D_x_t_det_F,path_save_b,'detachment')
                print_det_prof(test.C.x_sp,test.Det.tau_x_t_det_F,path_save_b,'tau')


                ipic = 0 
                ASCI_time_Vec(test.time,test_name,path_save_b)


            itest = itest+1

        print(itest)




@timer
class Test():
    def __init__(self,DB:Data_Base,Test_name:str):

        self.time = DB._read_variable(['/time','Myr', 'Time vector'],Test_name)
        
        self.time_M = (self.time[0:1:-1]+self.time[1:1:])/2

        self.IC = IC(DB,Test_name,self)

        self.C  = C(DB, Test_name,self)
                
        self.Det = Det(DB,Test_name,self)
        
        self.FS = FS(DB, Test_name,self)
        
        self.Ptr = Ptr(DB, Test_name,self)
            
        



# Class containing the coordinate system information 
class C():
    def __init__(self,DB:Data_Base,Test_name:str,T:Test):
        
        self.dict = {'xg': ['/Coordinate_System/x','km', 'Numerical Grid x'],
                     'yg': ['/Coordinate_System/y','km', 'Numerical Grid y'],
                     'zg': ['/Coordinate_System/y','km', 'Numerical Grid z'],
                     'xp': ['/Coordinate_System/xp','km', 'Refined Grid x'],
                     'yp': ['/Coordinate_System/yp','km', 'Refined Grid y'],
                     'zp': ['/Coordinate_System/zp','km', 'Refined Grid z'],
                     'x_sp': ['/Slab_Detachment/x_s','km','X coordinate trench, REFINED'],
                     'y_sp': ['/Slab_Detachment/y_b', 'km','Y coordinate trench, REFINED'],
                     'y_1': ['/Slab_Detachment/x1','km','Y coordinate center slab, phase'],
                     'y_2': ['/Slab_Detachment/x2','km','Y coordinate center slab, phase'],
                     'ind_x_trench_g':['n.a.','logical','index trench computational grid'],
                     'x_trench_p': ['n.a.','km','trench position along x direction in refined grid'],
                     'y_trench_p': ['n.a.','km','trench position along y direction in refined grid']
                     
        }
        
        self.xg = DB._read_variable(self.dict['xg'],Test_name)
        
        self.yg = DB._read_variable(self.dict['yg'],Test_name)
        
        self.zg = DB._read_variable(self.dict['zg'],Test_name)
        
        self.xp = DB._read_variable(self.dict['xp'],Test_name)
        
        self.yp = DB._read_variable(self.dict['yp'],Test_name)
        
        self.zp = DB._read_variable(self.dict['zp'],Test_name)
        
        self.x_sp = DB._read_variable(self.dict['x_sp'],Test_name)
        
        self.y_sp = DB._read_variable(self.dict['y_sp'],Test_name)
        
        # Compute the initial position of the slab 
        
        y_1 = DB._read_variable(self.dict['y_1'],Test_name)
        
        y_2 = DB._read_variable(self.dict['y_2'],Test_name)

        self.ind_x_trench_g,self.y_trench_p,self.x_trench_p=self._compute_initial_slab_position(y_1,y_2,T)
        
    def _compute_initial_slab_position(self,y_1,y_2,T:Test):
        """
        Function that retrieve the initial position of the slab and compute the reference 
        mid surface of the slab in x-y plane. 
        input: 
        y_1,y_2 => position top and bottom surface,
        T: Test data
        self: the FS sub class
        output: 
        y_trench_p: y position of the midsurface in the refined grid [float]
        x_trench_p: x position of the midsurface in the refined grid [float]
        ind_x_trench_g : indexes of the trench along x direction in the
                        computational grid [boolean array]

        """
    
        y1 = y_1[:,:,0]
        y2 = y_2[:,:,0]
        y1[y1 == -np.inf] = np.nan 
        y2[y2== -np.inf] = np.nan 
        y1_mean = np.nanmean(y1,1)
        y2_mean = np.nanmean(y2,1)
        y_trench_p = (y1_mean+y2_mean)/2
        x_trench_p  = self.xp[(self.xp>=np.min(T.IC.coordinate_Slab))& (self.xp<=
                                                                np.max(T.IC.coordinate_Slab))]
        ind_x_trench_g = (self.xg>=np.min(T.IC.coordinate_Slab))& (self.xg<=
                                     np.max(T.IC.coordinate_Slab))

        return ind_x_trench_g, y_trench_p,x_trench_p

        
# Class containing the Initial condition information  
class IC():
    def __init__(self,DB:Data_Base,Test_name:str,T:Test):
        self.dict = {'L0': ['/IC/L0','km', 'Length of Slab'],
                     'D0': ['/IC/D0','km', 'Thickness of Slab'],
                     'T_av': ['/IC/T_av','C', 'Average Temperature at -100 km'],
                     'Coord_Slab': ['/IG/Slab/B_main_coordinate','km','Data of the slab boundary type'],
                     'etarefS': ['/IC/eta_ref_S','Pas', 'Effective viscosity slab at reference condition'],
                     'etarefM': ['/IC/eta_ref_UM','Pas', 'Average effective viscosity of the mantle at tau0'],
                     'xiS': ['/IC/xiUS','n.d.', 'Dominance dislocation of the slab'],
                     'xiM': ['/IC/xiUM','n.d.','Dominance dislocation of the mantle'],
                     'tau0': ['/IC/tau0', 'Pa','Reference stress'],
                     'VnM' : ['/Phase_DB/Phase_5_/Rheology/Dislocation/V','m3/J','Activation volume mantle'],
                     'VnS' : ['/Phase_DB/Phase_6_/Rheology/Dislocation/V','m3/J','Activation volume Slab'],
                     'tau_co':['/Phase_DB/Phase_6_/Rheology/Plastic/ch','Pa','Stress limiter slab'],
        }
        self.L0 = DB._read_variable(self.dict['L0'],Test_name)
        self.D0 = DB._read_variable(self.dict['D0'],Test_name)
        self.T_av = DB._read_variable(self.dict['T_av'],Test_name)
        self.etarefS = DB._read_variable(self.dict['etarefS'],Test_name)
        self.etarefM = DB._read_variable(self.dict['etarefM'],Test_name)
        self.xiS = DB._read_variable(self.dict['xiS'],Test_name)
        self.xiM = DB._read_variable(self.dict['xiM'],Test_name)
        self.tau0 = DB._read_variable(self.dict['tau0'],Test_name)
        self.VnM = DB._read_variable(self.dict['VnM'],Test_name)
        self.VnS = DB._read_variable(self.dict['VnS'],Test_name)
        self.tau_co = DB._read_variable(self.dict['tau_co'],Test_name)
        coordinate = DB._read_variable(self.dict['Coord_Slab'],Test_name)
        self.coordinate_Slab = coordinate[0:1] 


# Class containing the information related to the detachment
class Det():
    def __init__(self,DB:Data_Base,Test_name:str,T:Test):
        self.dict = {'D': ['/Slab_Detachment/D','km', 'Thickness of the slab with time (xs-z)'],
                     'Psi': ['/Slab_Detachment/Psi','W/m3', 'Dissipative rate energy production'],
                     'T': ['/Slab_Detachment/T','C', 'Average Temperature of the slab with time (xs-z)'],
                     'tau_max': ['/Slab_Detachment/tau_max','MPa', 'Average Stress of the slab with time (xs-z)'],
                     'depth_vec': ['/Slab_Detachment/depth_vec','km', 'Depth of detachment '],
                     'det_vec': ['/Slab_Detachment/det_vec','Myr', 'Time of detachment '],
                     'tau_vec': ['/Slab_Detachment/tau_vec','MPa.', 'Stress at the detachment'],
                     'x_vec': ['/Slab_Detachment/x_vec','km','x position of detachment'],
                     'y_vec': ['/Slab_Detachment/y_vec', 'km','y position of detachment'],
                     'vel_tear' : ['/Slab_Detachment/average_tearing_velocity', 'cm/yrs','Velocity of detachment'],
                     'x_slab1'  : ['/Slab_Detachment/x1', 'km','Position slab x1'] , 
                     'x_slab2'  :  ['/Slab_Detachment/x2', 'km','Position slab x2'],
                     }
        self.D = DB._read_variable(self.dict['D'],Test_name)
        self.Psi = DB._read_variable(self.dict['Psi'],Test_name)
        self.T = DB._read_variable(self.dict['T'],Test_name)
        self.tau_max = DB._read_variable(self.dict['tau_max'],Test_name)
        self.depth_vec = DB._read_variable(self.dict['depth_vec'],Test_name)
        self.det_vec = DB._read_variable(self.dict['det_vec'],Test_name)
        self.tau_vec = DB._read_variable(self.dict['tau_vec'],Test_name)
        self.x_vec = DB._read_variable(self.dict['x_vec'],Test_name)
        self.y_vec = DB._read_variable(self.dict['y_vec'],Test_name)
        self.vel_tear = DB._read_variable(self.dict['vel_tear'],Test_name)
        self.y1 = DB._read_variable(self.dict['x_slab1'],Test_name)
        self.y2 = DB._read_variable(self.dict['x_slab2'],Test_name)
        # Derivative values 
        self.D_x_t_det = np.zeros([len(self.x_vec),len(T.time)],dtype = float)
        self.tau_x_t_det = np.zeros([len(self.x_vec),len(T.time)],dtype = float)
        # Filtered {these data are horrible to see without a moving average filter}
        self.D_x_t_det_F = np.zeros([len(self.x_vec),len(T.time)],dtype = float)
        self.tau_x_t_det_F = np.zeros([len(self.x_vec),len(T.time)],dtype = float)
        i_along_x        = np.zeros(len(self.x_vec),dtype = int) 
        #Find_index 
        for i in range(len(i_along_x)):
            ind = np.where(T.C.zp == self.depth_vec[i])
            if len(ind[0])>0:
                i_along_x[i]=ind[0][0]
            else:
                i_along_x[i]=-1 
        
        # find time evolution thickness, stress along the depth at which the detachment is occuring 
        
        self.time_evolution_necking(i_along_x,len(T.time))
    
    def time_evolution_necking(self,i_along_x,itime):
        """
        Function that simply select the nodes of the array that corresponds to the depth of detachment and saves the entire 
        timeseries.
        input:
        i_along_x = the depth index at which detachment occurs.
        itime = number of timestep
        output:
        self 

        """

        for i in range(len(i_along_x)):
            if i_along_x[i] != -1:
                self.D_x_t_det[i,:] = self.D[i,i_along_x[i],:]
                self.tau_x_t_det[i,:] = self.tau_max[i,i_along_x[i],:]
            else:
                self.D_x_t_det[i,:] = -np.inf
                self.tau_x_t_det[i,:] = -np.inf

        # Beautyfing the array. 
        for i in range(itime):
            self.D_x_t_det_F[:,i] = np.convolve(self.D_x_t_det[:,i], np.ones(30)/30, mode='same')
            self.tau_x_t_det_F[:,i] = np.convolve(self.tau_x_t_det[:,i], np.ones(30)/30, mode='same')

        
        
        return self 


        

        




# Class containing the information of the free surface
class FS():
    def __init__(self,DB:Data_Base,Test_name:str,T:Test):
        self.dict = {'H': ['/FS/Amplitude','km', 'Amplitude'],
                     'dH': ['/FS/dH','mm/yr', 'Rate of variation of Amplitude with time'],
                     'vz_M': ['/FS/vz_M','mm/yr', 'Filtered v_z of free surface'],
                     'Thickness': ['/FS/thickness','km', 'Thickness of the lithosphere '],
                     'tau_mean': ['/FS/mean_stress','MPa', 'Mean stress'],
                     'Topo': ['/FS/Topo','km', 'Topography'],
                     'eps' : ['/FS/mean_eps','1/s','strain rate'],
                     'vz' : ['/FS/vz','mm/yr','strain rate'],
       }
        self.H = DB._read_variable(self.dict['H'],Test_name)
        self.dH = DB._read_variable(self.dict['dH'],Test_name)
        self.vz_M = DB._read_variable(self.dict['vz_M'],Test_name)
        self.vz = DB._read_variable(self.dict['vz'],Test_name)
        self.Thickness = DB._read_variable(self.dict['Thickness'],Test_name)
        self.tau_mean = DB._read_variable(self.dict['tau_mean'],Test_name)
        self.Topo = DB._read_variable(self.dict['Topo'],Test_name)
        self.eps = DB._read_variable(self.dict['eps'],Test_name)
        # Derivative Data Set [filtered and post processed]
        self.dH_fil =  np.zeros(np.shape(self.dH),dtype=float)
        self.vz_fil =  np.zeros(np.shape(self.dH),dtype=float)
       
        self.dH_long_term = np.zeros(np.shape(self.dH[:,:,0]),dtype = float)
        self.dH_detachment = np.zeros(np.shape(self.dH[:,:,0]),dtype = float)
       
        self.Uplift_LT = np.zeros(np.shape(self.dH[:,:,0]),dtype=float)
        self.Uplift_det = np.zeros(np.shape(self.dH[:,:,0]),dtype=float)

        
        self.dH_fil = self.filter_array('dH')
        self.vz_fil = self.filter_array('vz')
        
        
        
        
    
    def filter_array(self,key):
        """
        Internal function that creates a filtered array
        Input argument: 
        key = the array that needs to be filtered. 
        self = the class 
        """
        # access to the variable 
        buf_array = eval(key,globals(),self.__dict__)
        # create new variable 
        new_array = np.zeros(np.shape(buf_array),dtype=float)
        # Prepare loop
        nit = len(new_array[0,0,:])
        # Loop 
        for i in range(nit):
            new_array[:,:,i]=scp.ndimage.median_filter(buf_array[:,:,i], size=(7))
    
        return new_array 

    def _compute_dH_tearing(self,T:Test):
        """
        Input: FS and T:test 
        ===============================
        Output: updated FS
        
        Function that selects:
        1) the topography at the beginning of the tearing 
        and its end, and create an array containing
        the total variation of topography
        2) the topography after the isostatic rebound,and the end of tearging
        and creating an array containing the total variation of topography
        during the long term history of the model
        Function that compute the associated average uplift. 
        """
        # Find when do the tearing starts: 
        i1 = np.nanmin(T.Det.det_vec[(T.C.x_sp>=100) & (T.C.x_sp<=1100) ])
        # Find when do the tearing end: 
        i2 = np.nanmax(T.Det.det_vec[(T.C.x_sp>=100) & (T.C.x_sp<=1100) ])
        # Filter isostatic 
        i_iso = np.where(T.time>1.0)
        i_iso = i_iso[0][0]

        dt_det = T.time[i2[0][0]]-T.time[i1[0][0]]
        dt_long = T.time[i2[0][0]]-T.time[i_iso]
        # Compute anomaly
        self.dH_detachment= self.Topo[:,:,i2]-self.Topo[:,:,i1]
        self.dH_long_term = self.Topo[:,:,i2]-self.Topo[:,:,i_iso]
       
        self.Uplift_det = self.dH_detachment/dt_det
        self.Uplift_LT  = self.dH_long_term/dt_long
        
        return self
    

# Class containing the Passive tracers information   
class Ptr(): 
    def __init__(self,DB:Data_Base,Test_name:str,T:Test):
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

        
def _merge_database(FileA:str,FileB:str,FileC:str,Dest_file:str):
        import h5py 

        with h5py.File(Dest_file,'w') as f_dest:
            with h5py.File(FileA,'r') as f_src:
                f_src.copy(f_src["/PR_r"],f_dest["/"],"PR_r")
                with h5py.File(FileB,'r') as f_src2:
                    f_src.copy(f_src2["/PR_200"],f_dest["/"],"PR_200")
                    with h5py.File(FileC,'r') as f_src3:
                        f_src3.copy(f_src3["/PR_600"],f_dest["/"],"PR_600")
                        
                        



#@timer          


# Auxilary function



def find_topography_uplift_trench(F:FS,C:C,ipic:int,x_trench,y_trench):
    Topography = F.H[:,:,ipic]
    dH         = F.dH[:,:,ipic]
    # Interpolate topography to x_trench y_trench 
    topo_trench = _interpolate_2D(Topography,C.xg,C.yg,x_trench,y_trench)
    dH_trench = _interpolate_2D(dH,C.xg,C.yg,x_trench,y_trench)
    tau_trench = _interpolate_2D(F.tau_mean[:,:,ipic],C.xg,C.yg,x_trench,y_trench)
    thickness = _interpolate_2D(F.Thickness[:,:,ipic],C.xg,C.yg,x_trench,y_trench)
    FB_trench = tau_trench*1e6*thickness*1e3 

    return topo_trench,dH_trench,FB_trench,x_trench,y_trench





def _interpolate_2D(Topo:float,xg:float,yg:float,xp:float,yp:float):
    topo_marker = np.zeros([len(xp)],dtype = float)
    for i in range(len(xp)):
        xx = xp[i]
        yy = yp[i]
        ix = find1Dnodes(xg,xx,len(xg))
        iy = find1Dnodes(yg,yy,len(yg))
        x1 = xg[ix]
        x2 = xg[ix+1]
        y1 = yg[iy]
        y2 = yg[iy+1]
        intp1=Topo[iy,ix]
        intp2=Topo[iy,ix+1]
        intp3=Topo[iy+1,ix+1]
        intp4=Topo[iy+1,ix]
        val_ = bilinearinterpolation(xx,yy,x1,x2,y1,y2,intp1,intp2,intp3,intp4)
        topo_marker[i]=val_
    return topo_marker

def find1Dnodes(cord,cordm,number):
    # I was fucking up a few stuff:
    #buf = cordm-cord 
    #min = np.min(np.abs(buf))
    for index in range(number):
        if (cord[index]>cordm):
            break 
    return index-1    

def bilinearinterpolation(xx,yy,x1,x2,y1,y2,intp1,intp2,intp3,intp4):
    wx=(xx-x1)/(x2-x1)
    wy=(yy-y1)/(y2-y1)

    # FUCKING BUG: intp4 -> 1-x*intp4
    R=intp1*(1-wx)*(1-wy)+intp2*wx*(1-wy)+intp3*wx*wy+intp4*wy*(1-wx)

    return R    

def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} km" if plt.rcParams["text.usetex"] else f"{s} km"

