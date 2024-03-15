
###
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
            if (test_name[1] != 'TSD2_HC_3') & (test_name[1] != 'TSD2_PR'):
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

                print(mean_v)

                # Collect starting time of detachment 

                self.Starting_tearing[itest] = np.nanmin(det_vec[(x_sp>=100) & (x_sp<=1100) ])

                self.Ending_tearing[itest]  = np.nanmax(det_vec[(x_sp>=100) & (x_sp<=1100) ])

                # Plot the difference of topography between starting and ending detachment 

                i1 = np.where(test.time==self.Starting_tearing[itest])[0][0]
                i2 = np.where(test.time==self.Ending_tearing[itest])[0][0]

                # Here I do not have a clear measure for assessing the impact of the tearing. 
                # The uplift associated to the tearing seems to starts before the actual geometrical tearing
                # So, during the first iteration I would like to plot the following picture: 
                # [A] Total dH during the geometrical tearing
                # [B] Total dH from 0.1->End of the tearing 
                # [C] Total dH with different time between 0.1->Beginning of tearing 
                # =============================================
                uplift_G,uplift_2,uplift_T,dH_A,dH_B,dH_C,dtA,dtB,dtC=_compute_dH_tearing(i1,i2,test.FS,test.C,test.Det,test.time)

                self.uplift[itest,0]= uplift_G
                self.uplift[itest,1]= uplift_2
                self.uplift[itest,2]= uplift_T

                self.dt[itest,0]    = dtA
                self.dt[itest,1]    = dtB 
                self.dt[itest,2]    = dtC 
                ipic = 0 
                ASCI_time_Vec(test.time,test_name,path_save_b)
                topo_t = np.zeros([len(test.C.x_sp),len(test.time)],dtype=float)
                dH_t = np.zeros([len(test.C.x_sp),len(test.time)],dtype=float)
                Lithos_FB = np.zeros([len(test.C.x_sp),len(test.time)],dtype=float)
                dF = np.zeros([len(test.C.x_sp),len(test.time)],dtype=float)

                

                path_save_c = os.path.join(path_save_b,'FreeSurface')
                if not os.path.isdir(path_save_c):
                    os.mkdir(path_save_c)                

                # Compute the initial y-coordinate of the slab
                x_trench, y_trench = _compute_initial_slab_position(test.C); 
                for i in range(len(test.time)-1):
                    ASCI_FILE_ALT(test.FS,ipic,test.time[ipic],test_name,path_save_b,test.C)
                    # Plot figure
                    # 1. Find topography at the trench {interpolate topography and uplift}
                    # 2. Create the figure
                    # ============ Post Process data ============================
                    # 1. -> Per each point belonging to the trench 
                    # -> compute anomaly of uplift with time 
                        # a. Compute the amplitude 
                        # b. Compute the wave length 
                        # c. plot 2D maps with time 
                        # d. collect data along profile for Paul 
                        # e. create table and figure 
                    topo_t[:,ipic],dH_t[:,ipic],Lithos_FB[:,ipic],x_trench,y_trench = find_topography_uplift_trench(test.FS,test.C,ipic,x_trench,y_trench)
                    if i > 0:
                        dF[:,ipic] = np.abs((Lithos_FB[:,ipic]-Lithos_FB[:,ipic-1])/((test.time[ipic]-test.time[ipic-1])*1e6))
                    _plot_detachment_topography(ipic,test.time[ipic],os.path.join(path_save_b,'uplift_detachment'),test.Det,dH_t,test.C.x_sp,r'$\dot{H}, [mm \cdot yr^{-1}]$',test.C,dF,r'$\dot{F^{lit}_z}, [N/m/yr]$')
                    _plot_detachment_topography(ipic,test.time[ipic],os.path.join(path_save_b,'Topo_detachment'),test.Det,topo_t,test.C.x_sp,r'${H}, [km]$',test.C,Lithos_FB,r'${F^{lit}_z}$ [10^{12} N/m]')
                    label = Label('$x$, [$km$]','$y$, [$km$]','$none$',r'$H$,[$km$]','no','Topography',[5,95])
                    _plot_2D_surface(test.time[i],test.FS,it,test.C,path_save_c,'H','cmc.oleron',label,ipic,x_trench,y_trench)
                    label = Label('$x$, [$km$]','$y$, [$km$]','$none$',r'$\bar{\dot{\varepsilon}}_{II}$,[$\frac{1}{s}$]','yes','Average lithospheric strain rate',[30,90])
                    _plot_2D_surface(test.time[i],test.FS,it,test.C,path_save_c,'eps','cmc.devon',label,ipic,x_trench,y_trench)
                    _scatter_dH_dF(dF[:,ipic],dH_t[:,ipic],ipic,test.time[ipic],os.path.join(path_save_b,'scatterdF'))
                    # find anomaly # 
                    find_anomaly_wave_length(test.C,test.FS,ipic, path_save_c)
                    

                
                    ipic+=1 
    
                
                path_save_c = os.path.join(path_save_b,'FreeSurface')
                


                
                label = Label('$x$, [$km$]','$y$, [$km$]','$Uplift$',r'$\bar{dH} [m]$','yes','Uplift',[30,90])
                
                scale = 1000 #km to meter
                _plot_Uplift([test.time[i1], test.time[i2]],dH_A*scale,test_name[1],test.C,path_save_c,'Uplift','nipy_spectral',label,'Geometric')
                _plot_Uplift([test.time[i1]-2, test.time[i2]],dH_B*scale,test_name[1],test.C,path_save_c,'Uplift','nipy_spectral',label,'2Myr')
                _plot_Uplift([0.5, test.time[i2]],dH_C*scale,test_name[1],test.C,path_save_c,'Uplift','nipy_spectral',label,'Beginning')
                
                scale = (1000*100)/1e6
                
                label = Label('$x$, [$km$]','$y$, [$km$]','$Uplift rate$',r'${\dot{dH}} [\frac{mm}{yrs}]$','yes','Uplift',[30,90])
                _plot_Uplift([test.time[i1], test.time[i2]],((dH_A)/dtA)*scale,test_name[1],test.C,path_save_c,'Uplift_r','cmc.lapaz',label,'Geometric')
                _plot_Uplift([test.time[i1]-2, test.time[i2]],((dH_B)/dtB)*scale,test_name[1],test.C,path_save_c,'Uplift_r','cmc.lapaz',label,'2Myr')
                _plot_Uplift([0.5, test.time[i2]],((dH_C/dtC)*scale),test_name[1],test.C,path_save_c,'Uplift_r','cmc.lapaz',label,'Beginning')

                # Delete the variable and start all over again. 
            itest = itest+1

        
        print(itest)
        # Print the global data base pictures 
        # Plot average velocity of tearing w.r.t. activation volume, 
        # Diamonds Pr = 200, circle = 400, square = 600 MPa stress limiter
        # Average temperature colorbar
        markers = ["d","o","s","P"]
        label_s = label_scatter(r'$V^{a}_{dis}$ $[10^6 \frac{m^3}{Pa}]$',
                                r'$\bar{v}_{tearing}$ $[\frac{cm}{yr}]$',
                                r'Average tearing velocity',
                                markers
                                ,'yes',
                                'cmc.bilbao',
                                r'$\bar{T}^{S}$ $[^{\circ}C]$')
        fields = ['Avolume','detachment_velocity','Temp','StressLimit']
        _scatter_plot_(self,path_save,label_s,fields,'Average_velocity_global_dataset')
        markers = ["d","o","s","P"]
        label_s = label_scatter(r'$\bar{v}_{tearing}$ $[\frac{cm}{yr}]$',
                                r'$\bar{\frac{dH}{dt}},$mm/yrs$',
                                r'Average tearing velocity',
                                markers
                                ,'yes',
                                'cmc.bilbao',
                                r'$\bar{T}^{S}$ $[^{\circ}C]$')
        fields = ['detachment_velocity','Uplift_rate','Temp','StressLimit']
        _scatter_plot_(self,path_save,label_s,fields,'Average_velocity_global_dataset')
        






@timer
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
                     'y_1': ['/Slab_Detachment/x1','km','Y coordinate center slab, phase'],
                     'y_2': ['/Slab_Detachment/x2','km','Y coordinate center slab, phase']
        }
        
        self.xg = DB._read_variable(self.dict['xg'],Test_name)
        
        self.yg = DB._read_variable(self.dict['yg'],Test_name)
        
        self.zg = DB._read_variable(self.dict['zg'],Test_name)
        
        self.xp = DB._read_variable(self.dict['xp'],Test_name)
        
        self.yp = DB._read_variable(self.dict['yp'],Test_name)
        
        self.zp = DB._read_variable(self.dict['zp'],Test_name)
        
        self.x_sp = DB._read_variable(self.dict['x_sp'],Test_name)
        
        self.y_sp = DB._read_variable(self.dict['y_sp'],Test_name)
        
        self.y_1 = DB._read_variable(self.dict['y_1'],Test_name)
        
        self.y_2 = DB._read_variable(self.dict['y_2'],Test_name)

        
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
                     'x_slab1'  : ['/Slab_Detachment/x1', 'km','Position slab x1'] , 
                     'x_slab2'  :  ['/Slab_Detachment/x2', 'km','Position slab x2'],
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
        self.x1 = DB._read_variable(self.dict['x_slab1'],Test_name)
        self.x2 = DB._read_variable(self.dict['x_slab2'],Test_name)


# Class containing the information of the free surface
class FS():
    def __init__(self,DB:Data_Base,Test_name:str):
        self.dict = {'H': ['/FS/Amplitude','km', 'Amplitude'],
                     'dH': ['/FS/dH','mm/yr', 'Rate of variation of Amplitude with time'],
                     'vz_M': ['/FS/vz_M','mm/yr', 'Filtered v_z of free surface'],
                     'Thickness': ['/FS/thickness','km', 'Thickness of the lithosphere '],
                     'tau_mean': ['/FS/mean_stress','MPa', 'Mean stress'],
                     'Topo': ['/FS/Topo','km', 'Topography'],
                     'eps' : ['/FS/mean_eps','1/s','strain rate']
       }
        self.H = DB._read_variable(self.dict['H'],Test_name)
        self.dH = DB._read_variable(self.dict['dH'],Test_name)
        self.vz_M = DB._read_variable(self.dict['vz_M'],Test_name)
        self.Thickness = DB._read_variable(self.dict['Thickness'],Test_name)
        self.tau_mean = DB._read_variable(self.dict['tau_mean'],Test_name)
        self.Topo = DB._read_variable(self.dict['Topo'],Test_name)
        self.eps = DB._read_variable(self.dict['eps'],Test_name)

 
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


class Label():
    def __init__(self,xlabel,ylabel,zlabel, cbar_label, log,title,lim): 
        self.xlabel = xlabel
        self.ylabel = ylabel 
        self.zlabel = zlabel
        self.cbar_label = cbar_label 
        self.log = log
        self.title = title
        self.min = lim[0]
        self.max = lim[1] 
        
def _merge_database(FileA:str,FileB:str,FileC:str,Dest_file:str):
        import h5py 

        with h5py.File(Dest_file,'w') as f_dest:
            with h5py.File(FileA,'r') as f_src:
                f_src.copy(f_src["/PR_r"],f_dest["/"],"PR_r")
                with h5py.File(FileB,'r') as f_src2:
                    f_src.copy(f_src2["/PR_200"],f_dest["/"],"PR_200")
                    with h5py.File(FileC,'r') as f_src3:
                        f_src3.copy(f_src3["/PR_600"],f_dest["/"],"PR_600")
                        
                        
def _compute_dH_tearing(i1:int,i2:int,Surf:FS,C:C,D:Det,time:float):
    
    i1_t_2Ma = np.max(np.where(time<time[i1]-2.0))
    i1_0_5  = np.min(np.where(time>0.5))
    # Compute anomaly
    dH_A = Surf.Topo[:,:,i2]-Surf.Topo[:,:,i1]
    dH_B = Surf.Topo[:,:,i2]-Surf.Topo[:,:,i1_t_2Ma]
    dH_C = Surf.Topo[:,:,i2]-Surf.Topo[:,:,i1_0_5]
    # 
    iA   = np.abs(dH_A)>np.mean(np.abs(dH_A))
    iB   = np.abs(dH_B)>np.mean(np.abs(dH_B))
    iC   = np.abs(dH_C)>np.mean(np.abs(dH_C))
    uplift_G = np.nanmean(dH_A[iA==1])
    uplift_2 = np.nanmean(dH_B[iB==1])
    uplift_T = np.nanmean(dH_C[iC==1])
    dtA      = time[i2]-time[i1]
    dtB      = time[i2]-time[i1_t_2Ma]
    dtC      = time[i2]-time[i1_0_5]
    
    return uplift_G,uplift_2,uplift_T,dH_A,dH_B,dH_C,dtA,dtB,dtC
    



class MidpointNormalize(mpl.colors.Normalize):
    def __init__(self, vmin, vmax, midpoint=0, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        normalized_min = max(0, 1 / 2 * (1 - abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax))))
        normalized_max = min(1, 1 / 2 * (1 + abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin))))
        normalized_mid = 0.5
        x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
        return np.ma.masked_array(np.interp(value, x, y))



class label_scatter():
    def __init__(self,xlabel:str,ylabel:str,title:str,markers:str,log:str,colormap:str,cbar_label:str):
        self.xlabel           = xlabel
        self.ylabel           = ylabel
        self.title            = title 
        self.markers          = markers 
        self.log              =  log 
        self.colormap         =  colormap
        self.cbar_label       = cbar_label 

def _plot_2D_surface(time:float,Data,Test_name,C:C,path_save:str,field:str,colorbar:str,label:Label,ipic,xs,ys):
    import cmcrameri as cmc 

    
    buf = eval(field,globals(),Data.__dict__)
    if label.log == 'yes':
        buf = np.log10(buf)
    if isinstance(Data,Det):
        buf[buf==-np.inf]= np.nan
    # Find the most reliable limits for the colorbar
    ptsave_c = os.path.join(path_save,field)
    if not os.path.isdir(ptsave_c):
            os.mkdir(ptsave_c)
    min = np.nanpercentile(buf,label.min)
    max = np.nanpercentile(buf,label.max)
    print('color maps limits are %2f  and %2f' %(min,max))
    if isinstance(Data,FS):
        val = np.zeros((len(C.yg),len(C.xg)),dtype=float)
        x = C.xg 
        y = C.yg
    elif isinstance(Data,Det):
        val = np.zeros((len(C.x_sp),len(C.zp)),dtype=float)
        x = C.x_sp 
        y = C.zp 
    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(15*cm, 15*cm))
    ax0 = fg.gca()
    if (min/abs(min)!= max/abs(max)):
        norm = MidpointNormalize(vmin=min, vmax=max, midpoint=0.0)
        cf =ax0.pcolormesh(x, y, val,norm=norm ,shading='gouraud')
    else: 
        cf =ax0.pcolormesh(x, y, val,shading='gouraud')
    cf1 = ax0.plot(xs,ys,linewidth=1.5,color='r')
    cbar = fg.colorbar(cf, ax=ax0,orientation='horizontal',extend="both",label=label.cbar_label)
    val = buf[:,:,ipic]
    tick = r"Time =  %s [$Myr$]" %("{:.3f}".format(time))
    fna='Fig'+"{:03d}".format(ipic)+'.png'
    fn = os.path.join(ptsave_c,fna)
   
    cf.set_array(val.ravel())
    cf.set_cmap(colorbar)

    cf.set_clim([min,max])
    cbar.vmin = min 
    cbar.vmax = max
    cbar.update_normal(cf)
    ax0.tick_params(axis='both', which='major', labelsize=14)
    ax0.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
    ax0.set_ylim(np.min(y),np.max(y))
    ax0.set_xlim(np.min(x),np.max(x))
    cbar.set_label(label.cbar_label)
    plt.title(tick,fontsize=15)
    fg.patch.set_facecolor('white')
    
    #plt.show()
        
    fg.savefig(fn,dpi=300)

#@timer          
def  _scatter_plot_(Data:Data_Base,path_save:str,label_scatter:label_scatter,fields:list,name_figure):
    import cmcrameri as cmc 

    
    x_f,y_f,z_f,m_f = fields # unpack the field 
    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(15*cm, 15*cm))
    ax0 = fg.gca()
    x =  eval(x_f,globals(),Data.__dict__)
    if y_f == 'Uplift_rate':
        du =  eval('uplift',globals(),Data.__dict__)*1000*1000
        dt = eval('dt',globals(),Data.__dict__)*1e6
        y  = du[:,0]/dt[:,0]
    else: 
        y = eval(y_f,globals(),Data.__dict__)
   # if label_scatter.log == 'yes':
        #y = np.log10(y)
   
    z =  eval(z_f,globals(),Data.__dict__)
    m =  eval(m_f,globals(),Data.__dict__)
    m_u = np.unique(m) 
    for im in range(len(m_u)):
        plt.scatter(x[m==m_u[im]],y[m==m_u[im]],60,z[m==m_u[im]],marker=label_scatter.markers[im],cmap = label_scatter.colormap,edgecolors='k')
        
    cbar = fg.colorbar( ax=ax0,orientation='horizontal',extend="both",label=label_scatter.cbar_label)
    ax0.set_clim([820,np.max(z)])
    cbar.vmin = 820 
    cbar.vmax = np.max(z)
    ax0.tick_params(axis='both', which='major', labelsize=14)
    ax0.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
    cbar.set_label(label_scatter.cbar_label)
    plt.title(label_scatter.title,fontsize=15)
    ax0.set_xlabel(label_scatter.xlabel, fontsize=14)
    ax0.set_ylabel(label_scatter.ylabel, fontsize=14)
    plt.show() 
    fg.patch.set_facecolor('white')
    name_fig = '%s.png' %(name_figure)
    fn = os.path.join(path_save,name_fig)  
    #plt.show()
            
    fg.savefig(fn,dpi=300)
    
    
def _plot_Uplift(time_v:float,dH,Test_name,C:C,path_save:str,field:str,colorbar:str,label:Label,type:str):
    import cmcrameri as cmc 
    

    buf = dH 
    buf[np.abs(dH)<np.mean(np.abs(dH))] = np.nan
    # Find the most reliable limits for the colorbar
    ptsave_c = os.path.join(path_save,field)
    if not os.path.isdir(ptsave_c):
            os.mkdir(ptsave_c)
    min = np.nanpercentile(buf,label.min)
    max = np.nanpercentile(buf,label.max)
    print('color maps limits are %2f  and %2f' %(min,max))
    x = C.xg 
    y = C.yg
    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(15*cm, 15*cm))
    ax0 = fg.gca()
        
    tick = r"Average uplift from %s to %s [Myrs]" %("{:.3f}".format(time_v[0]),"{:.3f}".format(time_v[1]))
    fna='Fig'+type+'.png'
    fn = os.path.join(ptsave_c,fna)
   
    cf =ax0.contourf(x, y, buf,)
    cbar = fg.colorbar(cf, ax=ax0,orientation='horizontal',extend="both",label=label.cbar_label)
    cf.set_cmap(colorbar)
    cbar.update_normal(cf)
    ax0.tick_params(axis='both', which='major', labelsize=14)
    ax0.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
    ax0.set_ylim(np.min(y),np.max(y))
    ax0.set_xlim(np.min(x),np.max(x))
    cbar.set_label(label.cbar_label)
    plt.title(tick,fontsize=15)
    fg.patch.set_facecolor('white')
    plt.draw()
    plt.show()
        
    fg.savefig(fn,dpi=300)
    
def ASCI_FILE_ALT(S,ipic,t_cur,Test_Name,ptsave,C:C):
            
    """
    Write a simple ascii file for the post processing of the free surface dat
    This is for the the free surface data, later on I will dedicate a bit of 
    more time on the usage of the passive tracers.     
    """
    file_name = str(ipic).zfill(7)+'__'+Test_Name[1]+'Free_surface_data.txt'
    
    ptsave_b=os.path.join(ptsave,'DataBase_FS')
    if not os.path.isdir(ptsave_b):
        os.mkdir(ptsave_b)
    
    filename = os.path.join(ptsave_b,file_name)
    Y,X = np.meshgrid(C.xg,C.yg)
    buf_x = X.ravel()
    buf_y = Y.ravel()
    vz_M    = S.vz_M[:,:,ipic]
    dH    = S.dH[:,:,ipic]
    H     = S.H[:,:,ipic]
    S        = np.array([buf_x*1e3,buf_y*1e3,vz_M.ravel(),dH.ravel()*1000,H.ravel()*1000])
    if(os.path.isfile(filename)):
        os.remove(filename)
    f = open(filename, 'a+')
    f.write('########################################\n')
    f.write('time [Myrs] time step []\n')
    f.write('x, y,v_z,dHdt, Topography\n')
    f.write('  [m],[m],[mm/yrs],[mm/yrs], [m]\n')
    f.write('########################################\n')
    f.write('time = %6f, timestep = %d\n' %(t_cur,ipic))
    f.write('\n')
    np.savetxt(f, np.transpose(S),fmt='%.6f', delimiter=' ', newline = '\n') 
    f.close()
    #print('Free surface data of the timestep %d, has been printed' %(ipic))
    
def ASCI_time_Vec(time,Test_Name,ptsave):
            
    """
    Write a simple ascii file for the post processing of the free surface dat
    This is for the the free surface data, later on I will dedicate a bit of 
    more time on the usage of the passive tracers.     
    """
    file_name = 'Time_Vector'+'__'+Test_Name[1]+'.txt'
    
    ptsave_b=os.path.join(ptsave,'DataBase_FS')
    if not os.path.isdir(ptsave_b):
        os.mkdir(ptsave_b)
    
    filename = os.path.join(ptsave_b,file_name)
    dt = np.diff(time)
    dt_s = 0.0*time
    dt_s[1:]=dt[:]
    S        = np.array([time,dt_s])

    if(os.path.isfile(filename)):
        os.remove(filename)
    f = open(filename, 'a+')
    f.write('########################################\n')
    f.write('Time_Vector\n')
    f.write('time dt\n')
    f.write('  [Myrs],[Myrs]\n')
    f.write('########################################\n')
    f.write('\n')
    np.savetxt(f, np.transpose(S),fmt='%.6f', delimiter=' ', newline = '\n') 
    f.close()
    
# Auxilary function


def _plot_detachment_topography(ipic,time_sim,ptsave_b,D:Det,field:float,x_s:float,field_name,C:C,FB,label_2):
    fna='Fig'+str(ipic)+'.png'
    fg = figure()
    tick=r'$Time = %s Myrs$' %(time_sim)
    if not os.path.isdir(ptsave_b):
        os.mkdir(ptsave_b)
    fn = os.path.join(ptsave_b,fna)
    ax1 = fg.add_axes([0.1, 0.7, 0.8, 0.2])
    ax0 = fg.add_axes([0.1, 0.05, 0.8, 0.5])
#    for ip in range(20):
#        it = ipic - ip
#        alpha_v= 0.8-(ip+1)*(1/29)
#        if ip == 0: 
#            cc = 'r'
#        else:
#            cc = 'b'
#        if (it == 0) & (ip == 0) :
#            ax1.plot(x_s,field[:,it],c = cc,alpha = alpha_v,linewidth=alpha_v)
#            break
#        if (ip >0 ) & (it == 0 ):
#            ax1.plot(x_s,field[:,it],c = cc,alpha = alpha_v,linewidth=alpha_v)
#            break 
#        else:
#            ax1.plot(x_s, field[:,it],c = cc,alpha = alpha_v,linewidth=alpha_v)

    ax1.plot(x_s,field[:,ipic],color='r',linewidth=1.2)
    ax1.set_ylabel(field_name)
    ax1.set_xlabel(r'$x_s, [km]$')
 
 
    
    #ax1.set_title(tick)
    ax1.set_xlim(0, 1200)           
    ax1.set_yscale('linear')    
    ax3= ax1.twinx()
    ax3.plot(x_s,FB[:,ipic]/1e12,color='k',linewidth=1.2) 
    ax3.set_ylabel(label_2)
    ax3.set_xlabel(r'$x_s, [km]$')
 
    
    buf = D.D[:,:,ipic]/100 # Hard coded value i know. 
    levels = np.linspace(np.round(0.1), np.round(0.85), num=10, endpoint=True, retstep=False, dtype=float)
    cf=ax0.pcolormesh(x_s,C.zp,np.transpose(buf),cmap='inferno',vmin = 0.1, vmax=0.85)
    cbar = fg.colorbar(cf,ax=ax0,orientation='horizontal')
    ax1.set_title(tick)
    ax0.tick_params(axis='both', which='major', labelsize=5)
    ax0.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
    ax0.set_ylabel(r'$z, [km]$')
    ax0.set_xlabel(r'$x_s, [Myrs]$')
    fg.tight_layout()    
    plt.draw()    # necessary to render figure before saving
    fg.savefig(fn,dpi=600)
    plt.close()

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


def _compute_initial_slab_position(C:C):
    
    y1 = C.y_1[:,:,0]
    y2 = C.y_2[:,:,0]
    y1[y1 == -np.inf] = np.nan 
    y2[y2== -np.inf] = np.nan 
    y1_mean = np.nanmean(y1,1)
    y2_mean = np.nanmean(y2,1)
    y_trench = (y1_mean+y2_mean)/2
    x_trench  = C.xp[(C.xp>=-600)& (C.xp<=600)]
    

    return x_trench, y_trench


def _scatter_dH_dF(dF,dH,ipic,time_sim,ptsave_b):
    
    fna='Fig'+str(ipic)+'.png'
    fg = figure()
    fn = os.path.join(ptsave_b,fna)

    tick=r'$Time = %s Myrs$' %(time_sim)
    if not os.path.isdir(ptsave_b):
        os.mkdir(ptsave_b)
    ax0 = fg.gca()
    ax0.set_title(tick)

    cf=ax0.scatter(dF,dH,s=10,color='#a4c330',edgecolors='k')
    ax0.tick_params(axis='both', which='major', labelsize=5)
    ax0.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
    ax0.set_xlabel(r'$\dot{F}, [N/m/yr]$')
    ax0.set_ylabel(r'$\dot{H}, [mm/yr]$')
    fg.tight_layout()    
    plt.draw()    # necessary to render figure before saving
    fg.savefig(fn,dpi=600)
    ax0.plot()

"""
find_anomaly_dH
function 

"""
def find_anomaly_wave_length(C:C, F:FS, ipic:int,path_save_c):
