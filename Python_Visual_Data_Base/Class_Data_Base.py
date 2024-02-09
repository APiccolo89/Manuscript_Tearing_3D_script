
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
            path_save_b = os.path.join(path_save,test_name[1])
            if not os.path.isdir(path_save_b):
                os.mkdir(path_save_b)
            
            print(test_name)
            if test_name[1] == 'TSD2_HC':
                print('Check')
            test = Test(self,test_name)
            # Collect initial data
            self.Avolume[itest] = test.IC.VnM*1e6
            
            self.StressLimit[itest] = test.IC.tau_co
            
            self.Temp[itest] = test.IC.T_av
            
            # Compute the average velocity of the slab tearing
            # Plot Slab surface
            #
            _detect_detachment(test.time,test.C,test.Det.D,test.Det.x1,test.Det.x2,path_save)
            
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

            # Plot figure of current test
            path_save_c = os.path.join(path_save_b,'FreeSurface')
            if not os.path.isdir(path_save_c):
                os.mkdir(path_save_c)
            
            label = Label('$x$, [$km$]','$y$, [$km$]','$none$',r'$\bar{\dot{\varepsilon}}_{II}$,[$\frac{1}{s}$]','yes','Average lithospheric strain rate',[30,90])
            
           # _plot_2D_surface(test.time,test.FS,it,test.C,path_save_c,'eps','cmc.devon',label)
            
            label = Label('$x$, [$km$]','$y$, [$km$]','$none$',r'$H$,[$km$]','no','Topography',[5,95])
            
           # _plot_2D_surface(test.time,test.FS,it,test.C,path_save_c,'H','cmc.oleron',label)

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

@timer
def _plot_2D_surface(time:float,Data,Test_name,C:C,path_save:str,field:str,colorbar:str,label:Label):
    import cmcrameri as cmc 
    
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "sans-serif",
        "font.sans-serif": "Helvetica",
    })
    
    timestep = len(time)
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
    cbar = fg.colorbar(cf, ax=ax0,orientation='horizontal',extend="both",label=label.cbar_label)
    for ipic in range(timestep): 
        val = buf[:,:,ipic]
        tick = r"Time =  %s [$Myr$]" %("{:.3f}".format(time[ipic]))
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
    
    
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "sans-serif",
        "font.sans-serif": "Helvetica",
    })
    
    x_f,y_f,z_f,m_f = fields # unpack the field 
    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(15*cm, 15*cm))
    ax0 = fg.gca()
    x =  eval(x_f,globals(),Data.__dict__)
    y =  eval(y_f,globals(),Data.__dict__)
   # if label_scatter.log == 'yes':
        #y = np.log10(y)
   
    z =  eval(z_f,globals(),Data.__dict__)
    m =  eval(m_f,globals(),Data.__dict__)
    m_u = np.unique(m) 
    s = ax0.scatter(x,y,60,z,marker=label_scatter.markers[0],cmap = label_scatter.colormap,edgecolors='k')
    cbar = fg.colorbar(s, ax=ax0,orientation='horizontal',extend="both",label=label_scatter.cbar_label)
    s.set_clim([820,np.max(z)])
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

def _detect_detachment(time:float,C:C,D:float,x1:float,x2:float,ptsave:str):
    
    itl = len(time)
    x_sp = C.x_sp
    z_p  = C.zp
    
    for i in range(itl)-1:
        D_ = D[:,:,i]
        x1_ = x1[:,:,i]
        x2_ = x2[:,:,i]
        