
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
from typing import Literal
from abc import ABCMeta, abstractmethod
#from Auxilary_function import *







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
class Data_Base(object):
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
                
        self.detachment_velocity = np.zeros([self.n_test],dtype = float)
        
        self.Starting_tearing   = np.zeros([self.n_test],dtype = float)
        
        self.Ending_tearing    = np.zeros([self.n_test],dtype = float)
        
        self.uplift            = np.zeros([self.n_test,3],dtype=float) 
        
        self.dt             = np.zeros([self.n_test,3],dtype=float) 
        
        self.Avolume = np.zeros([self.n_test],dtype = float)
        
        self.StressLimit = np.zeros([self.n_test],dtype = float)
        
        self.Temp      = np.zeros([self.n_test],dtype = float)
        
        self.tau_max   = np.zeros([self.n_test],dtype = float)
        
        self.latex_table = []
   
   
        
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
    
    
    
    """
    Function to save a smaller data base for the Karlsruhe and Glasgow group
    """
    def _write_h5_database(self,ptsave,TestName,T):
       data_name = "Data_base_KIT_GLA.hdf5"
       data_base_name = os.path.join(ptsave,data_name)
       print(data_base_name)
       f = h5py.File(data_base_name, 'a')
       node = "/"+TestName
       e = False
       if node in f.keys():
               f[node]
               e = True
       if e == False: 
           f.create_group(node)
       node_coordinate_system = node+"/"+"Coordinate_System"
       # function to save the initial coordinate system 
       f = self.save_test_data(f,node_coordinate_system,T.C)
       print("Size coordinate system is %2f" %(sys.getsizeof(f)))
       # function to save the slab detachment 
       node_S = node+"/"+"Det"
       f= self.save_test_data(f,node_S,T.Det)
       print("Size Det is %2f" %(sys.getsizeof(f)))
       # function to save the free surface
       node_FS = node+"/"+"FS"
       f= self.save_test_data(f,node_FS,T.FS)
       print("Size free surface is %2f" %(sys.getsizeof(f)))
       buf_name = node+"/time"
       if buf_name in f.keys():
           del f[buf_name]      # load the data
           f.create_dataset(buf_name,data = np.array(T.time))
       else:
           f.create_dataset(buf_name,data = np.array(T.time))
        
       print("Size free surface is %2f" %(sys.getsizeof(f)))

       f.close() 

    def save_test_data(self,f,path_DB,Type_DATA):

        keys_data=Type_DATA.__dict__.keys()
        # Loop over the several field of the sub classes 
        size_array = 0
        
        if isinstance(Type_DATA,FS):
            keys_data = ['dH','dH_fil','H','vz','vz_fil']
        
        for v in keys_data:
            if v != 'dict':
                buf_name = path_DB+"/"+v
                buf_cl  = eval(v,globals(),Type_DATA.__dict__)
                # Check if the node exist, and in case delete for overwritting it. 
                if buf_name in f.keys():
                    del f[buf_name]      # load the data                
                if (isinstance(buf_cl,str)) | (isinstance(buf_cl,list)):
                    f.create_dataset(buf_name,data=buf_cl)
                else:
                    f.create_dataset(buf_name,data=np.float32(buf_cl))
                size_array += sys.getsizeof(buf_cl)
        print('Size of the database %2f' %(size_array))
        return f

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
    @timer
    def _post_process_data(self,path:str,path_save:str,save_data:bool,print_topography:bool,check:bool):
        
        perform = True
        if check == True: 
            if 'General_DB' in self.Tests_Name[0][:]: 
                perform = False 
            else: 
                perform = True 
                
        if perform == True:
            itest = 0 
            self.latex_table = [['Test Name',r'$v_c$',r'$\tau_{lim}$',r'$V_{a,\mathrm{dis}}$',r'$v_{\mathrm{tearing}}$',r'$\dot{H}_{\mathrm{N}}$',r'$\dot{H}_{Te}$',r'$\dot{H}_{\mathrm{LT}}$']]
            for it in range(self.n_test):
                test_name = self.Tests_Name[it]
                if ((test_name[0] != 'PR_no') & (test_name[0] != 'General_DB')):
                    path_save_b = os.path.join(path_save,test_name[1])
                    if not os.path.isdir(path_save_b):
                        os.mkdir(path_save_b)

                    print(test_name)

                    test = Test(self.path,test_name)
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

                    self.uplift[itest,0] = test.FS.dH_uplift_mean[1,1]

                    self.uplift[itest,1] = test.FS.dH_uplift_mean[0,1]

                    self.uplift[itest,2]  = test.FS.dH_uplift_mean[2,1]

                    self.tau_max[itest] = test.Det.maxtau/(self.StressLimit[itest]/1e6)

                    print("Average uplift during necking stage is %.2f, during tearing stage is %.2f and the long term is %.2f"%(self.uplift[itest,0],self.uplift[itest,1],self.uplift[itest,2]))

                    ipic = 0 
                    if save_data == True:
                        self._write_h5_database(path_save,test_name[1],test)
                    if print_topography == True:
                        test._print_topographic_data_ASCI(test_name[1],path_save_b)

                    self.latex_table.append([test_name[1],f"{self.Temp[itest]:.1f}",f"{self.StressLimit[itest]/1e6:.1f}",f"{self.Avolume[itest]:.1f}",f"{mean_v:.1e}",f"{self.uplift[itest,1]:.1e}",f"{self.uplift[itest,0]:.1e}",f"{self.uplift[itest,2]:.1e}"])
                #self.create_latex_table(path_save,latex_table)
                itest = itest+1
            self._save_new_DB_voice()
        else: 
            self.read_DB_variable()
        print('Finished to collect the data')
        print('|= ============== =|')
        print(' =|==============|= ')
        print('_____Creating Table ____')
        self.create_latex_table(path_save)
        


            
    # Function to update the main database and saves the data that have been retrieved by the database reader
    # Save LateX Data Base 
    # function that saves a latex table for pubblication purposes 
    def create_latex_table(self,path_save):
        from texttable import Texttable
        import latextable
        table = Texttable()
        table.set_cols_align(["l", "c", "c","c","c","c","c","c"])
        table.set_cols_valign(["m", "m", "m","m","m","m","m","m"])
        table.add_rows(self.latex_table)
        print(table.draw())
        print('\nLatextable Output:')
        print(latextable.draw_latex(table, caption="Experiment list table", label="table:experiment list"))


        # Create rows 


    def _save_new_DB_voice(self):
        """
        Function that saves the post processed data into DB 
        """
        f = h5py.File(self.path, 'a')
        node = "/"+'General_DB'
        e = False
        if node in f.keys():
            f[node]
            e = True
        self.save_test_data(f,node,self)
        f.close() 
    
    def read_DB_variable(self): 
                
        f = h5py.File(self.path,'r')
        
        # path 2 variable
        
        keys_data=self.__dict__.keys()

        for k in keys_data:
            path2variable = "/%s/%s" %('General_DB',k)
            # Create array 
            self.__dict__[k] = np.array(f[path2variable])
        
        f.close()
        return self
        

        

@timer
class Test(Data_Base):
    def __init__(self,path:str,Test_name:str):

        super().__init__(path)
        
        self.Test_Name = Test_name[1]
        
        self.time = self._read_variable(['/time','Myr', 'Time vector'],Test_name)
        
        self.time_M = (self.time[0:1:-1]+self.time[1:1:])/2

        self.IC = IC(Test_name,self)

        self.C  = C(Test_name,self)
                
        self.Det = Det(Test_name,self)
        
        self.FS = FS(Test_name,self)
        
        self.Ptr = Ptr(Test_name,self)
        
        
    def print_topography_timeseries(self,path_save,coord_x):
        """
        Input:
        self. 
        path_save = destination where to save the profiles
        coord_x = what is the collect the profile
        Output: 
        -> txt file: 1D timeseries of the filtered uplift of the given profile along y
        -> txt file: 0D timeseries of the maximum of the filtered uplift for a given profile
        -> plot that illustrate where is the slab and the profile
        """
        def fmt(x):
            s = f"{x:.1f}"
            if s.endswith("0"):
                s = f"{x:.0f}"
            return rf"{s} km" if plt.rcParams["text.usetex"] else f"{s} km"


        save_folder = os.path.join(path_save,"Profiles")
        if not os.path.isdir(save_folder):
            os.mkdir(save_folder)
        
        save_folder = os.path.join(path_save,"Profiles",self.Test_Name)
        if not os.path.isdir(save_folder):
            os.mkdir(save_folder)

        
        file_name = os.path.join(save_folder,'Profile_%d'%int(coord_x))
        #[Prepare the time vector]
        
        time_v = self.time
        
        start_detachment = np.min(self.Det.det_vec[(self.C.x_sp>=100) & (self.C.x_sp<=1100)])
        end_detachment   = np.max(self.Det.det_vec[(self.C.x_sp>=100) & (self.C.x_sp<=1100)])

        det_ = np.asarray(time_v*0,dtype=int)
        det_[(time_v>=start_detachment) & (time_v<=end_detachment)] = 1
        #[Select the coordinate where to envaluate]
        # Select the indexes that are higher than this value 
        ind_x = np.where(self.C.xg >= coord_x)
        # Select the index that best represent the coord_x
        ind_x = ind_x[0][0]
        # Select the time series 
        time_series_1D = self.FS.dH_fil[:,ind_x,:]
        T,Y = np.meshgrid(time_v,self.C.yg)
        S = np.array([T.ravel(),Y.ravel(),time_series_1D.ravel()])
        if(os.path.isfile(file_name)):
            os.remove(file_name)
        f = open(file_name, 'a+')
        f.write('########################################\n')
        f.write('time [Myr], y [km], dH [mm/yr]\n')
        f.write('The profile is take perpendicular to the trench at x = %3f [km]'%coord_x)
        f.write('dH array that has been filtered with a median \n')
        f.write('filter (scipy) with a kernel size of 7.\n')
        f.write('########################################\n')
        f.write('\n')
        np.savetxt(f, np.transpose(S),fmt='%.6f', delimiter=' ', newline = '\n') 
        f.close()
        
        # Print maximum uplift for the profile 
        time_0D_series = np.zeros(len(time_v),dtype=float)
        for i in range(len(time_v)):
            time_0D_series[i] = np.max(time_series_1D[:,i])
        
        file_name = file_name+'_0D_time_series'
        
        if(os.path.isfile(file_name)):
            os.remove(file_name)
        f = open(file_name, 'a+')
        S = []
        S = np.array([time_v,time_0D_series,det_])
        f.write('########################################\n')
        f.write('time [Myr], max(dH) [mm/yr], Det (1=YES,0 = NO)\n')
        f.write('The profile is take perpendicular to the trench at x = %3f [km], and the value represents the maximum uplift'%coord_x)
        f.write('dH array that has been filtered with a median \n')
        f.write('filter (scipy) with a kernel size of 7.\n')
        f.write('########################################\n')
        f.write('\n')
        np.savetxt(f, np.transpose(S),fmt='%.6f', delimiter=' ', newline = '\n') 
        f.close()
        
        fg = figure()
        
        fn = os.path.join(save_folder,'Simplified_initial_topographyx_%d_km'%int(coord_x))
        ax = fg.gca()
        p1 = ax.contour(self.C.xg,self.C.yg,self.FS.H[:,:,10],levels = [-2.0,-1.5,-0.5,0.0,0.5,1.0,1.5,2.0],colors = 'k',linewidths=0.5)
        ax.clabel(p1, p1.levels, inline=True, fmt=fmt, fontsize=6)
        ax.axvline(x=coord_x,linewidth = 1.8,color='firebrick',label = r'Profile')
        ax.plot(self.C.x_trench_p,self.C.y_trench_p,linewidth = 1.5,linestyle = 'dashdot',label = r'Slab position',color = 'rebeccapurple')
        ax.set_xlabel(r'$x$/[km]')
        ax.set_ylabel(r'$y$/[km]')
        ax.legend(loc='upper right')
        ax.xaxis.set_tick_params(labelsize=18)
        ax.yaxis.set_tick_params(labelsize=18)
        
        fg.savefig(fn,dpi=600,transparent=False)

    def _print_topographic_data_ASCI(self,Testname:str,path_save):
        
        itime = len(self.time)
        
        # Print Time vector
        self.ASCI_time_Vec(Testname,path_save)
        
        for i in range(itime):
            self.ASCI_FILE_ALT(i,self.time[i],Testname,path_save)
        
        print('Topographic data of % has been printed.'%(Testname))

    def ASCI_FILE_ALT(self,ipic,t_cur,Test_Name,ptsave):

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
        Y,X = np.meshgrid(self.C.xg,self.C.yg)
        buf_x = X.ravel()
        buf_y = Y.ravel()
        vz_M    = self.FS.vz_M[:,:,ipic]
        dH    = self.FS.dH[:,:,ipic]
        H     = self.FS.H[:,:,ipic]
        dH_fil  = self.FS.dH_fil[:,:,ipic]
        vz_fil  = self.FS.vz_fil[:,:,ipic]
        S        = np.array([buf_x*1e3,buf_y*1e3,vz_M.ravel(),dH.ravel(),dH_fil.ravel(),vz_fil.ravel(),H.ravel()])
        if(os.path.isfile(filename)):
            os.remove(filename)
        f = open(filename, 'a+')
        f.write('########################################\n')
        f.write('time [Myrs] time step []\n')
        f.write('x, y,v_z,dHdt, dHdt_fil,v_z_fil ,Topography\n')
        f.write('  [m],[m],[mm/yrs],[mm/yrs],[mm/yrs],[mm/yrs], [m]\n')
        f.write('dH_fil and vz_fil: array that has been filtered with a median \n')
        f.write('filter (scipy) with a kernel size of 7.\n')
        f.write('########################################\n')
        f.write('time = %6f, timestep = %d\n' %(t_cur,ipic))
        f.write('\n')
        np.savetxt(f, np.transpose(S),fmt='%.6f', delimiter=' ', newline = '\n') 
        f.close()
        #print('Free surface data of the timestep %d, has been printed' %(ipic))

    def ASCI_time_Vec(self,Test_Name,ptsave):

        """
        Write a simple ascii file for the post processing of the free surface dat
        This is for the the free surface data, later on I will dedicate a bit of 
        more time on the usage of the passive tracers.     
        """
        file_name = 'Time_Vector'+'__'+Test_Name+'.txt'

        ptsave_b=os.path.join(ptsave,'DataBase_FS')
        if not os.path.isdir(ptsave_b):
            os.mkdir(ptsave_b)

        filename = os.path.join(ptsave_b,file_name)
        dt = np.diff(self.time)
        dt_s = 0.0*self.time
        dt_s[1:]=dt[:]
        S        = np.array([self.time,dt_s])

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

    



# Class containing the coordinate system information 
class C():
    def __init__(self,Test_name:str,T:Test):
        
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
        
        self.xg = T._read_variable(self.dict['xg'],Test_name)
        
        self.yg = T._read_variable(self.dict['yg'],Test_name)
        
        self.zg = T._read_variable(self.dict['zg'],Test_name)
        
        self.xp = T._read_variable(self.dict['xp'],Test_name)
        
        self.yp = T._read_variable(self.dict['yp'],Test_name)
        
        self.zp = T._read_variable(self.dict['zp'],Test_name)
        
        self.x_sp = T._read_variable(self.dict['x_sp'],Test_name)
        
        self.y_sp = T._read_variable(self.dict['y_sp'],Test_name)
        
        # Compute the initial position of the slab 
        
        y_1 = T._read_variable(self.dict['y_1'],Test_name)
        
        y_2 = T._read_variable(self.dict['y_2'],Test_name)

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
    def __init__(self,Test_name:str,T:Test):
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
        self.L0 = T._read_variable(self.dict['L0'],Test_name)
        self.D0 = T._read_variable(self.dict['D0'],Test_name)
        self.T_av = T._read_variable(self.dict['T_av'],Test_name)
        self.etarefS = T._read_variable(self.dict['etarefS'],Test_name)
        self.etarefM = T._read_variable(self.dict['etarefM'],Test_name)
        self.xiS = T._read_variable(self.dict['xiS'],Test_name)
        self.xiM = T._read_variable(self.dict['xiM'],Test_name)
        self.tau0 = T._read_variable(self.dict['tau0'],Test_name)
        self.VnM = T._read_variable(self.dict['VnM'],Test_name)
        self.VnS = T._read_variable(self.dict['VnS'],Test_name)
        self.tau_co = T._read_variable(self.dict['tau_co'],Test_name)
        coordinate = T._read_variable(self.dict['Coord_Slab'],Test_name)
        self.coordinate_Slab = coordinate[0:2] 


# Class containing the information related to the detachment
class Det():
    def __init__(self,Test_name:str,T:Test):
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
        self.D = T._read_variable(self.dict['D'],Test_name)
        self.Psi = T._read_variable(self.dict['Psi'],Test_name)
        self.T = T._read_variable(self.dict['T'],Test_name)
        self.tau_max = T._read_variable(self.dict['tau_max'],Test_name)
        self.depth_vec = T._read_variable(self.dict['depth_vec'],Test_name)
        self.det_vec = T._read_variable(self.dict['det_vec'],Test_name)
        self.tau_vec = T._read_variable(self.dict['tau_vec'],Test_name)
        self.x_vec = T._read_variable(self.dict['x_vec'],Test_name)
        self.y_vec = T._read_variable(self.dict['y_vec'],Test_name)
        self.vel_tear = T._read_variable(self.dict['vel_tear'],Test_name)
        self.y1 = T._read_variable(self.dict['x_slab1'],Test_name)
        self.y2 = T._read_variable(self.dict['x_slab2'],Test_name)
        # Derivative values 
        self.D_x_t_det = np.zeros([len(self.x_vec),len(T.time)],dtype = float)
        self.tau_x_t_det = np.zeros([len(self.x_vec),len(T.time)],dtype = float)
        self.Psi_det = np.zeros([len(self.x_vec),len(T.time)],dtype = float)
        self.T_det = np.zeros([len(self.x_vec),len(T.time)],dtype = float)
        # Filtered {these data are horrible to see without a moving average filter}
        self.deltaD = np.zeros([len(T.time)],dtype = float)
        self.meanD = np.zeros([len(T.time)],dtype = float)
        self.minD = np.zeros([len(T.time)],dtype = float)
        self.maxD = np.zeros([len(T.time)],dtype = float)
        self.maxtau = 0 

        condition = (T.C.x_sp > 100) & (T.C.x_sp<1100) & (self.depth_vec<-80)
        depth = self.depth_vec
        depth[depth>=-80]=np.nan
        ind = np.where((T.C.zp >= np.nanpercentile(self.depth_vec[condition==1],5)) &(T.C.zp <= np.nanpercentile(self.depth_vec[condition==1],95) ))
        
        self.time_evolution_necking(ind,(T.time),T.C.x_sp)
    
    def time_evolution_necking(self,i_along_x,time_d,x):
        """
        Function that simply select the nodes of the array that corresponds to the depth of detachment and saves the entire 
        timeseries.
        input:
        i_along_x = the depth index at which detachment occurs.
        itime = number of timestep
        output:
        self 

        """
        itime = len(time_d)
        
        start_detaching = np.nanmin(self.det_vec[(x>=100) & (x<=1100)])
        end_detaching = np.nanmax(self.det_vec[(x>=100) & (x<=1100)])
        
        for i in range(len(x)):
            for j in range(len(time_d)):
                self.D_x_t_det[i,j] = np.nanmean(self.D[i,i_along_x,j])
                self.tau_x_t_det[i,j] = np.nanmean(self.tau_max[i,i_along_x,j])
                self.T_det[i.j] = np.nanmean(self.T[i,i_along_x,j])
                self.Psi_det[i,j] = np.nanmean(self.Psi[i,i_along_x,j])

        max_tau = 0.0 
        for i in range(itime):
            a = self.D_x_t_det[:,i]
            a[a==-np.inf]=np.nan 
            a = a[(x>=100) & (x<=1100)]
            b = self.tau_x_t_det[:,i]
            if (np.nanmax(b) > max_tau) and (time_d[i]>1.0) :
                max_tau = np.nanmax(b)
                
            b=b[(x>=100) & (x<=1100)]
            n = int(np.floor(len(a[np.isnan(a)!=1]))/10)
            self.meanD[i]=np.nanmean(a)/100
            self.maxD[i]=np.nanmax(a)/100
            self.minD[i]=np.nanmin(a)/100
            
            if time_d[i] < start_detaching:
                self.deltaD[i] = np.nanmax(a)/100-np.nanmin(a)/100
            if time_d[i] >= start_detaching:
                self.minD[i] = 0.1

            a = []
            b = []
        self.maxtau = max_tau
        
        return self


# Class containing the information of the free surface
class FS():
    def __init__(self,Test_name:str,T:Test):
        self.dict = {'H': ['/FS/Amplitude','km', 'Amplitude'],
                     'dH': ['/FS/dH','mm/yr', 'Rate of variation of Amplitude with time'],
                     'vz_M': ['/FS/vz_M','mm/yr', 'Filtered v_z of free surface'],
                     'Thickness': ['/FS/thickness','km', 'Thickness of the lithosphere '],
                     'tau_mean': ['/FS/mean_stress','MPa', 'Mean stress'],
                     'eps' : ['/FS/mean_eps','1/s','strain rate'],
                     'vz' : ['/FS/vz','mm/yr','raw data'],
       }
        self.H = T._read_variable(self.dict['H'],Test_name)
        self.dH = T._read_variable(self.dict['dH'],Test_name)
        self.vz_M = T._read_variable(self.dict['vz_M'],Test_name)
        self.vz = T._read_variable(self.dict['vz'],Test_name)
        self.Thickness = T._read_variable(self.dict['Thickness'],Test_name)
        self.tau_mean = T._read_variable(self.dict['tau_mean'],Test_name)
        self.eps = T._read_variable(self.dict['eps'],Test_name)
        # Derivative Data Set [filtered and post processed]
        self.dH_fil =  np.zeros(np.shape(self.dH),dtype=float)
        self.vz_fil =  np.zeros(np.shape(self.dH),dtype=float)
       
        self.dH_long_term = np.zeros(np.shape(self.dH[:,:,0]),dtype = float)
        self.dH_detachment = np.zeros(np.shape(self.dH[:,:,0]),dtype = float)
       
        self.Uplift_LT = np.zeros(np.shape(self.dH[:,:,0]),dtype=float)
        self.Uplift_det = np.zeros(np.shape(self.dH[:,:,0]),dtype=float)

        self.total_uplift_LT = np.zeros(np.shape(self.dH[:,:,0]),dtype=float)
        self.total_uplift_NT = np.zeros(np.shape(self.dH[:,:,0]),dtype=float)
        self.total_uplift_Te = np.zeros(np.shape(self.dH[:,:,0]),dtype=float)
        self.dH_uplift_max = np.zeros([3,2],dtype=float)
        self.dH_uplift_min = np.zeros([3,2],dtype=float)
        self.dH_uplift_mean = np.zeros([3,2],dtype=float)
        self.time_series = []

        self.dH_fil = self.filter_array('dH')
        self.vz_fil = self.filter_array('vz_M')
        
        self.compute_total_uplift(T)
    
        self._compute_dH_tearing(T)
        
        
        
    
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
        i1 = np.where(T.time==np.nanmin(T.Det.det_vec[(T.C.x_sp>=100) & (T.C.x_sp<=1100) ]))
        # Find when do the tearing end: 
        i2 = np.where(T.time==np.nanmax(T.Det.det_vec[(T.C.x_sp>=100) & (T.C.x_sp<=1100) ]))
        # Filter isostatic 
        i_iso = np.where(T.time>1.0)
        i_iso = i_iso[0][0]

        dt_det = T.time[i2[0][0]]-T.time[i1[0][0]]
        dt_long = T.time[i2[0][0]]-T.time[i_iso]
        # Compute anomaly
        self.dH_detachment= self.H[:,:,i2]-self.H[:,:,i1]
        self.dH_long_term = self.H[:,:,i1]-self.H[:,:,i_iso]
       
        self.Uplift_det = self.dH_detachment/dt_det
        self.Uplift_LT  = self.dH_long_term/dt_long
        
        return self
    
    def compute_total_uplift(self,T:Test):
        
        def compute_cumulative_uplift(tm,H,t1,t2,C:C,ts:bool): 
            """
            Closure of the compute total uplift
            -> input argument: the time vector of the test 
            -> H the topography vector 
            -> t1 and t2 the time of the time interval 
            -> coordinate system 
            -> save or not save a timeseries 
            output if ts == False
            min_dH = 2 field array[0 = ]
            """
            tm_ = tm[t1:t2]
            H  =  H[:,:,t1:t2]
            dH = np.zeros(np.shape(H[:,:,0]),dtype=float)
            max_dH = -1000e6 
            min_dH = 1000e6
            dt_min = 0.0
            dt_max = 0.0
            mean_dH = np.zeros(2,dtype=float)
            max_dH  = np.zeros(2,dtype=float)
            min_dH  = np.zeros(2,dtype=float)
            if ts == True:
                time_series_cum = np.zeros([len(tm_)-1,2],dtype=float)
            for i in range(len(tm_)-1):
                dH += H[:,:,i+1]-H[:,:,i]
                if ts==True:
                    buf = _interpolate_2D(dH,C.xg,C.yg,C.x_trench_p,C.y_trench_p)
                    ind_g=(C.yg>=-100.0)&(C.yg<=100.0)
                    buf0 = dH[ind_g,:]
                    buf0 = buf0[:,C.ind_x_trench_g]
                    time_series_cum[i,1] = np.nanmean(buf0)
                    time_series_cum[i,0] = np.nanmean(buf)
                    buf= []
                    buf0 = []
            
            sc = 1000  #[Conversion km -> m]
            sc2 = 1000/1e6 #[conversion m/yr -> mm/yr]
            
            # Interpolate dH along trench direction 
            dH_trench = _interpolate_2D(dH,C.xg,C.yg,C.x_trench_p,C.y_trench_p)
            # Compute the average integral along trench direction 
            mean_dH[0] = np.nanmean(dH_trench)*sc
            max_dH[0] = np.nanmax(dH_trench)*sc
            min_dH[0] = np.nanmin(dH_trench)*sc
            dH *= sc 
            mean_dH[1] =mean_dH[0]*sc2/(np.max(tm_)-np.min(tm_))
            max_dH[1]  =  max_dH[0]*sc2/(np.max(tm_)-np.min(tm_))
            min_dH[1]  =  min_dH[0]*sc2/(np.max(tm_)-np.min(tm_))
            print('maximum uplift is %.2f / mm/yr ' %(max_dH[1]))
            print('min uplift is %.2f /mm/yr' %(min_dH[1]))
            if ts == False:
                return dH,max_dH,min_dH,mean_dH
            else:
                return dH,max_dH,min_dH,mean_dH,time_series_cum*sc
        
        
        
        simulation_st = np.where(T.time>1.0)
        simulation_st = simulation_st[0][0]
        i = np.where(T.time==np.nanmin(T.Det.det_vec[(T.C.x_sp>=100) & (T.C.x_sp<=1100) ]))
        start_tearing = i[0][0]
        i = np.where(T.time==np.nanmax(T.Det.det_vec[(T.C.x_sp>=100) & (T.C.x_sp<=1100) ]))
        end_tearing = i[0][0]
    # Compute the necking stage dH,max_dH,min_dH,mean_dH
        print('==================')
        print('Necking stage data')
        self.total_uplift_NT,self.dH_uplift_max[0,:],self.dH_uplift_min[0,:],self.dH_uplift_mean[0,:] = compute_cumulative_uplift(T.time,self.H,simulation_st,start_tearing,T.C,False)
        print('==================')
        

    # Compute the tearing stage dh ...
        print('==================')
        print('Tearing stage data')
        self.total_uplift_Te,self.dH_uplift_max[1,:],self.dH_uplift_min[1,:],self.dH_uplift_mean[1,:] = compute_cumulative_uplift(T.time,self.H,start_tearing,end_tearing,T.C,False)
        print('==================')
    # Compute the long term stuff
        print('==================')
        print('Long_term  data')
        self.total_uplift_LT,self.dH_uplift_max[2,:],self.dH_uplift_min[2,:],self.dH_uplift_mean[2,:],self.time_series = compute_cumulative_uplift(T.time,self.H,simulation_st,end_tearing,T.C,True)
        print('==================')

        return self 
        
        

# Class containing the Passive tracers information   
class Ptr(): 
    def __init__(self,Test_name:str,T:Test):
        self.dict = {'x': ['/PTrBas/x','km', 'x position'],
                     'y': ['/PTrBas/y','km', 'y position'],
                     'z': ['/PTrBas/z','km', 'z position'],
                     'P': ['/PTrBas/P','MPa', 'Pressure passive tracer '],
                     'T': ['/PTrBas/T','C', 'Temperature of passive tracer'],
                    }
        self.x = T._read_variable(self.dict['x'],Test_name)
        self.y = T._read_variable(self.dict['y'],Test_name)
        self.z = T._read_variable(self.dict['z'],Test_name)
        self.P = T._read_variable(self.dict['P'],Test_name)
        self.T = T._read_variable(self.dict['T'],Test_name)

        
def _merge_database(FileA:str,FileB:str,FileC:str,FileD:str,Dest_file:str):
        import h5py 

        with h5py.File(Dest_file,'w') as f_dest:
            with h5py.File(FileA,'r') as f_src:
                f_src.copy(f_src["/PR_r"],f_dest["/"],"PR_r")
                with h5py.File(FileB,'r') as f_src2:
                    f_src.copy(f_src2["/PR_200"],f_dest["/"],"PR_200")
                    with h5py.File(FileC,'r') as f_src3:
                        f_src3.copy(f_src3["/PR_600"],f_dest["/"],"PR_600")
                        with h5py.File(FileD,'r') as f_src4:
                            f_src4.copy(f_src4["/PR_no"],f_dest["/"],"PR_no")
                        
                        



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

def fmt2(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} $^\circ$C" if plt.rcParams["text.usetex"] else f"{s} degC"

