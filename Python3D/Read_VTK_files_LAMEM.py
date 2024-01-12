# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 16:15:42 2022

@author: Andrea Piccolo
"""
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
from Parser_File import * 
from numba import jit
from numba import jitclass, types, typed


######################vtk utilities###########################################
### SURF

"""
Each pvtr/pvts file must have a different reader IF and only IF the size of the vector are different
more over, all the vtk bullshit remain in memory and i should find a way to get rid
from the naming space once for all and make it locals variable
"""


reader_s=vtk.vtkXMLGenericDataObjectReader()
VTK_SET_s=reader_s.SetFileName
VTK_UPDATE_s=reader_s.Update
VTK_RELEASE_s=reader_s.UnRegister
VTK_POINTS_s=vtk.vtkPoints
VTK_OUT_s=reader_s.GetOutput
AllPoints_s=VTK_POINTS_s()
GetData_s=AllPoints_s.GetData
nodes_vtk_array= GetData_s()
### PASSIVE_TRACER
reader_ptt = vtk.vtkXMLGenericDataObjectReader()
VTK_SET_ptt=reader_ptt.SetFileName
VTK_UPDATE_ptt=reader_ptt.Update
VTK_RELEASE_ptt=reader_ptt.UnRegister
VTK_POINTS_ptt=vtk.vtkPoints
VTK_OUT_ptt=reader_ptt.GetOutput
AllPoints_ptt=VTK_POINTS_ptt()
GetData_ptt=AllPoints_ptt.GetData
nodes_vtk_array_ptt= GetData_ptt()    
### Dyn
reader=vtk.vtkXMLGenericDataObjectReader()
VTK_SET_d=reader.SetFileName
VTK_UPDATE_d=reader.Update
VTK_RELEASE_d=reader.UnRegister
VTK_POINTS_d=vtk.vtkPoints
VTK_OUT_d=reader.GetOutput
AllPoints_d=VTK_POINTS_d()
GetData_d=AllPoints_d.GetData
nodes_vtk_array_d= GetData_d()
### Reader Phase (I HATE VTK DIOC~!)
reader_ph=vtk.vtkXMLGenericDataObjectReader()
VTK_SET_ph=reader_ph.SetFileName
VTK_UPDATE_ph=reader_ph.Update
VTK_RELEASE_ph=reader.UnRegister
VTK_POINTS_ph=vtk.vtkPoints
VTK_OUT_ph=reader_ph.GetOutput
AllPoints_ph=VTK_POINTS_d()
GetData_ph=AllPoints_ph.GetData
nodes_vtk_array_ph= GetData_ph()
##########################################

"""
class dictionary
-> string list defining the 



Let's draft a new classes that read the field and so forth 
=> => Always assumed to be 2D, later on working on other stuff to introduce more 3D compatible option
Class A () 
input: 
    -> name of the file, its extension
    -> (specific option to handle the surface)
init:
    -> Parse file -> Dictionary as touple class 
    -> create variables and initialize them 
    i.e. given a list of variable, automatically generate the variable and its inverse dictionary
    
Update: 
    -> loop over the dictionary of the class, take the dictionary 
    -> update 
    -> taylor 
    
parse_file
    -> How many field are present within pvts, pvtr, and so forth 
    -> generate variables list => send to init class
 

Class B ()



-------------------------
grid_dictionary        :
    
    
    
    
    
    
    
    
    
free_surface_dictionary:
phase_dictionary       : 
------------------------



Class VTK_LAMEM_Data(Coordinate_System = C,  extension = pvtr, )
        
    def __init__(Coordinate_System, timestep0, extension = pvtr)
    self.var= self._parse_file():
    
    loop over variable 
    exec(self.var = np.....)
    => create variable within the class in a dynamic way and intiialize them 
    
"""
        
def _parse_line(line,key):
    """
    Do a regex search against all defined regexes and
    return the key and match result of the first matching regex

    """
    rx_dict = {
        'timestep': re.compile(r'timestep=(\".*\") '),
        'file': re.compile(r'file=(\".*\")'),
        'WholeExtent': re.compile(r'WholeExtent=(\".*\")')
        }
    rx = rx_dict.get(key)
    match = rx.search(line)
    if match:
        return match
    return None



   # File list 
   

def _parse_grid(Filename,p):
    F=open(Filename,'r')
    d=F.readlines()
    F.seek(0)
    for line in d:
        key = 'WholeExtent'
        match = _parse_line(line,key)
        if match:
            nodes_num=(match.group(1))
            match=[]
            break
           
    F.close()
    nodes_num=nodes_num[1:-1]
    if p == 1:
        nodes_num = nodes_num[:-15]


    nodes_num=[int(s) for s in nodes_num.split(' ')]
    nx = nodes_num[1]
    ny = nodes_num[3]
    nz = nodes_num[5]
    return nx, ny, nz




def _file_list(fname):
   
   
    Time=[]    # Time of the numerical experiment
    Flist=[]   # File list 
   
    F=open(fname,'r')
    d=F.readlines()
    F.seek(0)
    for line in d:
        key = 'timestep'
        match = _parse_line(line,key)
        if match:
            Time.append(match.group(1))
            match =[] 
        key = 'file'
        match = _parse_line(line,key)
        if match:
            Flist.append(match.group(1))
    F.close()
    time=np.zeros(len(Time))
    n_tim_step=[]
    i=0
    for itime in Time:
        tm=itime[1:-1]
        time[i]=np.float(tm)
        i+=1

    for istep in Flist:
        index = [x.start() for x in re.finditer('_', istep)]
        istp=istep[index[0]+1:index[1]]
        n_tim_step.append(istp)
        
        
        
    return time, Flist, n_tim_step




class Coordinate_System():
    def __init__(self,filename,ptsave,xlim = (0,0), ylim = (0,0),zlim = (0,0)): #Private
        
        self.xt,self.yt,self.zt,self.nx,self.ny,self.nz = self._Initial_Grid_(filename[0],0)
        
        self.xP,self.yP,self.zP,self.nxP,self.nyP,self.nzP = self._Initial_Grid_(filename[1],1) # _phase_grid
        
        self.x,self.y,self.z,self.ind_x,self.ind_y,self.ind_z = self.__zoom_grid_(xlim,ylim,zlim,0)
        
        xlim_new = (np.min(self.x),np.max(self.x))
        ylim_new = (np.min(self.y),np.max(self.y))
        zlim_new = (np.min(self.z),np.max(self.z))
        self.xp,self.yp,self.zp,self.ind_xp,self.ind_yp,self.ind_zp = self.__zoom_grid_(xlim_new,ylim_new,zlim_new,1)

        
    def __zoom_grid_(self,xlim,ylim,zlim,p): #Private
        """
        
        Parameters 
        -----------
        xlim = touple of x lim
        ylim = touple of y lim 
        zlim = touple of z lim 
        
        Returns
        
        x,y,z zoomed grid
        
        """
        x = self.xt
        y = self.yt
        z = self.zt
        
        if p == 1:
            x = self.xP
            y = self.yP
            z = self.zP
            x = (x[1:]+x[0:-1])/2
            y = (y[1:]+y[0:-1])/2
            z = (z[1:]+z[0:-1])/2
        
        if(xlim[1]-xlim[0]>0.0):
        
            ind_x = (x>=xlim[0])&(x<=xlim[1])
            x = x[ind_x]
            
        else:
            
            ind_x = []
            
        if(ylim[1]-ylim[0]>0.0):
        
            ind_y = (y>=ylim[0])&(y<=ylim[1])
            y = y[ind_y]
            
        else:
            
            ind_y = []
            
        if(zlim[1]-zlim[0]>0.0):    
            
            ind_z = (z>=zlim[0])&(z<=zlim[1])   
            z = z[ind_z]
        
        else: 
        
            ind_z = []
            
            
            
        
        
        return x,y,z,ind_x,ind_y,ind_z
    
    
    def _Initial_Grid_(self,Filename_0,p):
        """
        

        Parameters
        ----------
        Filename_0 : [string, path]
            First timestep to retrieve the main data of the grid
        Read the vtk file from the grid values. This allow to retrieve the main
        grid data. 
        Modify the function as such it creates the proper layer for the Phase 
        routine
        
        Returns
        -------
        xd : [float, array]
            Taylor x axis
        yd : [float, array]
            y axis
        zd : [float, array]
            z axis
        nx : [int]
            number of node along x
        ny : [int]
            number of node along x
        nz : [int]
            number of node along x

        """
        if p == 0: 
            VTK_SET_d(Filename_0)
            VTK_UPDATE_d()
            VTK_OUT_d().GetPoints(AllPoints_d)  
            nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array_d)
        else:
            VTK_SET_ph(Filename_0)
            VTK_UPDATE_ph()
            VTK_OUT_ph().GetPoints(AllPoints_ph)   
            nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array_ph)
            

        xd,yd,zd= nodes_nummpy_array[:,0] , nodes_nummpy_array[:,1] , nodes_nummpy_array[:,2]
        nx, ny, nz=_parse_grid(Filename_0,p)
        if p ==1: 
            """
            Necessary to introduce such correction
            """
            nx += 1 
            nz += 1
            ny += 1 
        
        xd = xd.reshape(nz,ny,nx)
        yd = yd.reshape(nz,ny,nx)
        zd = zd.reshape(nz,ny,nx)
        
        xd = xd[0,0,:]
        yd = yd[0,:,0]
        zd = zd[:,0,0]
        
       
        
    
        return xd,yd,zd,nx,ny,nz

class Phase():
    def __init__(self,C,Phase_dic):
        ind_x = C.ind_xp
        ind_y = C.ind_yp
        ind_z = C.ind_zp
        tx    = np.sum(ind_x)
        ty    = np.sum(ind_y)
        tz    = np.sum(ind_z)
        self.Phase = np.zeros([tz,ty,tx],dtype = int)
        self.Phase_dic = Phase_dic
    def _update_phase(self,fileph,C):
        VTK_SET_ph(fileph)
        VTK_UPDATE_ph()
        P_vtk = reader_ph.GetOutput().GetCellData().GetArray("phase")
        ph = vtk_to_numpy(P_vtk)
        ph=ph.reshape(C.nzP-1,C.nyP-1,C.nxP-1) # The grid needs the node, but phases value are "cell" whise data (i.e. elementwise)
        ind_x = C.ind_xp
        ind_y = C.ind_yp
        ind_z = C.ind_zp
        ph = ph[:,ind_y,:]
        ph = ph[ind_z,:,:]
        ph = ph[:,:,ind_x]
        
        self.Phase = ph 
        
        return self 
    
    #def interpolate_dynamic_properties(self,V,C):
        
    
    
    def _plot_phase_field(self,C,Val,ptsave,ipic,t_cur,FS,tdimen):
                
        ptsave_b=os.path.join(ptsave,'phase')
        time_sim = "{:.3f}".format(t_cur)
        time_dimen = "{:.3f}".format(tdimen)
        tick=r'$Time$ = %s $Myrs$, $t^{\dagger}$ = %s' %(time_sim,time_dimen)
        
        ph = self.Phase 
        
        x  = C.xp 
        z  = C.zp 
        
        xd = C.x 
        zd = C.z 
        
        list = []
        labels =[]
        for v,k in self.Phase_dic.items():
            list.append(k[1])
            labels.append(k[0])
        
        cmap = colors.ListedColormap(list)
        
        
        
        if not os.path.isdir(ptsave_b):
            os.mkdir(ptsave_b)
        fig = plt.figure()
        fna='Fig'+str(ipic)+'.png'
        fn=os.path.join(ptsave_b,fna)
        #grid = plt.GridSpec(2, 3, wspace=0.4, hspace=0.3)
        ax1 = fig.add_axes([0.1, 0.5, 0.8, 0.4],
                   xticklabels=[])
        ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.4],
                           ylim=(-5.0, 5.0), xlim = (np.min(C.x),np.max(C.x)))  
        cf=ax1.pcolormesh(x, z, ph,cmap=cmap,vmin=0,vmax=len(labels), shading='flat')
        cf0 = ax1.contour(C.x,C.z,Val.OP,levels = [0.9],colors = "k",alpha=0.6,linewidths=1.0)
        cf1 = ax1.contour(C.x,C.z,Val.C1,levels = [0.9],colors = "k",alpha=0.6,linewidths=1.0)
        cf1 = ax1.contour(C.x,C.z,Val.C2,levels = [0.9],colors = "k",alpha=0.6,linewidths=1.0)
        cf0 = ax1.contour(C.x,C.z,Val.CC1,levels = [0.9],colors = "k",alpha=0.6,linewidths=1.0)
        cf0 = ax1.contour(C.x,C.z,Val.CC2,levels = [0.9],colors = "k",alpha=0.6,linewidths=1.0)

        X,Y = np.meshgrid(C.x,C.z)
        
        cf2 = ax1.quiver(X[::20,::20], Y[::20,::20], Val.vx[::20,::20]/(Val.vm[::20,::20]), Val.vz[::20,::20]/(Val.vm[::20,::20]),units='width',width=0.0011)
        ax1.set_aspect(0.5*(np.max(x)-np.min(x)) / float(np.max(z)-np.min(z)))
        ax1.tick_params(axis='both', which='major', labelsize=5)
        ax1.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
        ax1.set_title(tick)
        ax1.set_ylim(np.min(z),np.max(z))
        ax1.set_xlim(np.min(x),np.max(x))
        topo = FS.Amplitude
        iteration = np.linspace(9,0,num=10)
        for ip in range(10):
            it = ipic - ip
            alpha_v= 0.8-ip*(1/12)
            if ip == 0: 
                cc = 'r'
            else:
                cc = 'b'
            if (it == 0) & (ip == 0) :
                ax2.plot(C.x, topo[:,0],c = cc,alpha = alpha_v,linewidth=alpha_v)
                break
            if (ip >0 ) & (it == 0 ):
                ax2.plot(C.x, topo[:,it],c = cc,alpha = alpha_v,linewidth=alpha_v)
                break 
            else: 
                ax2.plot(C.x, topo[:,it],c = cc,alpha = alpha_v,linewidth=alpha_v)
        ax2.plot(C.x, topo[:,ipic],c = 'r',alpha = 1.0,linewidth=1.2)
        
        ax2.set_aspect(0.5*(np.max(x)-np.min(x)) / float(10))
        plt.grid(True)
        plt.xlabel('x, [km]')
        plt.ylabel('H, [km]')
        ax2.tick_params(axis='both', which='major', labelsize=5)
        ax2.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
        
        ###############################################################       
        plt.draw()    # necessary to render figure before saving
        fig.savefig(fn,dpi=300,transparent=False)
        fig.clear()
        plt.close()
        

    
        
class VAL():
    def __init__(self,C,dictionary):
        
        ind_x = C.ind_x
        ind_y = C.ind_y
        ind_z = C.ind_z
        tx    = np.sum(ind_x)
        ty    = np.sum(ind_y)
        tz    = np.sum(ind_z)
        
        self.dict = dictionary
        self.OP    = np.zeros([tz,ty,tx],dtype=float)
        self.C1    = np.zeros([tz,ty,tx],dtype=float)
        self.C2    = np.zeros([tz,ty,tx],dtype=float)
        self.CC1    = np.zeros([tz,ty,tx],dtype=float)
        self.CC2    = np.zeros([tz,ty,tx],dtype=float)
        self.Lit    = np.zeros([tz,ty,tx],dtype=float)
        self.Sed    = np.zeros([tz,ty,tx],dtype=float)
        self.T     = np.zeros([tz,ty,tx],dtype=float)
        self.vx,self.vy,self.vz,self.vm = np.zeros([tz,ty,tx],dtype=float),np.zeros([tz,ty,tx],dtype=float),np.zeros([tz,ty,tx],dtype=float),np.zeros([tz,ty,tx],dtype=float)
        self.dx,self.dy, self.dz,self.dm = np.zeros([tz,ty,tx],dtype=float),np.zeros([tz,ty,tx],dtype=float),np.zeros([tz,ty,tx],dtype=float),np.zeros([tz,ty,tx],dtype=float)
        self.tauxx,self.tauzz,self.tauxz = np.zeros([tz,ty,tx],dtype=float),np.zeros([tz,ty,tx],dtype=float),np.zeros([tz,ty,tx],dtype=float)
        self.epsxx,self.epszz,self.epsxz = np.zeros([tz,ty,tx],dtype=float),np.zeros([tz,ty,tx],dtype=float),np.zeros([tz,ty,tx],dtype=float)
        self.vis  = np.zeros([tz,ty,tx],dtype=float)
        self.nu   = np.zeros([tz,ty,tx],dtype=float) 
        self.eps = np.zeros([tz,ty,tx],dtype=float)
        self.tau =  np.zeros([tz,ty,tx],dtype=float) 
        self.Rho   =np.zeros([tz,ty,tx],dtype=float)   
        self.gamma   =np.zeros([tz,ty,tx],dtype=float)
        self.Psi   =np.zeros([tz,ty,tx],dtype=float)           
        self.LGV         = ["tau","nu",'vz','vm',"gamma","eps","T","Psi"]
        self.Label       = [ r"$\tau^{\dagger}_{II} []$", r"$\Psi []$",r'$v_z [cm/yr]$',r'$v_m [cm/yr]$',"$\gamma [n.d.]$","$log_{10}(\dot{\epsilon^{\dagger}}_{II})$ $[]$","$T^{\dagger} []$"]
        self.Colormap    = ["cmc.bilbao","cmc.devon","cmc.broc","cmc.bilbao","cmc.nuuk","cmc.lapaz","cmc.lapaz","cmc.turku","cmc.cork","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao"]
        self.Val         = [(0.5,6.0),
                            ("min","max"),
                            ("min","max"),
                            ("min","max"),
                            (10**(-3),10**1),
                            (1e-4,1e2),
                            (0,1),
                           ("min","max"),
                            ("min","max"),
                            ("min","max"),
                            ("min","max")]
       
    def _update_Val(self,Filename_dyn,C,IC):
        
        VTK_SET_d(Filename_dyn)
        VTK_UPDATE_d()
        diction = self.dict 
        
        for k in diction:
            key = str(k)
            f   = str(self.dict[k])
            buf=reader.GetOutput().GetPointData().GetArray(f)
            buf=vtk_to_numpy(buf)
            if(key == "tau"):
                buf = buf#/(IC.tau0/1e6)
            if(key == "T"):
                buf = (buf)#-IC.T_av)/(IC.Tc)
            if(key == 'eps'):
                buf = buf#/IC.epsc
            if((key == "velocity") | (key == "disp")):
                if key == "velocity":
                    self.vx,self.vy,self.vz,self.vm = self.taylor_grid(C,key,buf)
                else: 
                    self.dx,self.dy,self.dz,self.dm = self.taylor_grid(C,key,buf)
            elif((key == "stress_T") | (key == "eps_T")):
               
               if (key == "stress_T"):
                   self.tauxx,self.tauzz,self.tauxz = self.taylor_grid(C,key,buf)
               else:
                   self.epsxx,self.epszz,self.epsxz = self.taylor_grid(C,key,buf)
            else: 
                eval(key,globals(),self.__dict__)[:,:,:] =self.taylor_grid(C,key,buf) 
            
        return self 
        
    def taylor_grid(self,C,key,buf): 
        
        nx = C.nx
        ny = C.ny
        nz = C.nz
        ind_x = C.ind_x
        ind_y = C.ind_y 
        ind_z = C.ind_z
        
        if ((key == "velocity") | (key == "disp")):
        
            t1,t2,t3= buf[:,0] , buf[:,1] , buf[:,2]
            t1 = t1.reshape([nz,ny,nx])
            t2 = t2.reshape([nz,ny,nx])
            t3 = t3.reshape([nz,ny,nx])
            
            tm=(t1**2+t2**2+t3**2)**(0.5)
            
            t1= t1[ind_z,:,:]
            t1= t1[:,ind_y,:]
            t1= t1[:,:,ind_x]

            
            
            
            t2 = t2[ind_z,:,:]
            t2 = t2[:,ind_y,:]
            t2 = t2[:,:,ind_x]
            
            
            t3 = t3[ind_z,:,:]
            t3 = t3[:,ind_y,:]
            t3 = t3[:,:,ind_x]
            
            tm = tm[ind_z,:,:]
            tm = tm[:,ind_y,:]
            tm = tm[:,:,ind_x]

            
            return t1,t2,t3,tm
        
        elif((key == "eps_T") | (key == "stress_T")):
        
            txx,txy,txz,txy,tyy,tyz,txz,tyz,tzz = buf[:,0],buf[:,1],buf[:,2],buf[:,3],buf[:,4],buf[:,5],buf[:,6],buf[:,7],buf[:,8]
            
            txx = txx.reshape([nz,ny,nx])
            tzz = tzz.reshape([nz,ny,nx])
            txz = txz.reshape([nz,ny,nx])
            
            
            txx = txx[ind_z,:,:]
            txx = txx[:,ind_y,:]
            txx = txx[:,:,ind_x]

            tzz = tzz[ind_z,:,:]
            tzz = tzz[:,ind_y,:]
            tzz = tzz[:,:,ind_x]

            
            txz = txz[ind_z,:,:]
            txz = txz[:,ind_y,:]
            txz = txz[:,:,ind_x]

            
            return txx,tzz,txz
        else:
            
            buf=buf.reshape([nz,ny,nx])
            buf = buf[ind_z,:,:]
            buf = buf[:,ind_y,:]
            buf = buf[:,:,ind_x]

            
            return buf 
        
    def _plot_maps_V(self,t_cur,z,x,ptsave,ipic):


        import cmcrameri.cm as cmc
        from matplotlib.colors import LogNorm
        ptsave_b=os.path.join(ptsave,"WhMaps")
        
        if not os.path.isdir(ptsave_b):
            os.mkdir(ptsave_b)
        
        values = self.LGV 
        index = self.Label
        cmaps = self.Colormap 
        LIM   = self.Val
        
        
        time_sim = "{:.3f}".format(t_cur)
      
        ic = 0  
        val = np.zeros((len(z),len(x)),dtype=float)
        fg = figure()
        ax0 = fg.gca()
        fna='Fig'+"{:03d}".format(ipic)+'.png'       
        cf =ax0.pcolormesh(x, z, val, shading='gouraud')
        cbar = fg.colorbar(cf, ax=ax0,orientation='horizontal',extend="both")
        cf0 = ax0.contour(x,z,self.OP,levels = [0.9],colors = "k",linewidths=0.75)
        cf1 = ax0.contour(x,z,self.C1,levels = [0.9],colors = "k",linewidths=0.75)
        cf1 = ax0.contour(x,z,self.C2,levels = [0.9],colors = "k",linewidths=0.75)
        cf0 = ax0.contour(x,z,self.CC1,levels = [0.9],colors = "k",linewidths=0.75)
        cf0 = ax0.contour(x,z,self.CC2,levels = [0.9],colors = "k",linewidths=0.75)

       
        for name in values:
        
            cmap2 = eval(cmaps[ic])

            val = eval(name,globals(),self.__dict__)
            log=0 
            if (name == "eps" )|(name =="gamma"):
                log = 1
                
            tick = r"t = %s [Myrs]" %(time_sim)
            
            ptsave_c=os.path.join(ptsave_b,name)
            
            if not os.path.isdir(ptsave_c):
            
                os.mkdir(ptsave_c)
        
            lm    = LIM[ic]
            
            
            lim_m = lm[0]
            lim_M = lm[1]
            
            val[abs(val) == np.inf] = np.nan
            
            if(lm[0]=="min"):
            
                lim_m = np.nanmin(val)

            if (lm[1]=="max"):
               
                lim_M = np.nanmax(val)
                
                if lm[0] != "min":
                     lim_m = lm[0]
                     
            if (np.isnan(lim_M)) | (np.isnan(lim_m))| (lim_m==lim_M): 
                lim_m = 0.01
                lim_M = +0.1
                

            fna='Fig'+"{:03d}".format(ipic)+'.png'
            fn = os.path.join(ptsave_c,fna)
           
            cf.set_array(val.ravel())
            cf.set_cmap(cmaps[ic])
            if log == 1:
                cf.norm = colors.LogNorm(vmin=lim_m, vmax=lim_M)
            else: 
                cf.norm=colors.Normalize(vmin=lim_m, vmax=lim_M)
                
            cf.set_clim([lim_m,lim_M])
            cbar.vmin = lim_m 
            cbar.vmax = lim_M
            cbar.update_normal(cf)
            ax0.set_aspect(0.5*(np.max(x)-np.min(x)) / float(np.max(z)-np.min(z)))
            ax0.tick_params(axis='both', which='major', labelsize=5)
            ax0.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
            ax0.set_ylim(np.min(z),np.max(z))
            ax0.set_xlim(np.min(x),np.max(x))
            cbar.set_label(index[ic])
            ax0.set_title(tick)

            #plt.draw()    # necessary to render figure before saving
            
            fg.savefig(fn,dpi=300,transparent=False)
            
            
            val = [] 

            ic +=1 
        fg.clear
        plt.close()
        
      
        
class FS():
  
    def __init__(self,C,nstep):
    
        tx  = np.sum(C.ind_x)
        ty  = np.sum(C.ind_y)
        
        self.vx = np.zeros((ty,tx,nstep),dtype=float)
        self.vy = np.zeros((ty,tx,nstep),dtype=float)
        self.vz = np.zeros((ty,tx,nstep),dtype=float)
        self.vm = np.zeros((ty,tx,nstep),dtype=float)
        self.vx = np.zeros((ty,tx,nstep),dtype=float)
        self.Topo = np.zeros((ty,tx,nstep),dtype=float)
        self.Amplitude = np.zeros((ty,tx,nstep),dtype=float)
        
    def _Update_(self,Filename_s,C,ipic):
        
        ind_x = C.ind_x
        ind_y = C.ind_y
        nx    = C.nx
        ny    = C.ny
        x     = C.x 
        y     = C.y 
        
        self.vx[:,:,ipic],self.vy[:,:,ipic],self.vz[:,:,ipic],self.vm[:,:,ipic] = self._Read_Field_Surf(Filename_s,"velocity [cm/yr]",ind_x,ind_y,x,y,nx,ny)
        self.Topo[:,:,ipic]               = self._Read_Field_Surf(Filename_s,"topography [km]",ind_x,ind_y,x,y,nx,ny)
        self.Amplitude[:,:,ipic]          = self._Read_Field_Surf(Filename_s,"amplitude [km]",ind_x,ind_y,x,y,nx,ny)
    
    def _Read_Field_Surf(self, Filename_s,Field,ind_x,ind_y,x,y,nx,ny):
    
        VTK_SET_s(Filename_s)
        VTK_UPDATE_s()
        buf=reader_s.GetOutput().GetPointData().GetArray(Field)
        buf=vtk_to_numpy(buf)
      
        
        if (Field == "velocity [cm/yr]"):
        
            vxS = buf[:,0]
            vyS = buf[:,1] 
            vzS = buf[:,2]
            vxS = vxS.reshape([ny,nx])
            vyS = vxS.reshape([ny,nx])
            vzS = vzS.reshape([ny,nx])
            vxS = vxS[ind_y,:]
            vxS = vxS[:,ind_x]

            vyS = vyS[ind_y,:]
            vyS = vyS[:,ind_x]
            
            vzS = vzS[ind_y,:]
            vzS = vzS[:,ind_x]
            vmS = (vxS**2+vzS**2)**0.5
            
            return vxS,vyS,vzS,vmS
        
        else:
            
            buf =  buf.reshape([ny,nx])
            buf = buf[ind_y,:]
            buf = buf[:,ind_x]

            return buf 
    
    def ASCI_FILE(self,ipic,t_cur,Test_Name,ptsave,C:Coordinate_System):
        
        """
        Write a simple ascii file for the post processing of the free surface data
        This is for the the free surface data, later on I will dedicate a bit of 
        more time on the usage of the passive tracers.     
        """
        file_name = str(ipic).zfill(7)+'__'+Test_Name+'Free_surface_data.txt'
    
        ptsave_b=os.path.join(ptsave,'DataBase_FS')
        if not os.path.isdir(ptsave_b):
            os.mkdir(ptsave_b)
        
        filename = os.path.join(ptsave_b,file_name)
        Y,X = np.meshgrid(C.y,C.x)
        buf_x = X.ravel()
        buf_y = Y.ravel()
        vx    = self.vx[:,:,ipic]
        vy    = self.vy[:,:,ipic]
        vz    = self.vz[:,:,ipic]
        vm    = self.vm[:,:,ipic]
        H     = self.Amplitude[:,:,ipic]
        
        S        = np.array([buf_x,buf_y,vx.ravel(),vy.ravel(),vz.ravel(),vm.ravel(),H.ravel()])
    
        if(os.path.isfile(filename)):
            os.remove(filename)

        f = open(filename, 'a+')

        f.write('########################################\n')
        f.write('time [Myrs] time step []\n')
        f.write('x, y, v_x,v_y ,v_z, v_m, Topography\n')
        f.write('  [km],  [cm/yrs], [cm/yrs],[cm/yrs], [km]\n')
        f.write('########################################\n')
        f.write('time = %6f, timestep = %d\n' %(t_cur,ipic))
        f.write('\n')
        np.savetxt(f, np.transpose(S),fmt='%.6f', delimiter=' ', newline = '\n') 
        f.close()
        print('Free surface data of the timestep %d, has been printed' %(ipic))
        
class Passive_Tracers():
    def __init__(self,Filename_ptr):
        self.n_marker= self._ReadPassive_(Filename_ptr,'none')
        self.x = np.zeros((self.n_marker),dtype=float)
        self.y = np.zeros((self.n_marker),dtype=float)
        self.z = np.zeros((self.n_marker),dtype=float)
        self.Ph =np.zeros((self.n_marker),dtype=int)
        self.T = np.zeros((self.n_marker),dtype=float)
        self.P =np.zeros((self.n_marker),dtype=float)
        self.ID = np.zeros((self.n_marker),dtype=int)
        
    def _update_PTracer(self,Filename_ptr):
        self = self._ReadPassive_(Filename_ptr,'Coordinate')
        self.Ph = self._ReadPassive_(Filename_ptr,'Phase')
        self.T = self._ReadPassive_(Filename_ptr,'Pressure [MPa]')
        self.P = self._ReadPassive_(Filename_ptr,'Temperature [C]')
        self.ID = self._ReadPassive_(Filename_ptr,'ID')
    
    def _ReadPassive_(self,Filename_ptr,Field):      
        VTK_SET_ptt(Filename_ptr)
        VTK_UPDATE_ptt()
        VTK_OUT_ptt().GetPoints()
        n_marker = VTK_OUT_ptt().GetNumberOfPoints()
        if Field == 'none':
            print('Number of Passive tracers: ', n_marker)
            return n_marker
        if(Field == 'Coordinate'):
            for i in range(self.n_marker):
                x,y,z = VTK_OUT_ptt().GetPoint(i)
                self.x[i] = x
                self.y[i] = y 
                self.z[i] = z
            return self 
        else:
            buf_ptt=vtk_to_numpy(VTK_OUT_ptt().GetPointData().GetArray(Field))
            return buf_ptt

@jit(nopython=True)   
def find1Dnodes(cord,cordm,number):
    # I was fucking up a few stuff:
    #buf = cordm-cord 
    #min = np.min(np.abs(buf))
    for index in range(number):
        if (cord[index]>cordm):
            break 
    return index-1    
    
@jit(nopython=True)
def findnodes(GrID,ix,iy,iz):
    intp1=GrID[iz,iy,ix]
    intp2=GrID[iz,iy,ix+1]
    intp3=GrID[iz,iy+1,ix+1]
    intp4=GrID[iz,iy+1,ix]
    intp5=GrID[iz+1,iy,ix]
    intp6=GrID[iz+1,iy,ix+1]
    intp7=GrID[iz+1,iy+1,ix+1]
    intp8=GrID[iz+1,iy+1,ix]
    return intp1,intp2,intp3,intp4,intp5,intp6,intp7,intp8

@jit(nopython=True)
def linearinterpolation(xx,x1,x2,intp1,intp2):
    wx=(xx-x1)/(x2-x1)
    R=intp1*(1-wx)+intp2*wx

    return R

@jit(nopython=True)
def bilinearinterpolation(xx,yy,x1,x2,y1,y2,intp1,intp2,intp3,intp4):
    wx=(xx-x1)/(x2-x1)
    wy=(yy-y1)/(y2-y1)

    # FUCKING BUG: intp4 -> 1-x*intp4
    R=intp1*(1-wx)*(1-wy)+intp2*wx*(1-wy)+intp3*wx*wy+intp4*wy*(1-wx)

    return R    

@jit(nopython=True)
def trilinearinterpolation(xx,yy,zz,x1,x2,y1,y2,z1,z2,intp1,intp2,intp3,intp4,intp5,intp6,intp7,intp8):
    wx=(xx-x1)/(x2-x1)
    wy=(yy-y1)/(y2-y1)
    wz=(zz-z1)/(z2-z1)    
    i1=bilinearinterpolation(xx,yy,x1,x2,y1,y2,intp1,intp2,intp3,intp4)
    i2=bilinearinterpolation(xx,yy,x1,x2,y1,y2,intp5,intp6,intp7,intp8)
    R=linearinterpolation(zz,z1,z2,i1,i2)

    return R    







