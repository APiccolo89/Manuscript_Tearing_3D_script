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
from Class_Data_Base import _interpolate_2D
from Class_Data_Base import timer
import matplotlib.cm as colmap

from time import perf_counter 
plugin='DICOM'
import cmcrameri as cmc
import scipy 
import os
import latex
import imageio
sys.path.insert(1, '../Python3D/')    
import  Read_VTK_files_LAMEM as RW
import  Slab_detachment      as SD 
from mycolorpy import colorlist as mcp

"""
Short disclaimer: The functions here are to plot the relevant data using 
the small database that I created from the real numerical simulation. During the 
several months of work and revision, the functions are not exactly well documented. 
Before the publication I will remove most of the deadly sins that I accumulated.
"""


os.environ["PATH"] += os.pathsep + '/usr/local/texlive/2015/bin/x86_64-darwin'
print(os.getenv("PATH"))

#matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#plt.rcParams.update(params)
matplotlib.rcParams['text.usetex'] = True

"""
Function that I found in stackoverflow: 
https://stackoverflow.com/questions/33737736/matplotlib-axis-arrow-tip
"""

def define_colorbar(cf,ax,lim:list,ticks:list,label:str):
        
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    
    cbaxes = inset_axes(ax, borderpad=1.3,  width="100%", height="15%", loc=3)   
    cbar=plt.colorbar(cf,cax=cbaxes, ticks=ticks, orientation='horizontal',extend="both")
    print(lim[0])
    print(lim[1])
    cbar.set_label(label=label,size=fnt_g.label_) 
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_coords(0.5,-0.15)
    cbar.ax.xaxis.set_tick_params(pad=0.1)
    cbar.ax.xaxis.set_label_position('bottom')
    cbar.ax.tick_params(labelsize=fnt_g.axis_)
    
    return cbaxes,cbar


def arrowed_spines(fig, ax):

    xmin, xmax = ax.get_xlim() 
    ymin, ymax = ax.get_ylim()

    # removing the default axis on all sides:
    for side in ['bottom','right','top','left']:
        ax.spines[side].set_visible(False)

    # removing the axis ticks
    #plt.xticks([]) # labels 
    #plt.yticks([])
    #ax.xaxis.set_ticks_position('none') # tick markers
    #ax.yaxis.set_ticks_position('none')

    # get width and height of axes object to compute 
    # matching arrowhead length and width
    dps = fig.dpi_scale_trans.inverted()
    bbox = ax.get_window_extent().transformed(dps)
    width, height = bbox.width, bbox.height

    # manual arrowhead width and length
    hw = 1./20.*(ymax-ymin) 
    hl = 1./20.*(xmax-xmin)
    lw = 1. # axis line width
    ohg = 0.3 # arrow overhang

    # compute matching arrowhead length and width
    yhw = hw/(ymax-ymin)*(xmax-xmin)* height/width 
    yhl = hl/(xmax-xmin)*(ymax-ymin)* width/height

    # draw x and y axis
    ax.arrow(xmin, 0, xmax-xmin, 0., fc='k', ec='k', lw = lw, 
             head_width=hw, head_length=hl, overhang = ohg, 
             length_includes_head= True, clip_on = False) 

    ax.arrow(0, ymax, 0., ymin-ymax, fc='k', ec='k', lw = lw, 
             head_width=yhw, head_length=yhl, overhang = ohg, 
             length_includes_head= True, clip_on = False)

class fnt_g():
    point_to_cm = 0.035277777777778
    label_ = 0.45/point_to_cm 
    axis_  = 0.4/point_to_cm
    title_ = 0.45/point_to_cm
    legend_= 0.35/point_to_cm

class MidpointNormalize(mpl.colors.Normalize):
    def __init__(self, vmin, vmax, midpoint=0, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        normalized_min = max(0, 1/ 2 * (1 - abs((self.midpoint - self.vmin)/ (self.midpoint - self.vmax))))
        normalized_max = min(1, 1/ 2 * (1 + abs((self.vmax - self.midpoint)/ (self.midpoint - self.vmin))))
        normalized_mid = 0.5
        x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
        return np.ma.masked_array(np.interp(value, x, y))
            
def figure1_par(path_save):
    
    """

    """
        
    def create_axis(ax, letter:str):
        plt.setp(ax.spines.values(), linewidth=1.4)
        ax.tick_params(axis="y",direction="in")
        ax.tick_params(axis="x",direction="in")
        ax.tick_params(width=1.2)
        ax.xaxis.set_tick_params(labelsize=fnt_g.axis_)
        ax.yaxis.set_tick_params(labelsize=fnt_g.axis_)
        ax.set_xlim([-480,480])
        ax.set_ylim([-400,50])
        props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)
        ax.text(0.90, 0.95, letter, transform=ax.transAxes, fontsize=fnt_g.label_,
            verticalalignment='top', bbox=props,color='white')
        #ax.set_xticks([])
        ax.yaxis.set_label_coords(-0.01,0.4)
        ax.set_yticks([-400,-200,0, 50])
        ax.set_ylabel(r'$z$ [km]',fontsize=fnt_g.label_)
        l_abelr = ['$-400$','','$0$','']
        ax.set_yticklabels(l_abelr)
        ax.tick_params(left=True,labelleft=True,labelbottom=False,right=True,top=True) 

    
        return ax
    
    def _import_profile():
        name_file = '../Initial_Setup/Timestep_00000001_1.00000000e-07/SDet_phase.pvtr'
        ptsave = ''
        Filename_0 = ['../Initial_Setup/Timestep_00000001_1.00000000e-07/SDet.pvtr',name_file]
        C = RW.Coordinate_System(Filename_0,ptsave,(-700.0,700.0),(-1500,500),(-600.0,50.0))
        x = C.xp 
        z = C.zp 
        phase_dictionary={
            "0" : ["Air","white"],
            "1" : ["Upper Crust", "orangered"],
            "2" : ["Lower Crust", "bisque"],
            "3" : ["Continental Lithosperhic mantle","orange"],
            "4" : ["Continental Lithosperhic mantle2","deepskyblue"],
            "5" : ["Upper Mantle", "aliceblue"],
            "6" : ["Slab","navy"],
            "7" : ["Oceanic Crust","c"],
            "8" : ["Sediments O crust","rosybrown"],
            "9" : ["Weak Zone","linen"],
            "10" : ["Lower  Crust 2","skyblue"],
            "11" : ["Upper Crust2","plum"],
            "12" : ["Prism","salmon"],
            "13" : ["Passive Margin","olive"],
            "14" : ["Background O. Lithosphere","palegreen"],
            "15" : ["Background O. Crust","peachpuff"],
            "16" : ["Metamorphosed Prism","palevioletred"],
            "17" : ["Flysh","black"],
        }

        Ph    = SD.Phase_det(C,phase_dictionary)
        
        # Find Sections: 
        Ph._update_phase(name_file,C)  
        ind_0 = np.where(C.xp>=-590.0)
        ind_1 = np.where(C.xp>=0.0)
        ind_2 = np.where(C.xp>=590.0)
        ind_0 = ind_0[0][0]
        ind_1 = ind_1[0][0]
        ind_2 = ind_2[0][0]
        
        yp,zp= np.meshgrid(C.yp,C.zp)

        
        Ph1 = Ph.Phase[:,:,ind_0]
        Ph2 = Ph.Phase[:,:,ind_1]
        Ph3 = Ph.Phase[:,:,ind_2]
        
        # Small correction: once I had a weak zone. Then I realized that 
        # for such setup is causing a slab reatreat. I decided to substitute the
        # phase with the continental lithosphere (Phase 3). I changed the direction
        # of the slab and I forgot to update it: Phase 3 and Phase 4 are equivalent 
        # in all the properties: THIS IS NOT AFFECTING THE RESULTS. I simply homogeneise 
        # the visualisation. 
        
        Ph1[(Ph1==3) & (yp<-100.0)] = 4
        Ph2[(Ph2==3) & (yp<-100.0)] = 4
        Ph3[(Ph3==3) & (yp<-100.0)] = 4
        
        
        list_c=[]
        for v,k in phase_dictionary.items():
            list_c.append(k[1])
        
        cmap = colors.ListedColormap(list_c)
        

        
        return C.yp,C.zp,Ph1,Ph2,Ph3,cmap 
    
    # Prepare axis 
    
    fn = os.path.join(path_save,'Initial_Phase.svg')
    fn2 = os.path.join(path_save,'Initial_Phase.png')

    
    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(16*cm, 12*cm))  
    
    bx = 0.54
    by = 0.12
    sx = 0.45
    dx = 0.1
    sy = 0.27
    dy = 0.02
    ax0 = fg.add_axes([bx, by+2*dy+2*sy, sx, sy])
    ax1 = fg.add_axes([bx, by+dy+sy, sx, sy])
    ax2 = fg.add_axes([bx, by, sx, sy])
    ax0 = create_axis(ax0,r'[b]')
    ax1 = create_axis(ax1,r'[c]')
    ax2 = create_axis(ax2,r'[d]')
    

    yp,zp,P1,P2,P3,cmap_phase = _import_profile()
    
    cf4 = ax0.pcolormesh(yp, zp, (P1),cmap=cmap_phase,vmin=0,vmax=17, shading='gouraud')
    cf5 = ax1.pcolormesh(yp, zp, (P2),cmap=cmap_phase,vmin=0,vmax=17, shading='gouraud')
    cf6 = ax2.pcolormesh(yp, zp, (P3),cmap=cmap_phase,vmin=0,vmax=17, shading='gouraud')
    
    ax2.tick_params(left=True,labelleft=True,labelbottom=True,right=True,top=True) 
    ax2.set_xlabel(r'$y$ [km]',fontsize=fnt_g.label_)
    
    fg.savefig(fn2,dpi=600,transparent=True)
    
def figure_S2(path_save):
        
    def create_axis(ax, top_row:bool,letter:str,right:bool):
        plt.setp(ax.spines.values(), linewidth=1.4)
        ax.tick_params(axis="y",direction="in")
        ax.tick_params(axis="x",direction="in")
        ax.tick_params(width=1.2)
        ax.xaxis.set_tick_params(labelsize=fnt_g.axis_)
        ax.yaxis.set_tick_params(labelsize=fnt_g.axis_)
        props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)

        
        if top_row == True & right == True:
            ax.spines['bottom'].set_color('white')
            ax.spines['top'].set_color('black') 
            ax.spines['right'].set_color('white')
            ax.spines['left'].set_color('black')
            ax.xaxis.tick_top()
            ax.tick_params(labelbottom=False,labeltop=True)
            ax.set_ylim(-130,0)
            ax.set_xlim(20,1350)
            arrowed_spines(fg, ax)
            ax.set_xticks([0.0,600.0,1200.0])
            ax.set_xticklabels([r'$0$','','$1200$']) 
            ax.set_yticks([-100,-30,0])
            ax.text(0.87, 0.90, letter, transform=ax.transAxes, fontsize=fnt_g.label_,
                verticalalignment='top', bbox=props,color='white')
        elif (top_row == False) & (right == True): 
            ax.set_xticks([-100.0,-50,0.0])
            ax.set_yticks([-400,-200,0])
            ax.set_xticklabels([r'$100$','','$0$'])
            ax.text(0.70, 0.95, letter, transform=ax.transAxes, fontsize=fnt_g.label_,
                verticalalignment='top', bbox=props,color='white')            
        else: 
            ax.text(0.70, 0.95, letter, transform=ax.transAxes, fontsize=fnt_g.label_,
                verticalalignment='top', bbox=props,color='white')
        #ax.set_xticks([])
            ax.set_yticks([-600,-300,0])    
        
        return ax
    # Initial variables 
    z = np.linspace(start=0.0,stop=150.0,num=100)
    d = np.linspace(start=0.0,stop=100.0,num=100)
    L0 = 400.0 
    T_continent = compute_geoterhm(35.0,20,1350,600,z)
    T_Oplate    = compute_OC_geoterhm(20,1350,z,1e-6,30e6*60*60*365.25*24)
    T_Oplate    = (T_Oplate[1:]+T_Oplate[0:-1])/2
    T_slab_v0,l,lc   = compute_McKenzie_Field(20,1350,2.5,d,L0,-100)
    T_slab_v1,l,lc   = compute_McKenzie_Field(20,1350,5.0,d,L0,-100)
    T_slab_v2,l,lc   = compute_McKenzie_Field(20,1350,10.0,d,L0,-100)
    
    
    
    # Figure s.s.
    
    fn = os.path.join(path_save,'Figure_S2.svg')
    fn2 = os.path.join(path_save,'Figure_S3.png')

    
    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(10*cm, 12*cm))  
    bx = 0.1
    by = 0.15
    sx1 = 0.34
    sx2 = 0.2
    dx = 0.1
    dx1 = 0.1
    sy1 = 0.25
    sy2 = 0.45
    dy = 0.05
    ax0 = fg.add_axes([bx, by+dy+sy2, sx1, sy1])
    ax1 = fg.add_axes([bx+dx+sx1, by+dy+sy2, sx1, sy1])
    ax2 = fg.add_axes([bx, by, sx2, sy2])
    ax3 = fg.add_axes([bx+dx1+sx2, by, sx2, sy2])
    ax4 = fg.add_axes([bx+2*dx1+2*sx2, by, sx2, sy2])
    ax5 = fg.add_axes([bx+sx2/2, 0, 2*sx2+dx1*2, by])
   
    # Detailed axis 
    """
    1. I tried to group the common properties of the axis and assign with the closure
    above, however, I cannot do miracles, and I introduce a few command here to finalize 
    the layout. 
    2. Are there any clever way? Because I feel stupid atm. 
    """
    ax0 = create_axis(ax0, True,'[a]',True)
    ax1 = create_axis(ax1, True,'[b]',True)
    ax2 = create_axis(ax2, False,'[c]',True)
    ax3 = create_axis(ax3, False,'[d]',True)
    ax4 = create_axis(ax4, False,'[e]',True)



    ax1.tick_params(left=True,right=False,labelleft=False) 
    ax3.tick_params(left=True,right=True,labelleft=False) 
    ax4.tick_params(left=True,right=True,labelleft=False) 

    
    
    ax0.set_xlabel(r'$T$ [$^{\circ}\mathrm{C}]$',fontsize=fnt_g.label_)
    ax0.set_ylabel(r'$z$ [km]',fontsize=fnt_g.label_)
    ax1.set_xlabel(r'$T$ [$^{\circ}\mathrm{C}]$',fontsize=fnt_g.label_)
    ax2.set_xlabel(r'$d$ [km]',fontsize=fnt_g.label_)
    ax2.set_ylabel(r'$\ell$ [km]',fontsize=fnt_g.label_)
    ax3.set_xlabel(r'$d$ [km]',fontsize=fnt_g.label_)
    ax4.set_xlabel(r'$d$ [km]',fontsize=fnt_g.label_)

    l_abel = ['$400$','','$0$']
    ax2.set_yticklabels(l_abel)
    ax0.yaxis.set_label_coords(-0.03,0.5)
    ax2.yaxis.set_label_coords(-0.03,0.5)
    ax0.xaxis.set_label_coords(0.5,1.2)
    ax1.xaxis.set_label_coords(0.5,1.2)
    ax2.xaxis.set_label_coords(0.5,-0.03)
    ax3.xaxis.set_label_coords(0.5,-0.03)
    ax4.xaxis.set_label_coords(0.5,-0.03)
    
    # Plot
    
    ax0.plot(T_continent[z<100],-z[z<100],linewidth = 2.0, color = 'red')
    z_alt = (z[1:]+z[:-1])/2
    
    ax1.plot(T_Oplate[z_alt<100],-z_alt[z_alt<100],linewidth = 2.0, color = 'blue')
    t_field = np.linspace(start=20,stop=1350,num=13)
    
    cf1=ax2.contourf(-d,-l,np.transpose(T_slab_v0),cmap = 'cmc.lipari',levels=t_field)
    ax2.axhline(y=-lc,color='lightgreen',linewidth=1.2)
    
    cf2=ax3.contourf(-d,-l,np.transpose(T_slab_v1),cmap = 'cmc.lipari',levels=t_field)
    ax3.axhline(y=-lc,color='lightgreen',linewidth=1.2)
    
    cf3=ax4.contourf(-d,-l,np.transpose(T_slab_v2),cmap = 'cmc.lipari',levels=t_field)
    ax4.axhline(y=-lc,color='lightgreen',linewidth=1.2)
    

    cbaxes,cbar = define_colorbar(cf3,ax5,[20,1350],[20,700,1350],r'T/$^{\circ}C$')
    ax5.axis('off')


    
    #fg.savefig(fn,dpi=600,transparent=True)
    fg.savefig(fn2,dpi=600,transparent=True)

def compute_geoterhm(m_d,TS,TP,T_m,z):
    gr1 = (T_m-TS)/(m_d-0)
    gr2 = (TP-T_m)/(100-m_d)
    T = 0.0*z 
    T[z<=m_d] = TS+z[z<m_d]*gr1 
    T[(z>m_d) & (z<=100.0)]  = T_m+(z[(z>m_d) & (z<=100.0)]-m_d)*gr2 
    T[z>100.0] = TP 
    return T 

def compute_OC_geoterhm(TS,TP,z,kappa,age):
    erf = 0.0*z
    z_c = z[z<=100]*1000
    T = z*0.0
    erf[z<=100] = scipy.special.erf(z_c/(2*(kappa*age)**0.5))
    T[z<=100] = TS+(TP-TS)*erf[z<=100]
    T[z>100]  = TP 
    return T 

def compute_McKenzie_Field(TS,TP,vl,d,L0,dc):

    k = 3
    rho = 3300 
    Cp  = 1050 
    D0 = max(d)
    l  = np.linspace(start = 0, stop = L0,num=50)
    vl = vl/100/(365.25*60*60*24)
    Re = (Cp*vl*rho*D0*1000)/(2*k)
    lc = compute_lc(L0,90*np.pi/180,dc,l)
    L,D = np.meshgrid(l,d)
    L = L/D0 
    D = (D)/D0
    T_mk = 0.0*L 
    T_hc = 0.0*L 
    weight = 0.0*L
    T_hc = compute_OC_geoterhm(TS,TP,D*D0,1e-6,30e6*60*60*24*365.25)
    n = 26 
    Sigma = np.zeros(np.shape(L),dtype = float)

    
    for ir in range(n):
        i = ir+1
        a = (-1)**(i)/(i*np.pi)
        b = (Re-(Re**2+i**2*np.pi**2)**(0.5))*L
        c = np.sin(i*np.pi*(1-abs(D)))
        e = np.exp(b)
        Sigma += a*e*c
    T_mk = (TP)+2*(TP-TS)*Sigma
    T_mk[T_mk<0] = TS
    
    weight = 0+(0.8-0.1)/((lc/D0+0.04)-0.01)*(L-0.1)
    weight[weight>1] = 1.0
    weight[weight<0] = 0.1
    T =  T_mk*weight + T_hc*(1-weight)
    return T,l,lc
  
def compute_lc(L0,theta,dc,l): 
    z_in = -50.0 
    dl = np.abs(np.mean(np.diff(l)))
    lc = 0.0 
    for i in range(len(l)-1):
        sint,cost = compute_ribe_theta(l[i],l[i+1],theta)
        
        z_in = z_in-dl*sint 
        if z_in < dc:
            lc = (l[i]+l[i+1])/2.0
            break 
    return lc 
        
def compute_ribe_theta(l0,l1,theta): 
    Lb  = 150.0 
    if l0>Lb:
        theta_m = theta
    else:
        theta_l0 = theta*l0**2*((3*Lb-2*l0))/(Lb**3);
        theta_l1 = theta*l1**2*((3*Lb-2*l0))/(Lb**3);
        theta_m = (theta_l0+theta_l1)/2
    
    return np.sin(theta_m),np.cos(theta_m)


@timer
def _make_gif(test,ptsave_b):

    test_gf = os.path.join(ptsave_b,test.Test_Name)
    if not os.path.isdir(test_gf): 
        os.mkdir(test_gf)    
    
    cm = 1/2.54  # centimeters in inches
    file_list = []

    check = 1

    for ipic in range(len(test.time)):
        if test.time[ipic]>=1.0:
            if check == 1: 
                lim_1 = np.nanmin(test.FS.dH[:,:,ipic::])
                lim_2 = np.nanmax(test.FS.dH[:,:,ipic::])
                lim_psi = [round(np.nanpercentile(np.log10(test.Det.Psi[:,:,ipic::]),40),0),round(np.nanpercentile(np.log10(test.Det.Psi[:,:,ipic::]),100),0)]
                check = 0


            dH_trench = np.zeros(len(test.C.x_trench_p),dtype=float)

            # interpolate dh along the trench (filtered)
            dH_trench = _interpolate_2D(test.FS.dH[:,:,ipic],test.C.xg,test.C.yg,test.C.x_trench_p,test.C.y_trench_p)


            fna='Fig'+str(ipic)+'.png'
            fg = figure(figsize=(14*cm, 10*cm))  
            tick=r'$Time = %.3f$/Myr' %(test.time[ipic])
            fn = os.path.join(test_gf,fna)
            ax1 = fg.add_axes([0.12, 0.7, 0.7, 0.2])
            ax0 = fg.add_axes([0.12, 0.1, 0.7, 0.56])
            #axcb = fg.add_axes([0.1, 0.01, 0.6, 0.10])
            ax1b = ax1.twinx()
            ax1.tick_params(axis="y",direction="in")
            ax1.tick_params(axis="x",direction="in")
            ax1.tick_params(left=True,right=False,labelbottom=False)
            ax1b.tick_params(axis="y",direction="in")
            ax1b.tick_params(axis="x",direction="in")
            ax1b.tick_params(left=True,right=False,labelbottom=False)
    
            ax1.spines['bottom'].set_color('k')
            ax1.spines['top'].set_color('w') 
            ax1.spines['right'].set_color('k')
            ax1.spines['left'].set_color('k')

            ax1b.spines['bottom'].set_color('k')
            ax1b.spines['top'].set_color('w') 
            ax1b.spines['left'].set_color('k')
            ax1b.spines['right'].set_color('k')
            if (np.isnan(lim_1) == False) and  (np.isnan(lim_2) == False):
    #secax_y0b1= ax0b.secondary_yaxis(1.2, functions=(celsius_to_anomaly, anomaly_to_celsius))
    #secax_y0b1set_ylabel(r'$T - \ overline{T}\ [^oC]$')
                ax1.plot(test.C.x_sp,dH_trench,color='k',linewidth=1.2)
                ax1.axhline(y=0.0, color = 'k', linestyle=':', linewidth = 0.4)
                ax1.set_yticks([round(lim_1,2),round(lim_2,2)])
            else:
                ax1.set_yticks([-1,1])

            ax1.set_ylabel(r'$\dot{H}$ [$\frac{\mathrm{mm}}{\mathrm{yr}}$]',size=fnt_g.label_)
            ax1b.plot(test.C.x_sp,test.Det.tau_x_t_det[:,ipic]/(test.IC.tau_co/1e6),color='forestgreen',linewidth=1.2) 
            ax1b.plot(test.C.x_sp,test.Det.D_x_t_det[:,ipic]/100,color='firebrick',linewidth=1.2) 
            ax1b.set_ylabel(r'$\tau^{\dagger}_{max}$ []',size=fnt_g.label_)
            ax1b.set_yticks([0.0,1.0]) 
            secax_y0b = ax1b.secondary_yaxis(1.12)
            secax_y0b.set_ylabel(r'$D^{\dagger}$ []',size=fnt_g.label_)
            secax_y0b.set_ylim([0.0,1.0])
            secax_y0b.spines['right'].set_color('white')
    
            ax1b.yaxis.label.set_color('forestgreen')
            ax1b.tick_params(axis='y', colors='k')
            secax_y0b.yaxis.label.set_color('firebrick')
            secax_y0b.tick_params(axis='y',color='white')
            ax1.yaxis.set_label_coords(-0.1,0.5)

            ax1b.yaxis.set_label_coords(1.03,0.5)
            secax_y0b.yaxis.set_label_coords(1.4,0.5)

           
            secax_y0b.yaxis.set_tick_params(labelsize=fnt_g.axis_)
            secax_y0b.set_yticklabels(['',''])


            
            plt.setp(ax0.spines.values(), linewidth=1.4)
            plt.setp(ax1.spines.values(), linewidth=1.4)
            plt.setp(ax1b.spines.values(), linewidth=1.4)

            
            buf = np.log10(test.Det.Psi[:,:,ipic]) # Hard coded value i know. 
            levels = np.linspace(np.round(0.1), np.round(0.75), num=10, endpoint=True, retstep=False, dtype=float)
            cmap2 = colmap.get_cmap("cmc.lipari",lut=10)
            cmap2.set_under("k")
            cf=ax0.pcolormesh(test.C.x_sp,test.C.zp,np.transpose(buf),cmap=cmap2,vmin = lim_psi[0], vmax=lim_psi[1])
            cbar = fg.colorbar(cf,ax=ax0,orientation='horizontal',label=r'$\dot{\Psi}$ [$\mathrm{\frac{W}{m^3}}$]',extend="both",shrink=0.5)
            
            depth = test.Det.depth_vec
            depth[depth>=-80]=np.nan
            
            condition = (test.C.x_sp > 100) & (test.C.x_sp <1100)
            cf2=ax0.plot(test.C.x_sp[condition==1],depth[condition==1],color='red',linewidth=1.0)
            #dummy_var = test.C.x_sp[condition==1]/test.C.x_sp[condition==1]
            #cf2=ax0.fill_between(test.C.x_sp[condition==1],dummy_var*np.nanmin(depth[condition==1]),dummy_var*np.nanmax(depth[condition==1]),color='forestgreen',alpha=0.2,linewidth=1.2)
            #p1 = ax0.contour(test.C.x_sp,test.C.zp,np.transpose(test.Det.T[:,:,ipic]),levels = [800,900,1000,1100,1200],colors = 'k',linewidths=0.5)
            #ax0.clabel(p1, p1.levels, inline=True, fmt=fmt2, fontsize=6)
            
            ax1.set_title(tick)
            ax0.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
            ax0.set_ylabel(r'$z$ [km]')
            ax0.set_xlabel(r'$x_{trech}$ [km]')
            ax0.tick_params(axis="y",direction="in")
            ax0.tick_params(axis="x",direction="in")
            ax1.tick_params(axis="y",direction="in")
            ax1.tick_params(axis="x",direction="in")
            ax1b.tick_params(axis="y",direction="in")
            ax1b.tick_params(axis="x",direction="in")
            ax1.tick_params(left=True,right=False,labelbottom=False) 
            ax1b.tick_params(left=False,right=True,labelbottom=False) 

            ax0.tick_params(left=True,right=True,labelbottom=True) 
            ax0.set_xticks([200,600,1000])
            ax0.set_yticks([-80,-200,-400])
            ax0.set_ylim([-500,-60])
            ax1.set_xticks([200,600,1000])
            ax0.set_xticklabels([r'$200$','','$1200$'])
            ax0.xaxis.set_label_coords(0.5,-0.02)
            ax0.yaxis.set_label_coords(-0.01,0.5)
            ax1.yaxis.set_label_coords(-0.02,0.5)
            ax1b.yaxis.set_label_coords(1.02,0.5)
            #plt.draw()    # necessary to render figure before saving
            fg.savefig(fn,dpi=600)
            plt.close() 
            file_list.append(fn)

        # make gif
    
def figure3_5(A,path_figure,figure_name,time,Flag):
    # figure name
    cm = 1/2.54  # centimeters in inches

    fn = os.path.join(path_figure,'%s.png'%(figure_name))
    
    fg = figure(figsize=(14*cm, 10*cm))  
    # Prepare all the axis 
    
    bx = 0.12
    by = 0.12
    dx = 0.0 
    dy = 0.05 
    sx = 0.7 
    sy = 0.2 
    
    ax0 = fg.add_axes([bx, by+dy*2+2*sy, sx, sy]) # t1 
    ax1 = fg.add_axes([bx, by+dy+sy, sx, sy])     # t2 
    ax2 = fg.add_axes([bx, by, sx, sy])           # t3 
    
    # Prepare data 
    
    t0  = time[0]
    t1  = time[1]
    t2  = time[2] 
    
    tick0 = r'$t = %.2f$ [Myr]' %(t0)
    tick1 = r'$t = %.2f$ [Myr]' %(t1)
    tick2 = r'$t = %.2f$ [Myr]' %(t2)
    
    ind0 = np.where(A.time>=t0)
    ind1 = np.where(A.time>=t1)
    ind2 = np.where(A.time>=t2)
    
    ind0 = ind0[0][0]
    ind1 = ind1[0][0]
    ind2 = ind2[0][0]
    
    dH0 = _interpolate_2D(A.FS.vz_fil[:,:,ind0],A.C.xg,A.C.yg,A.C.x_trench_p,A.C.y_trench_p)
    dH1 = _interpolate_2D(A.FS.vz_fil[:,:,ind1],A.C.xg,A.C.yg,A.C.x_trench_p,A.C.y_trench_p)
    dH2 = _interpolate_2D(A.FS.vz_fil[:,:,ind2],A.C.xg,A.C.yg,A.C.x_trench_p,A.C.y_trench_p)

    ax0b = ax0.twinx()
    ax1b = ax1.twinx()
    ax2b = ax2.twinx()



    plt.setp(ax0.spines.values(), linewidth=1.4)
    plt.setp(ax1.spines.values(), linewidth=1.4)
    plt.setp(ax2.spines.values(), linewidth=1.4)

    plt.setp(ax0b.spines.values(), linewidth=1.4)
    plt.setp(ax1b.spines.values(), linewidth=1.4)
    plt.setp(ax2b.spines.values(), linewidth=1.4)


    # Do axis 0 
    # Estetic stuff 
    ax0.tick_params(axis="y",direction="in")
    ax0.tick_params(axis="x",direction="in")
    ax0.tick_params(left=True,right=False,labelbottom=False)
    ax0b.tick_params(axis="y",direction="in")
    ax0b.tick_params(axis="x",direction="in")
    ax0b.tick_params(left=False,right=True,labelbottom=False)
    
    ax0.spines['bottom'].set_color('k')
    ax0.spines['top'].set_color('w') 
    ax0.spines['right'].set_color('k')
    ax0.spines['left'].set_color('k')

    ax0b.spines['bottom'].set_color('k')
    ax0b.spines['top'].set_color('w') 
    ax0b.spines['left'].set_color('k')
    ax0b.spines['right'].set_color('k')
    
    # Create the shading of the axis for highlighting the tearing 
    D_0 = A.Det.D_x_t_det[:,ind0]
    D_1 = A.Det.D_x_t_det[:,ind1]
    D_2 = A.Det.D_x_t_det[:,ind2]

    i0 = A.C.x_sp[(A.C.x_sp>=100) & (A.C.x_sp<=1100) & (np.isnan(D_0)==1)]
    i1 = A.C.x_sp[(A.C.x_sp>=100) & (A.C.x_sp<=1100) & (np.isnan(D_1)==1)]
    i2 = A.C.x_sp[(A.C.x_sp>=100) & (A.C.x_sp<=1100) & (np.isnan(D_2)==1)]
    x0 = []
    x1 = []
    x2 = []
    if len(i0) >0:
        x0 = [np.nanmin(i0),np.nanmax(i0)]
    if len(i1) >0:
        x1 = [np.nanmin(i1),np.nanmax(i1)]
    if len(i2) > 0:
        x2 = [np.nanmin(i2),np.nanmax(i2)]
    lim_1 = np.min([dH0,dH1,dH2])
    lim_2 = np.max([dH0,dH1,dH2])

    #secax_y0b = ax0b.secondary_yaxis(1.2, functions=(celsius_to_anomaly, anomaly_to_celsius))
    #secax_y0b.set_ylabel(r'$T - \overline{T}\ [^oC]$')
    ax0.plot(A.C.x_sp,dH0,color='k',linewidth=1.2)
    ax0.axhline(y=0.0, color = 'k', linestyle=':', linewidth = 0.4)
    ax0.set_yticks([round(lim_1,2),round(lim_2,2)])
    ax0.set_ylabel(r'$\dot{H}$ [$\frac{\mathrm{mm}}{\mathrm{yr}}]$',size=fnt_g.label_)
    ax0b.plot(A.C.x_sp,A.Det.tau_x_t_det[:,ind0]/(A.IC.tau_co/1e6),color='brown',linewidth=1.2) 
    ax0b.plot(A.C.x_sp,A.Det.D_x_t_det[:,ind0]/100,color='orange',linewidth=1.2) 
    ax0b.set_ylabel(r'$\tau^{\dagger}_{max}$ []',size=fnt_g.label_)
    ax0b.set_yticks([0.0,1.0]) 
    secax_y0b = ax0b.secondary_yaxis(1.12)
    secax_y0b.set_ylabel(r'$D^{\dagger}$ []',size=fnt_g.label_)
    secax_y0b.set_ylim([0.0,1.0])
    secax_y0b.spines['right'].set_color('white')
    
    ax0b.yaxis.label.set_color('brown')
    ax0b.tick_params(axis='y', colors='k')
    secax_y0b.yaxis.label.set_color('orange')
    secax_y0b.tick_params(axis='y',color='white')
    
    #ax0.xaxis.set_label_coords(0.5,-0.02)
    ax0.yaxis.set_label_coords(-0.1,0.5)

    ax0b.yaxis.set_label_coords(1.03,0.5)
    secax_y0b.yaxis.set_label_coords(1.4,0.5)

    ax0.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax0.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax0b.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax0b.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    secax_y0b.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    secax_y0b.set_yticklabels(['',''])
    
    #ax1 {I can do a function for doing this boring stuff, I am repeating a lot, but hey, easy}
    
    ax1.tick_params(axis="y",direction="in")
    ax1.tick_params(axis="x",direction="in")
    ax1.tick_params(left=True,right=False,labelbottom=False)
    ax1b.tick_params(axis="y",direction="in")
    ax1b.tick_params(axis="x",direction="in")
    ax1b.tick_params(left=False,right=True,labelbottom=False)
    
    ax1.spines['bottom'].set_color('k')
    ax1.spines['top'].set_color('w') 
    ax1.spines['right'].set_color('k')
    ax1.spines['left'].set_color('k')

    ax1b.spines['bottom'].set_color('k')
    ax1b.spines['top'].set_color('w') 
    ax1b.spines['left'].set_color('k')
    ax1b.spines['right'].set_color('k')
    
    #secax_y0b = ax0b.secondary_yaxis(1.2, functions=(celsius_to_anomaly, anomaly_to_celsius))
    #secax_y0b.set_ylabel(r'$T - \overline{T}\ [^oC]$')
    ax1.plot(A.C.x_sp,dH1,color='k',linewidth=1.2)
    ax1.axhline(y=0.0, color = 'k', linestyle=':', linewidth = 0.4)
    ax1.set_yticks([round(lim_1,2),round(lim_2,2)])
    ax1.set_ylabel(r'$\dot{H}$ [$\frac{\mathrm{mm}}{\mathrm{yr}}$]',size=fnt_g.label_)

    ax1b.plot(A.C.x_sp,A.Det.tau_x_t_det[:,ind1]/(A.IC.tau_co/1e6),color='brown',linewidth=1.2) 
    ax1b.plot(A.C.x_sp,A.Det.D_x_t_det[:,ind1]/100,color='orange',linewidth=1.2) 
    ax1b.set_ylabel(r'$\tau^{\dagger}_{max}$ [] ',size=fnt_g.label_)
    ax1b.set_yticks([0.0,1.0]) 
    secax_y1b = ax1b.secondary_yaxis(1.12)
    secax_y1b.set_ylabel(r'$D^{\dagger}$ [] ',size=fnt_g.label_)
    secax_y1b.set_ylim([0.0,1.0])
    secax_y1b.spines['right'].set_color('white')
    
    ax1b.yaxis.label.set_color('brown')
    ax1b.tick_params(axis='y', colors='k')
    secax_y1b.yaxis.label.set_color('orange')
    secax_y1b.tick_params(axis='y',color='white')
    
    #ax0.xaxis.set_label_coords(0.5,-0.02)
    ax1.yaxis.set_label_coords(-0.1,0.5)

    ax1b.yaxis.set_label_coords(1.03,0.5)
    secax_y1b.yaxis.set_label_coords(1.4,0.5)

    ax1.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1b.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1b.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    secax_y1b.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    secax_y1b.set_yticklabels(['',''])

    #ax2 

    ax2.tick_params(axis="y",direction="in")
    ax2.tick_params(axis="x",direction="in")
    ax2.tick_params(left=True,right=False,labelbottom=True)
    ax2b.tick_params(axis="y",direction="in")
    ax2b.tick_params(axis="x",direction="in")
    ax2b.tick_params(left=False,right=True,labelbottom=False)
    
    ax2.spines['bottom'].set_color('k')
    ax2.spines['top'].set_color('w') 
    ax2.spines['right'].set_color('k')
    ax2.spines['left'].set_color('k')

    ax2b.spines['bottom'].set_color('k')
    ax2b.spines['top'].set_color('w') 
    ax2b.spines['left'].set_color('k')
    ax2b.spines['right'].set_color('k')
    
    #secax_y0b = ax0b.secondary_yaxis(1.2, functions=(celsius_to_anomaly, anomaly_to_celsius))
    #secax_y0b.set_ylabel(r'$T - \overline{T}\ [^oC]$')
    ax2.plot(A.C.x_sp,dH2,color='k',linewidth=1.2)
    ax2.axhline(y=0.0, color = 'k', linestyle=':', linewidth = 0.4)
    ax2.set_yticks([round(lim_1,2),round(lim_2,2)])
    ax2.set_ylabel(r'$\dot{H}$ [$\frac{\mathrm{mm}}{\mathrm{yr}}$]',size=fnt_g.label_)
    ax2b.plot(A.C.x_sp,A.Det.tau_x_t_det[:,ind2]/(A.IC.tau_co/1e6),color='brown',linewidth=1.2) 
    ax2b.plot(A.C.x_sp,A.Det.D_x_t_det[:,ind2]/100,color='orange',linewidth=1.2) 
    ax2b.set_ylabel(r'$\tau^{\dagger}_{max}$ [] ',size=fnt_g.label_)
    ax2b.set_yticks([0.0,1.0]) 
    secax_y2b = ax2b.secondary_yaxis(1.12)
    secax_y2b.set_ylabel(r'$D^{\dagger}$ [] ',size=fnt_g.label_)
    secax_y2b.set_ylim([0.0,1.0])
    secax_y2b.spines['right'].set_color('white')
    
    ax2b.yaxis.label.set_color('brown')
    ax2b.tick_params(axis='y', colors='k')
    secax_y2b.yaxis.label.set_color('orange')
    secax_y2b.tick_params(axis='y',color='white')
    
    ax2.xaxis.set_label_coords(0.5,-0.02)
    ax2.yaxis.set_label_coords(-0.1,0.5)

    ax2b.yaxis.set_label_coords(1.03,0.5)
    secax_y2b.yaxis.set_label_coords(1.4,0.5)

    ax2.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax2.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax2b.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax2b.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    secax_y2b.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    secax_y2b.set_yticklabels(['',''])
    ax0.set_xticks([200,600,1000])
    ax1.set_xticks([200,600,1000])
    ax2.set_xticks([200,600,1000])
    
    if Flag == True:
        m0 = np.round(np.nanmin(dH0))-1
        M0 = np.round(np.nanmax(dH0))+1
        m1 = np.round(np.nanmin(dH1))-1
        M1 = np.round(np.nanmax(dH1))+1
        m2 = np.round(np.nanmin(dH2))-1
        M2 = np.round(np.nanmax(dH2))+1
    
    
        ax0.set_yticks([m0,M0])
        ax1.set_yticks([m1,M1])
        ax2.set_yticks([m2,M2])
    
        ax0.set_ylim([m0,M0])
        ax0.set_ylim([m0,M0])
        ax1.set_ylim([m1,M1])
        ax1.set_ylim([m1,M1])
        ax2.set_ylim([m2,M2])
        ax2.set_ylim([m2,M2])

        


    ax2.set_xticklabels([r'$200$','','$1000$'])
    ax2.set_xlabel(r'$x_{trench}$ [km]',size=fnt_g.label_)

    props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)
    ax0.text(0.025, 1.15, tick0, transform=ax0.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.025, 1.15, tick1, transform=ax1.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
    ax2.text(0.025, 1.15, tick2, transform=ax2.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')

    ax0.text(0.9, 1.15, '[a]', transform=ax0.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.9, 1.15,'[b]', transform=ax1.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
    ax2.text(0.9, 1.15, '[c]', transform=ax2.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')

    
    ax0.spines[['top']].set_visible(False)
    ax0b.spines[['left','top']].set_visible(False)
    ax1.spines[['top']].set_visible(False)
    ax1b.spines[['left','top']].set_visible(False)
    ax2.spines[['top']].set_visible(False)
    ax2b.spines[['left','top']].set_visible(False)

    #ax2.yaxis.set_label_coords(1.02,0.5)
    fg.savefig(fn,dpi=600)

def figure_S3(A,path_save):

    def fmt(x):
        s = f"{x:.1f}"
        if s.endswith("0"):
            s = f"{x:.0f}"
        return rf"{s} km" if plt.rcParams["text.usetex"] else f"{s} km"
        
    fg = figure()
        
    fn = os.path.join(path_save,'Figure_S3')
    ax = fg.gca()
    p1 = ax.contourf(A.C.xg,A.C.yg,A.FS.H[:,:,10],levels = [-3.5,-3.0,-2.5,-2.0,-1.5,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5],cmap = 'cmc.oleron',linewidths=0.5)
    #ax.clabel(p1, p1.levels, inline=True, fmt=fmt, fontsize=fnt_g.axis_)
    ax.plot(A.C.x_trench_p,A.C.y_trench_p,linewidth = 2.0,linestyle = 'dashdot',label = r'Slab position',color = 'firebrick')
    ax.set_xlabel(r'$x$ [km]',fontsize=fnt_g.label_)
    ax.set_ylabel(r'$y$ [km]',fontsize=fnt_g.label_)
    ax.legend(loc='upper right')
    ax.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    cbar0 = fg.colorbar(p1, ax=ax,orientation='horizontal',extend="both")
    cbar0.ax.tick_params(labelsize=fnt_g.legend_)
    cbar0.set_label(label=r'${{H}}$ [$\mathrm{km}$]',size=fnt_g.label_) 
        
    fg.savefig(fn,dpi=600,transparent=False)

def figure_S9(A,path_figure,figure_name,time):
    # figure name
    from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
    def setup_axix(ax,ylabel):
        plt.setp(ax.spines.values(), linewidth=1.4)
        
        ax.tick_params(axis="y",direction="in")
        ax.tick_params(axis="x",direction="in")
        ax.tick_params(left=True,right=False,labelbottom=False)
        ax.spines['bottom'].set_color('k')
        ax.spines['top'].set_color('w') 
        ax.spines['right'].set_color('k')
        ax.spines['left'].set_color('k')

        ax.yaxis.set_label_coords(-0.1,0.5)
        ax.xaxis.set_tick_params(labelsize=fnt_g.axis_)
        ax.yaxis.set_tick_params(labelsize=fnt_g.axis_)
        
        ax.set_xticks([200,600,1000])
        ax.tick_params(left=True,right=False,labelbottom=True)
        ax.spines[['top','right']].set_visible(False)
        ax.set_ylabel(ylabel,fontsize=fnt_g.label_)
        return ax 
    
    
    
    cm = 1/2.54  # centimeters in inches

    fn = os.path.join(path_figure,'%s.png'%(figure_name))
    
    fg = figure(figsize=(14*cm, 10*cm))  
    # Prepare all the axis 
    
    bx = 0.2
    by = 0.12
    dx = 0.0 
    dy = 0.05 
    sx = 0.7 
    sy = 0.2 
    
    ax0 = fg.add_axes([bx, by+dy*2+2*sy, sx, sy]) # t1 
    ax1 = fg.add_axes([bx, by+dy+sy, sx, sy])     # t2 
    ax2 = fg.add_axes([bx, by, sx, sy])           # t3 
    
    # Prepare data 
    
    t0  = time[0]

    
    ind0 = np.where((A.time>1.0) & (A.time<=t0))
    
    

    ax0=setup_axix(ax0,r'$\Psi$ [$\mathrm{W}\mathrm{m}^3$]')
    ax1=setup_axix(ax1,r'$T$ [$^{\circ}$C]')
    ax2=setup_axix(ax2,r'$\tau^{\dagger}_{\mathrm{max}}$ []')
    it = 0 
    for i in ind0[0]:
        alpha = 0.05+(0.55/(len(ind0[0])-1))*it
        lw = 0.05+(0.55/(len(ind0[0])-1))*it
        print(alpha)
        ax0.plot(A.C.x_sp,(A.Det.Psi_det[:,i]),color='k',linewidth=lw,alpha=alpha) 
        ax1.plot(A.C.x_sp,A.Det.T_det[:,i],color='firebrick',linewidth=lw,alpha=alpha) 
        ax2.plot(A.C.x_sp,A.Det.tau_x_t_det[:,i]/(A.IC.tau_co/1e6),color='forestgreen',linewidth=lw,alpha=alpha) 
        if i == np.max(ind0):
            ax0.plot(A.C.x_sp,np.log10(A.Det.Psi_det[:,i]),color='k',linewidth=1.2,alpha=1.0) 
            ax1.plot(A.C.x_sp,A.Det.T_det[:,i],color='firebrick',linewidth=1.2,alpha=1.0) 
            ax2.plot(A.C.x_sp,A.Det.tau_x_t_det[:,i]/(A.IC.tau_co/1e6),color='forestgreen',linewidth=1.2,alpha=1.0) 
    
        
        it = it+1 
        

    ax0.set_yscale('log')
    ax0.set_ylim([10**(-9),10**(-3)])
    ax0.yaxis.set_minor_locator(AutoMinorLocator())
    ax0.tick_params(which='minor', length=4, color='k')
    ax1.set_ylim([850,1000])
    ax0.set_xticklabels(['','',''])
    ax1.set_xticklabels(['','',''])
    ax2.set_xlabel(r'$x_{\mathrm{trench}}$ [km]',fontsize=fnt_g.label_)


    props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)
    ax0.text(0.025, 1.15, '[a]', transform=ax0.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.025, 1.15, '[b]', transform=ax1.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
    ax2.text(0.025, 1.15, '[c]', transform=ax2.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
    
   

    #ax2.yaxis.set_label_coords(1.02,0.5)
    fg.savefig(fn,dpi=600)
    
def figure_S4_S7(A,path_figure,figure_name,letters_):

    """
    This is a picture that I wanted to do after reading the complete suite of articles of Fox. The exhumation rate that he provides has a maximum 
    resolution of 2 Myr, which entails, that I can port my data using this resolution. 
    """

    # Preparing the variables: 
    time = A.time
    time_1D = A.FS.time_1D_series_c
    time_0D = A.FS.time_0D_series_c 
    
    time_1D_d = np.zeros(np.shape(time_1D),dtype = float)
    time_0D_d = np.zeros(np.shape(time_0D),dtype = float)
    max_up_x = np.zeros(len(time),dtype=float)
    
    time_max = np.floor(np.max(time))
    time_pv = np.arange(1,time_max+2,2)
    
    # Compute the discrete field of uplift

    
    ia = np.where(A.C.x_sp>=100)
    ia = ia[0][0]
    
    ib = np.where(A.C.x_sp>=500)
    ib = ib[0][0]
    
    ic = np.where(A.C.x_sp>=1000)
    ic = ic[0][0]

    t_v = []
    ua = []
    ub =[]
    uc =[]
    for i in range(len(time_pv)-1): 
    
        t0 = time_pv[i]
        t1 = time_pv[i+1]
        t_buf = (time>=t0) & (time<=t1)
        
        for ix in range(len(A.C.x_trench_p)):
            time_1D_d[ix,t_buf==1] = np.nanmean(time_1D[ix,t_buf==1])
            time_0D_d[t_buf==1]    = np.nanmean(time_0D[t_buf==1])
        t_v.append(t0)
        t_v.append(t1)
        
        ua.append(np.nanmean(time_1D[ia,t_buf==1]))
        ua.append(np.nanmean(time_1D[ia,t_buf==1]))
        
        ub.append(np.nanmean(time_1D[ib,t_buf==1]))
        ub.append(np.nanmean(time_1D[ib,t_buf==1]))
        
        uc.append(np.nanmean(time_1D[ic,t_buf==1]))
        uc.append(np.nanmean(time_1D[ic,t_buf==1]))

        
        
        
    
    cm = 1/2.54  # centimeters in inches

    fn = os.path.join(path_figure,'%s.png'%(figure_name))
    
    fg = figure(figsize=(14*cm, 12*cm))
    time_max = np.floor(np.nanmax(A.Det.det_vec[(A.C.x_sp>=100) & (A.C.x_sp<=1100)])+1)
    bx = 0.1
    by = 0.12
    dx = 0.14 
    dy = 0.02 
    sx = 0.30 
    sy1 = 0.8
    sy2 = 0.25
    if time_max>20:
        time_max = 30
        ticks = [4,8,12,16,20]
        ticks_label = [r"4",r"8",r"12",r"16",r"20"]
    if len(letters_)>0:
        time_max = 30

        
    # First panel with the pcolormap
    ax0 = fg.add_axes([bx, by, sx, sy1]) # t1 
    
    # 0D timeseries discrete
    ax1 = fg.add_axes([bx+sx+dx, by+2*sy2+2*dy, sx, sy2])           # t3 
    ax2 = fg.add_axes([bx+sx+dx, by+dy+sy2, sx, sy2])     # t2 
    ax3 = fg.add_axes([bx+sx+dx, by, sx, sy2])     # t2 
    
    p1 = ax0.pcolormesh(A.C.x_sp,time,np.transpose(time_1D_d[:-1, :-1]),shading='auto',cmap = 'cmc.lajolla',vmin = np.nanpercentile(time_1D_d[time_1D_d != 0.0],10), vmax=np.nanpercentile(time_1D_d[time_1D_d != 0.0],90))
    #ax.clabel(p1, p1.levels, inline=True, fmt=fmt, fontsize=fnt_g.axis_)
    ax0.axvline(100,linewidth = 2.0,linestyle = 'dashdot',label = r'[b]',color = 'blue')
    ax0.axvline(500,linewidth = 2.0,linestyle = 'dashdot',label = r'[c]',color = 'forestgreen')
    ax0.axvline(1000,linewidth = 2.0,linestyle = 'dashdot',label = r'[d]',color = 'palevioletred')
    ax0.legend(loc='upper center', bbox_to_anchor=(0.45, 1.12),ncol=3, columnspacing=0.05,handletextpad=0.1, shadow=True,fontsize=8)
    ax0.set_xticks([100,500,1000])

    ax0.set_xlabel(r'$x_{trench}$ [km]',fontsize=fnt_g.label_)
    ax0.set_ylabel(r'$time$ [Myr]',fontsize=fnt_g.label_)
    ax0.set_ylim(1,time_max)
   # if time_max == 20:
    #    ax0.set_yticks(ticks)
     #   ax0.set_yticklabels(ticks_label)

        

    
    ax0.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax0.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    cbar0 = fg.colorbar(p1, ax=ax0,orientation='horizontal',extend="both")
    cbar0.ax.tick_params(labelsize=fnt_g.legend_)
    cbar0.set_label(label=r'${{\dot{H}}}$ [$\mathrm{mm/yr}$]',size=fnt_g.label_) 
    
    pb1 = ax1.plot(time,time_1D[ia,:],color = 'indigo',linewidth=0.7,label='Raw Data')    
    pb=ax1.plot(t_v,ua,color = 'salmon',linewidth=1.4,label='Filtered Data')
    ax1.legend(loc='upper center', bbox_to_anchor=(0.45, 1.25),ncol=2, columnspacing=0.05,handletextpad=0.1, shadow=True,fontsize=8)

    #pb2 = ax1.plot(time,time_0D_d,color='gainsboro',linewidth=0.6,label='Average')
    ax1.set_xlim(1,time_max)
    ax1.set_ylim(np.nanmin(time_1D_d[time_1D_d != 0.0]),np.nanmax(time_1D_d[time_1D_d != 0.0]))
    
    pc1 = ax2.plot(time,time_1D[ib,:],color = 'indigo',linewidth=0.7,label='Raw Data')
    pc=ax2.plot(t_v,ub,color = 'salmon',linewidth=1.4,label='Filtered Data')
    #pc2 = ax2.plot(time,time_0D_d,color='gainsboro',linewidth=0.6,label='Average')
    ax2.set_xlim(1,time_max)
    ax2.set_ylim(np.nanmin(time_1D_d[time_1D_d != 0.0]),np.nanmax(time_1D_d[time_1D_d != 0.0]))

    pd1 = ax3.plot(time,time_1D[ic,:],color = 'indigo',linewidth=0.7,label='Raw Data')
    pd=ax3.plot(t_v,uc,color = 'salmon',linewidth=1.4,label='Filtered Data')
    #pd2 = ax3.plot(time,time_0D_d,color='gainsboro',linewidth=0.6,label='Average')
    ax3.set_xlim(1,time_max)
    ax3.set_ylim(np.nanmin(time_1D_d[time_1D_d != 0.0]),np.nanmax(time_1D_d[time_1D_d != 0.0]))
    
    plt.setp(ax0.spines.values(), linewidth=1.4)
    plt.setp(ax1.spines.values(), linewidth=1.4)
    plt.setp(ax2.spines.values(), linewidth=1.4)
    plt.setp(ax3.spines.values(), linewidth=1.4)
    
    ax1.spines[['top','right']].set_visible(False)
    ax2.spines[['top','right']].set_visible(False)
    ax3.spines[['top','right']].set_visible(False)
    
    ax3.set_xlabel(r'$time$ [Myr]',fontsize=fnt_g.label_)
    ax1.set_ylabel(r'$\dot{H}$ [mm/yr]',fontsize=fnt_g.label_)
    ax2.set_ylabel(r'$\dot{H}$ [mm/yr]',fontsize=fnt_g.label_)
    ax3.set_ylabel(r'$\dot{H}$ [mm/yr]',fontsize=fnt_g.label_)
    ax1.tick_params(left = True, right = False , labelleft = True , 
        labelbottom = False, bottom = True) 
    ax2.tick_params(left = True, right = False , labelleft = True , 
        labelbottom = False, bottom = True) 

    if len(letters_)==0:
        letters_ = ['[a]','[b]','[c]','[d]']


    props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)
    ax0.text(0.025, 0.95, letters_[0], transform=ax0.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.025, 0.95, letters_[1], transform=ax1.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
    ax2.text(0.025, 0.95, letters_[2], transform=ax2.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
    ax3.text(0.025, 0.95, letters_[3], transform=ax3.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
   

    
    fg.savefig(fn,dpi=600)


"""
Figure 3/4: 
Stress and Thickness profile of a given test. 
-> input parameter 
A = Test database
path_figure = path for the figure
figure_name = figure name
= 
output the figure. 
"""

def figure6(A,B,path_figure,figure_name,det):
    
    """
    make axis for each of the subplot
    
    """
    def create_axis(ax,name_fig:str,stage:str):
        ax.tick_params(axis="y",direction="in")
        ax.tick_params(axis="x",direction="in")
        
        plt.setp(ax.spines.values(), linewidth=1.4)
    
        ax.tick_params(width=1.2)
        
        ax.set_xlabel(r'$x$ [km]',fontsize=fnt_g.label_)
        ax.set_xlabel(r'$x$ [km]',fontsize=fnt_g.label_)
        ax.set_ylabel(r'$y$ [km]',fontsize=fnt_g.label_)
        ax.yaxis.set_label_coords(-0.01,0.6)
        ax.xaxis.set_label_coords(0.5,-0.01)
        ax.xaxis.set_tick_params(labelsize=fnt_g.axis_)
        ax.yaxis.set_tick_params(labelsize=fnt_g.axis_)
        
        ax.set_xticks([-500,500])
        ax.set_yticks([-500,500])
        #ax.set_xlim([-500,500])
        #ax.set_ylim([-500,500])
        
        ax.tick_params(left = True, right = True , labelleft = True , 
                labelbottom = True, bottom = True) 
    


        props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)
        ax.text(0.9, 0.95, name_fig, transform=ax.transAxes, fontsize=fnt_g.label_,
            verticalalignment='top', bbox=props,color='white')
        ax.text(0.01, 0.95, stage, transform=ax.transAxes, fontsize=8,
            verticalalignment='top', bbox=props,color='white')
        return ax    

    
    def create_figure(ax,a,b,c,lim,t_x,t_y):
        #lim = [-500, 2000]
        norm = MidpointNormalize(vmin=lim[0], vmax=lim[1], midpoint=10.0)
        levels = np.linspace(np.round(lim[0]),np.round(lim[1]),60)
        cf=ax.pcolormesh(a,b,c,cmap = colormap,norm=norm)#,levels = levels)
        pl = ax.plot(t_x,t_y,linewidth=1.3, color = 'k')
        
        return ax,cf,pl
    
    # figure name
    fn = os.path.join(path_figure,'%s.png'%(figure_name))
    
    # Prepare variables
    c0 = A.FS.total_uplift_NT
    c1 = B.FS.total_uplift_NT
    c2 = A.FS.total_uplift_Te
    c3 = B.FS.total_uplift_Te
    c4 = A.FS.total_uplift_LT
    c5 = B.FS.total_uplift_LT
    lim_0 = np.round([np.nanpercentile(c0,1),np.nanmean(c0),np.nanpercentile(c0,99)])
    lim_1 = np.round([np.nanpercentile(c1,1),np.nanmean(c1),np.nanpercentile(c1,99)])
    lim_2 = np.round([np.nanpercentile(c2,1),np.nanmean(c2),np.nanpercentile(c2,99)])
    lim_3 = np.round([np.nanpercentile(c3,1),np.nanmean(c3),np.nanpercentile(c3,99)])
    lim_4 = np.round([np.nanpercentile(c4,1),np.nanmean(c4),np.nanpercentile(c4,99)])
    lim_5 = np.round([np.nanpercentile(c5,1),np.nanmean(c5),np.nanpercentile(c5,99)])
    lim_def_0=np.round([np.min([lim_0[0],lim_1[0],lim_2[0],lim_3[0],lim_4[0],lim_5[0]]),np.max([lim_0[2],lim_1[2],lim_2[2],lim_3[2],lim_4[2],lim_5[2]])])
    a = A.C.xg
    b = A.C.yg 
    t_x = A.C.x_trench_p
    t_y = A.C.y_trench_p
    i10,i20 = np.where(A.time==np.nanmin(A.Det.det_vec)),np.where(A.time==np.nanmax(A.Det.det_vec))
    i11,i21 = np.where(B.time==np.nanmin(B.Det.det_vec)),np.where(B.time==np.nanmax(B.Det.det_vec))
    string_title0 = r'Fast Tearing'
    string_title1 = r'Slow Tearing'
    # Prepare figure layout 
    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(18*cm, 11*cm))  
    bx = 0.07
    by = 0.175
    sx = 0.40
    dx = 0.03
    sy = 0.25
    dy = 0.01
    ax0 = fg.add_axes([bx, by+2*dx+2*sy, sx, sy])
    ax1 = fg.add_axes([bx+sx+dx, by+2*dx+2*sy, sx, sy])
    ax2 = fg.add_axes([bx, by+1*dx+1*sy, sx, sy])
    ax3 = fg.add_axes([bx+sx+dx, by+1*dx+1*sy, sx, sy])
    ax4 = fg.add_axes([bx, by, sx, sy])
    ax5 = fg.add_axes([bx+sx+dx, by, sx, sy])
    ax6 = fg.add_axes([bx+sx/2, 0, sx+dx, by])

    
    colormap = 'cmc.oleron'
    
    # Axis Y 
    ax0 = create_axis(ax0,'$[a]$','Necking Stage, $\Delta t = %.2f$ [Myr]'%(A.time[i10[0][0]]))
    ax1 = create_axis(ax1,'$[d]$','Necking Stage, $\Delta t = %.2f$ [Myr]'%(B.time[i11[0][0]]))
    ax0.tick_params(left = True, right = True , labelleft = True , 
                labelbottom = False, bottom = True) 
    ax1.tick_params(left = True, right = True , labelleft = False , 
                labelbottom = False, bottom = True) 
    ax0.set_xlabel(r'',fontsize=fnt_g.label_)
    ax1.set_xlabel(r'',fontsize=fnt_g.label_)
    ax1.set_ylabel(r'',fontsize=fnt_g.label_)
    # Axis Y
    ax2 = create_axis(ax2,'$[b]$','Tearing Stage, $\Delta t = %.2f$ [Myr]'%(A.time[i20]-A.time[i10[0][0]]))
    ax3 = create_axis(ax3,'$[e]$','Tearing Stage, $\Delta t = %.2f$ [Myr]'%(B.time[i21]-B.time[i11[0][0]]))
    ax2.tick_params(left = True, right = True , labelleft = True , 
                labelbottom = False, bottom = True) 
    ax3.tick_params(left = True, right = True , labelleft = False , 
                labelbottom = False, bottom = True) 
    ax2.set_xlabel(r'',fontsize=fnt_g.label_)
    ax3.set_xlabel(r'',fontsize=fnt_g.label_)
    ax3.set_ylabel(r'',fontsize=fnt_g.label_)

    
    # Axis X and Y
    ax4 = create_axis(ax4,'$[c]$','Long Term, $\Delta t = %.2f$ [Myr]'%(A.time[i20]))
    ax5 = create_axis(ax5,'$[f]$','Long Term, $\Delta t = %.2f$ [Myr]'%(B.time[i21]))
    ax5.set_ylabel(r'',fontsize=fnt_g.label_)
    ax5.tick_params(left = True, right = True , labelleft = False , 
                labelbottom = True, bottom = True) 


    ax0,cf0,pl0 =create_figure(ax0,a,b,c0,lim_def_0,t_x,t_y)
    ax1,cf1,pl1 =create_figure(ax1,a,b,c1,lim_def_0,t_x,t_y)
    ax2,cf2,pl2 =create_figure(ax2,a,b,c2,lim_def_0,t_x,t_y)
    ax3,cf3,pl3 =create_figure(ax3,a,b,c3,lim_def_0,t_x,t_y)
    ax4,cf4,pl4 =create_figure(ax4,a,b,c4,lim_def_0,t_x,t_y)
    ax5,cf5,pl5 =create_figure(ax5,a,b,c5,lim_def_0,t_x,t_y)

       
    cbaxes,cbar = define_colorbar(cf5,ax6,lim_def_0,[lim_def_0[0],300,lim_def_0[1]],r'${{\Delta H}}$ [$\mathrm{m}$]')
    ax6.axis('off')
     

    fg.savefig(fn,dpi=600)

def make_figure6_sup(A,B,path_figure,figure_name,det):
    
    """
    make axis for each of the subplot
    
    """
    def create_axis(ax,name_fig:str,stage:str):
        ax.tick_params(axis="y",direction="in")
        ax.tick_params(axis="x",direction="in")
        
        plt.setp(ax.spines.values(), linewidth=1.4)
    
        ax.tick_params(width=1.2)
        
        ax.set_xlabel(r'$x$ [km]',fontsize=fnt_g.label_)
        ax.set_xlabel(r'$x$ [km]',fontsize=fnt_g.label_)
        ax.set_ylabel(r'$y$ [km]',fontsize=fnt_g.label_)
        ax.yaxis.set_label_coords(-0.01,0.6)
        ax.xaxis.set_label_coords(0.5,-0.01)
        ax.xaxis.set_tick_params(labelsize=fnt_g.axis_)
        ax.yaxis.set_tick_params(labelsize=fnt_g.axis_)
        
        ax.set_xticks([-500,500])
        ax.set_yticks([-500,500])
        #ax.set_xlim([-500,500])
        #ax.set_ylim([-500,500])
        
        ax.tick_params(left = True, right = True , labelleft = True , 
                labelbottom = True, bottom = True) 
    


        props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)
        ax.text(0.9, 0.95, name_fig, transform=ax.transAxes, fontsize=fnt_g.label_,
            verticalalignment='top', bbox=props,color='white')
        ax.text(0.01, 0.95, stage, transform=ax.transAxes, fontsize=8,
            verticalalignment='top', bbox=props,color='white')
        return ax    

    
    def create_figure(ax,a,b,c,lim,t_x,t_y):
        #lim = [-500, 2000]
        norm = MidpointNormalize(vmin=lim[0], vmax=lim[1], midpoint=0.0)
        levels = np.linspace(np.round(lim[0]),np.round(lim[1]),60)
        cf=ax.pcolormesh(a,b,c,cmap = colormap,norm=norm)#,levels = levels)
        pl = ax.plot(t_x,t_y,linewidth=1.3, color = 'k')
        
        return ax,cf,pl
    
    # figure name
    fn = os.path.join(path_figure,'%s.png'%(figure_name))
    
    i10,i20 = np.where(A.time==np.nanmin(A.Det.det_vec)),np.where(A.time==np.nanmax(A.Det.det_vec))
    i11,i21 = np.where(B.time==np.nanmin(B.Det.det_vec)),np.where(B.time==np.nanmax(B.Det.det_vec))

    
    # Prepare variables
    c0 = A.FS.total_uplift_NT/A.time[i10[0][0]]
    c1 = B.FS.total_uplift_NT/B.time[i11[0][0]]
    c2 = A.FS.total_uplift_Te/(A.time[i20]-A.time[i10[0][0]])
    c3 = B.FS.total_uplift_Te/(B.time[i21]-B.time[i11[0][0]])
    c4 = A.FS.total_uplift_LT/A.time[i20]
    c5 = B.FS.total_uplift_LT/B.time[i21]
    """
    Convert from m/Myr -> mm/yr 1m = 1000 mm, 1Myr = 1e6 yr 
    """
    
    c0 = c0*(1000/1e6)
    c1 = c1*(1000/1e6)
    c2 = c2*(1000/1e6)
    c3 = c3*(1000/1e6)
    c4 = c4*(1000/1e6)
    c5 = c5*(1000/1e6)
    
    
    lim_0 = np.round([np.nanpercentile(c0,20),np.nanmean(c0),np.nanpercentile(c0,90)])
    lim_1 = np.round([np.nanpercentile(c1,20),np.nanmean(c1),np.nanpercentile(c1,90)])
    lim_2 = np.round([np.nanpercentile(c2,20),np.nanmean(c2),np.nanpercentile(c2,90)])
    lim_3 = np.round([np.nanpercentile(c3,20),np.nanmean(c3),np.nanpercentile(c3,90)])
    lim_4 = np.round([np.nanpercentile(c4,20),np.nanmean(c4),np.nanpercentile(c4,90)])
    lim_5 = np.round([np.nanpercentile(c5,20),np.nanmean(c5),np.nanpercentile(c5,90)])
    lim_def_0=np.round([np.min([lim_0[0],lim_1[0],lim_2[0],lim_3[0],lim_4[0],lim_5[0]]),np.max([lim_0[2],lim_1[2],lim_2[2],lim_3[2],lim_4[2],lim_5[2]])])
    a = A.C.xg
    b = A.C.yg 
    t_x = A.C.x_trench_p
    t_y = A.C.y_trench_p
    string_title0 = r'Fast Tearing'
    string_title1 = r'Slow Tearing'
    # Prepare figure layout 
    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(18*cm, 11*cm))  
    bx = 0.07
    by = 0.175
    sx = 0.40
    dx = 0.03
    sy = 0.25
    dy = 0.01
    ax0 = fg.add_axes([bx, by+2*dx+2*sy, sx, sy])
    ax1 = fg.add_axes([bx+sx+dx, by+2*dx+2*sy, sx, sy])
    ax2 = fg.add_axes([bx, by+1*dx+1*sy, sx, sy])
    ax3 = fg.add_axes([bx+sx+dx, by+1*dx+1*sy, sx, sy])
    ax4 = fg.add_axes([bx, by, sx, sy])
    ax5 = fg.add_axes([bx+sx+dx, by, sx, sy])
    ax6 = fg.add_axes([bx+sx/2, 0, sx+dx, by])

    
    colormap = 'cmc.roma'
    
    # Axis Y 
    ax0 = create_axis(ax0,'$[a]$','Necking Stage, $\Delta t = %.2f$ [Myr]'%(A.time[i10[0][0]]))
    ax1 = create_axis(ax1,'$[d]$','Necking Stage, $\Delta t = %.2f$ [Myr]'%(B.time[i11[0][0]]))
    ax0.tick_params(left = True, right = True , labelleft = True , 
                labelbottom = False, bottom = True) 
    ax1.tick_params(left = True, right = True , labelleft = False , 
                labelbottom = False, bottom = True) 
    ax0.set_xlabel(r'',fontsize=fnt_g.label_)
    ax1.set_xlabel(r'',fontsize=fnt_g.label_)
    ax1.set_ylabel(r'',fontsize=fnt_g.label_)
    # Axis Y
    ax2 = create_axis(ax2,'$[b]$','Tearing Stage, $\Delta t = %.2f$ [Myr]'%(A.time[i20]-A.time[i10[0][0]]))
    ax3 = create_axis(ax3,'$[e]$','Tearing Stage, $\Delta t = %.2f$ [Myr]'%(B.time[i21]-B.time[i11[0][0]]))
    ax2.tick_params(left = True, right = True , labelleft = True , 
                labelbottom = False, bottom = True) 
    ax3.tick_params(left = True, right = True , labelleft = False , 
                labelbottom = False, bottom = True) 
    ax2.set_xlabel(r'',fontsize=fnt_g.label_)
    ax3.set_xlabel(r'',fontsize=fnt_g.label_)
    ax3.set_ylabel(r'',fontsize=fnt_g.label_)

    
    # Axis X and Y
    ax4 = create_axis(ax4,'$[c]$','Long Term, $\Delta t = %.2f$ [Myr]'%(A.time[i20]))
    ax5 = create_axis(ax5,'$[f]$','Long Term, $\Delta t = %.2f$ [Myr]'%(B.time[i21]))
    ax5.set_ylabel(r'',fontsize=fnt_g.label_)
    ax5.tick_params(left = True, right = True , labelleft = False , 
                labelbottom = True, bottom = True) 


    ax0,cf0,pl0 =create_figure(ax0,a,b,c0,lim_def_0,t_x,t_y)
    ax1,cf1,pl1 =create_figure(ax1,a,b,c1,lim_def_0,t_x,t_y)
    ax2,cf2,pl2 =create_figure(ax2,a,b,c2,lim_def_0,t_x,t_y)
    ax3,cf3,pl3 =create_figure(ax3,a,b,c3,lim_def_0,t_x,t_y)
    ax4,cf4,pl4 =create_figure(ax4,a,b,c4,lim_def_0,t_x,t_y)
    ax5,cf5,pl5 =create_figure(ax5,a,b,c5,lim_def_0,t_x,t_y)

       
    cbaxes,cbar = define_colorbar(cf5,ax6,lim_def_0,[lim_def_0[0],0,lim_def_0[1]],r'${{\Delta H}}$ [$\mathrm{m}$]')
    ax6.axis('off')
     

    fg.savefig(fn,dpi=600)




def figure8(DB,path_save,figure_name):
    """
    Major update: better to put everything in a row otherwise it is horrible in a paper
    
    """
    # figure name
    fn = os.path.join(path_save,'%s.png'%(figure_name))
    
    # Prepare variables
    vel_tearing = (DB.detachment_velocity)
    AVol        = DB.Avolume 
    T           = DB.Temp 
    SLim        = DB.StressLimit/1e6
    
    # Prepare axis of the figures 
    
    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(18*cm, 9*cm))  
    bx = 0.12
    by = 0.15
    sx = 0.27
    dx = 0.02
    sy = 0.69
    dy = 0.01
    ax0 = fg.add_axes([bx, by, sx, sy])
    ax1 = fg.add_axes([bx+dx+sx, by, sx, sy])
    ax2 = fg.add_axes([bx+2*dx+2*sx, by, sx, sy]) 
    
    colors = ['royalblue','goldenrod','tomato']
    label_fig = [r'$v_c = 10$ [$\frac{\mathrm{cm}}{\mathrm{yr}}$]',r'$v_c = 5.0$ [$\frac{\mathrm{cm}}{\mathrm{yr}}$]',r'$v_c = 2.5$ [$\frac{\mathrm{cm}}{\mathrm{yr}}$]']

    T_u = np.sort(np.unique(T))
    T_u = T_u[T_u>0.0]
    
    # find the global limit of the axis:
    min_ax = 1.0
    max_ax = np.nanmax(vel_tearing[(T > 0)])
    max_ax = max_ax+max_ax*0.2

    for i in range(len(T_u)):
        ax0.scatter(AVol[(SLim==200.0) & (T == T_u[i])],vel_tearing[(SLim==200.0) & (T == T_u[i])],s=50,c=colors[i],edgecolor = 'k',label=label_fig[i])
    
    for i in range(len(T_u)):
        ax1.scatter(AVol[(SLim==400.0) & (T == T_u[i])],vel_tearing[(SLim==400.0) & (T == T_u[i])],s=50,c=colors[i],edgecolor = 'k',label=label_fig[i])
    
    for i in range(len(T_u)):
        ax2.scatter(AVol[(SLim==600.0) & (T == T_u[i])],vel_tearing[(SLim==600.0) & (T == T_u[i])],s=50,c=colors[i],edgecolor = 'k',label=label_fig[i])
    
    ax1.legend(loc='upper center', bbox_to_anchor=(0.45, 1.20),ncol=3, columnspacing=0.10,handletextpad=0.01, shadow=True,fontsize=8)
    
    #ax0.set_ytick(0.1,0.5,0.9)
    
    ax0.set_ylim(min_ax,max_ax)
    ax1.set_ylim(min_ax,max_ax)
    ax2.set_ylim(min_ax,max_ax)
    
    ax0.tick_params(axis="y",direction="in")
    ax0.tick_params(axis="x",direction="in")
    ax0.tick_params(left=True,right=True,labelbottom=True,labelleft = True) 
    ax1.tick_params(left=True,right=True,labelbottom=True,labelleft = False) 
    ax2.tick_params(left=True,right=True,labelbottom=True,labelleft = False) 

    ax1.tick_params(axis="y",direction="in")
    ax1.tick_params(axis="x",direction="in")
    ax2.tick_params(axis="y",direction="in")
    ax2.tick_params(axis="x",direction="in")
    
    plt.setp(ax0.spines.values(), linewidth=1.4)
    plt.setp(ax1.spines.values(), linewidth=1.4)
    plt.setp(ax2.spines.values(), linewidth=1.4)


    ax0.tick_params(width=1.2)
    ax1.tick_params(width=1.2)
    ax2.tick_params(width=1.2)
    
    ax0.set_ylabel(r'$v_{tearing}$ [$\frac{\mathrm{cm}}{\mathrm{yr}}$]',fontsize=fnt_g.label_)

    ax0.set_xlabel(r'$V_{a,dis}$ [$\mu \frac{\mathrm{m}^3}{\mathrm{mol}}$]',fontsize=fnt_g.label_)
    ax1.set_xlabel(r'$V_{a,dis}$ [$\mu \frac{\mathrm{m}^3}{\mathrm{mol}}$]',fontsize=fnt_g.label_)
    ax2.set_xlabel(r'$V_{a,dis}$ [$\mu \frac{\mathrm{m}^3}{\mathrm{mol}}$]',fontsize=fnt_g.label_)

    
    ax0.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax0.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax2.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax2.yaxis.set_tick_params(labelsize=fnt_g.axis_)

    ax0.set_yscale('log')
    ax1.set_yscale('log')
    ax2.set_yscale('log')


    
    props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)
    ax0.text(0.87, 0.94, '$[a]$', transform=ax0.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.87, 0.94, '$[b]$', transform=ax1.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    ax2.text(0.87, 0.94, '$[c]$', transform=ax2.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    
    ax0.text(0.05, 0.10, r'$\tau_{lim} = 200 [MPa]$', transform=ax0.transAxes, fontsize=fnt_g.axis_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.05, 0.10, r'$\tau_{lim} = 400 [MPa]$', transform=ax1.transAxes, fontsize=fnt_g.axis_,
        verticalalignment='top', bbox=props,color='white')
    ax2.text(0.05, 0.10, r'$\tau_{lim} = 600 [MPa]$', transform=ax2.transAxes, fontsize=fnt_g.axis_,
        verticalalignment='top', bbox=props,color='white')
    

    
    fg.savefig(fn,dpi=600)

def figure9(DB,path_save,figure_name):
    from scipy import stats
    # figure name
    def linear_fuc(x,slope,intercept):
        return (x*slope+intercept)
    

    def linear_internal_regression(a,b,min_v,max_v):
        slope, intercept, r, p, std_err = stats.linregress(a, b)
        x = np.linspace((min_v),(max_v),num=1000)
        eq = linear_fuc(x,slope,intercept)
        return x,eq,eq+std_err,eq-std_err,slope,intercept

    
    # Prepare variables
    vel_tearing = (DB.detachment_velocity)
    AVol        = DB.Avolume 
    T           = DB.Temp 
    SLim        = DB.StressLimit/1e6
    
    
    # Prepare axis of the figures 
    
    cm = 1/2.54  # centimeters in inches
     
    bx = 0.15
    by = 0.15
    sx = 0.40
    dx = 0.02
    sy = 0.6
    dy = 0.01
    type = ['a','b','c']
    Uplift_discrete = DB.Uplift_Te_discrete
    for ip in range(1):
        fg = figure(figsize=(18*cm, 9*cm)) 
        ax0 = fg.add_axes([bx, by, sx, sy])
        ax1 = fg.add_axes([bx+sx+dx,by,sx,sy])

        fn = os.path.join(path_save,'%s%s.png'%(figure_name,type[ip]))

        uplift = DB.uplift[:,ip]
        colors = [mcp.gen_color(cmap="cmc.lipari",n=6)]
        colors=colors[0][:]
        label_fig = [r'$V_{a,dis} = 8$ [$\frac{\mathrm{cm^3}}{\mathrm{mol}}$]',r'$V_{a,dis} = 10$ [$\frac{\mathrm{cm^3}}{\mathrm{mol}}$]',r'$V_{a,dis} = 11$ [$\frac{\mathrm{cm^3}}{\mathrm{mol}}$]',r'$V_{a,dis} = 12$ [$\frac{\mathrm{cm^3}}{\mathrm{mol}}$]',r'$V_{a,dis} = 13 $ [$\frac{\mathrm{cm^3}}{\mathrm{mol}}$]',r'$V_{a,dis} = 15$ [$\frac{\mathrm{cm^3}}{\mathrm{mol}}$]']
        a=np.log10(vel_tearing[(AVol > 0)])
        b=np.log10(uplift[(AVol > 0) ])
        c=np.log10(Uplift_discrete[(AVol > 0) ])



        AVol_u = np.sort(np.unique(AVol))
        AVol_u = AVol_u[AVol_u>0.0]
        for i in range(len(AVol_u)):
            ax0.scatter(vel_tearing[(AVol == AVol_u[i])],uplift[(AVol == AVol_u[i])],c=colors[i],s=50,edgecolor = 'k',label=label_fig[i])
           
            ax1.scatter(vel_tearing[(AVol == AVol_u[i])],Uplift_discrete[(AVol == AVol_u[i])],c=colors[i],s=50,edgecolor = 'k',label=label_fig[i])

        v,u,ea1,ea2,slope,intercept=linear_internal_regression(a,b,np.min(a),np.max(a))
        ax0.plot(10**v,10**u,color = 'forestgreen',linewidth=1.2)
        u_alp=linear_fuc(np.log10([7,45]),slope,intercept)
        u_nlm=linear_fuc(np.log10([7,130]),slope,intercept)
        u_app=linear_fuc(np.log10([14]),slope,intercept)
        
        #ax0.scatter([10,45],10**u_alp,s=50,color="forestgreen")

        v,u2,ea1,ea2,slope2,intercept2=linear_internal_regression(a,c,np.min(a),np.max(a))
        ax1.plot(10**v,10**u2,color = 'forestgreen',linewidth=1.2)





        ax0.legend(loc='upper center', bbox_to_anchor=(1.1, 1.30),ncol=3, columnspacing=0.02,handletextpad=0.005, shadow=True,fontsize=8)
        ax0.tick_params(axis="y",direction="in")
        ax1.tick_params(axis="y",direction="in")
        ax0.tick_params(axis="x",direction="in")
        ax1.tick_params(axis="x",direction="in")

        ax0.tick_params(left=True,right=True,labelbottom=True) 
        ax0.set_ylim(0.01,100)
        ax1.set_ylim(0.01,100)
        ax0.set_xlim(3,6000)
        ax1.set_xlim(3,6000)


        plt.setp(ax0.spines.values(), linewidth=1.4)
        plt.setp(ax1.spines.values(), linewidth=1.4)
        ax1.tick_params(left=True,right=True,labelbottom=True,labelleft = False) 

    

        ax0.tick_params(width=1.2)
    

        ax0.set_xlabel(r'$v_{\mathrm{tearing}}$ [$\frac{\mathrm{cm}}{\mathrm{yr}}$]',fontsize=fnt_g.label_)
        ax0.set_ylabel(r'$\dot{H}_{\mathrm{mean}}$ [$\frac{\mathrm{mm}}{\mathrm{yr}}$]',fontsize=fnt_g.label_)
        ax1.set_xlabel(r'$v_{\mathrm{tearing}}$ [$\frac{\mathrm{cm}}{\mathrm{yr}}$]',fontsize=fnt_g.label_)


        ax0.xaxis.set_tick_params(labelsize=fnt_g.axis_)
        ax0.yaxis.set_tick_params(labelsize=fnt_g.axis_)
        
        
        ax1.xaxis.set_tick_params(labelsize=fnt_g.axis_)
        ax1.yaxis.set_tick_params(labelsize=fnt_g.axis_)

        ax0.set_xscale('log')
        ax1.set_xscale('log')
        props = dict(boxstyle='round', facecolor='k',edgecolor='none', alpha=0.8)
        props2 = dict(boxstyle='round', facecolor='w',edgecolor='none', alpha=0.8)

        ax0.text(0.05, 0.96, '$[a]$', transform=ax0.transAxes, fontsize=fnt_g.label_,
            verticalalignment='top', bbox=props,color='white')
        ax1.text(0.05, 0.96, '$[b]$', transform=ax1.transAxes, fontsize=fnt_g.label_,
            verticalalignment='top', bbox=props,color='white')
    
        ax0.text(0.25, 0.1, '$log_{10}(\dot{H})=%.1f log_{10}(v_{tearing})%.1f$'%(slope,intercept), transform=ax0.transAxes, fontsize=9,
            verticalalignment='top', bbox=props,color='k')
        ax1.text(0.25, 0.1, '$log_{10}(\dot{H})=%.1f log_{10}(v_{tearing})%.1f$'%(slope2,intercept2), transform=ax1.transAxes, fontsize=9,
            verticalalignment='top', bbox=props2,color='k')



        if ip == 0:
            ax0.set_yscale('log')
            ax1.set_yscale('log')
        

        fg.savefig(fn,dpi=600)
        plt.close()
        ax0 = []
        
def figure_S10(DB,path_save,figure_name):
    """
    Major update: better to put everything in a row otherwise it is horrible in a paper
    
    """
    # figure name
    fn = os.path.join(path_save,'%s.png'%(figure_name))
    
    # Prepare variables
    depth_tearing = (DB.depth_tearing)
    AVol        = DB.Avolume 
    T           = DB.Temp 
    SLim        = DB.StressLimit/1e6
    
    # Prepare axis of the figures 
    
    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(18*cm, 9*cm))  
    bx = 0.12
    by = 0.15
    sx = 0.27
    dx = 0.02
    sy = 0.69
    dy = 0.01
    ax0 = fg.add_axes([bx, by, sx, sy])
    ax1 = fg.add_axes([bx+dx+sx, by, sx, sy])
    ax2 = fg.add_axes([bx+2*dx+2*sx, by, sx, sy]) 
    
    colors = ['royalblue','goldenrod','tomato']
    label_fig = [r'$v_s = 10$ [$\frac{\mathrm{cm}}{\mathrm{yr}}$]',r'$v_s = 5.0$ [$\frac{\mathrm{cm}}{\mathrm{yr}}$]',r'$v_s = 2.5$ [$\frac{\mathrm{cm}}{\mathrm{yr}}$]']

    T_u = np.sort(np.unique(T))
    T_u = T_u[T_u>0.0]
    
    # find the global limit of the axis:
    min_ax = -180.0
    max_ax = -100.0
    
    for i in range(len(T_u)):
        ax0.scatter(AVol[(SLim==200.0) & (T == T_u[i])],depth_tearing[(SLim==200.0) & (T == T_u[i]),0],s=50,c=colors[i],edgecolor = 'k',label=label_fig[i])
    
    for i in range(len(T_u)):
        ax1.scatter(AVol[(SLim==400.0) & (T == T_u[i])],depth_tearing[(SLim==400.0) & (T == T_u[i]),0],s=50,c=colors[i],edgecolor = 'k',label=label_fig[i])
    
    for i in range(len(T_u)):
        ax2.scatter(AVol[(SLim==600.0) & (T == T_u[i])],depth_tearing[(SLim==600.0) & (T == T_u[i]),0],s=50,c=colors[i],edgecolor = 'k',label=label_fig[i])
    
    ax1.legend(loc='upper center', bbox_to_anchor=(0.45, 1.20),ncol=3, columnspacing=0.10,handletextpad=0.01, shadow=True,fontsize=8)
    
    #ax0.set_ytick(0.1,0.5,0.9)
    
    ax0.set_ylim(min_ax,max_ax)
    ax1.set_ylim(min_ax,max_ax)
    ax2.set_ylim(min_ax,max_ax)
    
    ax0.tick_params(axis="y",direction="in")
    ax0.tick_params(axis="x",direction="in")
    ax0.tick_params(left=True,right=True,labelbottom=True,labelleft = True) 
    ax1.tick_params(left=True,right=True,labelbottom=True,labelleft = False) 
    ax2.tick_params(left=True,right=True,labelbottom=True,labelleft = False) 

    ax1.tick_params(axis="y",direction="in")
    ax1.tick_params(axis="x",direction="in")
    ax2.tick_params(axis="y",direction="in")
    ax2.tick_params(axis="x",direction="in")
    
    plt.setp(ax0.spines.values(), linewidth=1.4)
    plt.setp(ax1.spines.values(), linewidth=1.4)
    plt.setp(ax2.spines.values(), linewidth=1.4)


    ax0.tick_params(width=1.2)
    ax1.tick_params(width=1.2)
    ax2.tick_params(width=1.2)
    
    ax0.set_ylabel(r'$d_{tearing}$ [km]',fontsize=fnt_g.label_)

    ax0.set_xlabel(r'$V_{a,dis}$ [$\mu \frac{\mathrm{m}^3}{\mathrm{mol}}$]',fontsize=fnt_g.label_)
    ax1.set_xlabel(r'$V_{a,dis}$ [$\mu \frac{\mathrm{m}^3}{\mathrm{mol}}$]',fontsize=fnt_g.label_)
    ax2.set_xlabel(r'$V_{a,dis}$ [$\mu \frac{\mathrm{m}^3}{\mathrm{mol}}$]',fontsize=fnt_g.label_)

    
    ax0.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax0.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax2.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax2.yaxis.set_tick_params(labelsize=fnt_g.axis_)

    #ax0.set_yscale('log')
    #ax1.set_yscale('log')
    #ax2.set_yscale('log')


    
    props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)
    ax0.text(0.87, 0.94, '$[a]$', transform=ax0.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.87, 0.94, '$[b]$', transform=ax1.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    ax2.text(0.87, 0.94, '$[c]$', transform=ax2.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    
    ax0.text(0.05, 0.10, r'$\tau_{lim} = 200 [MPa]$', transform=ax0.transAxes, fontsize=fnt_g.axis_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.05, 0.10, r'$\tau_{lim} = 400 [MPa]$', transform=ax1.transAxes, fontsize=fnt_g.axis_,
        verticalalignment='top', bbox=props,color='white')
    ax2.text(0.05, 0.10, r'$\tau_{lim} = 600 [MPa]$', transform=ax2.transAxes, fontsize=fnt_g.axis_,
        verticalalignment='top', bbox=props,color='white')
    

    
    fg.savefig(fn,dpi=600)
    
def make_figure9_sup_Depth(DB,path_save,figure_name):
    
    # figure name
    
    # Prepare variables
    depth_tearing = (DB.depth_tearing)
    AVol        = DB.Avolume 
    T           = DB.Temp 
    SLim        = DB.StressLimit/1e6
    
    
    # Prepare axis of the figures 
    
    cm = 1/2.54  # centimeters in inches
     
    bx = 0.15
    by = 0.15
    sx = 0.40
    dx = 0.02
    sy = 0.6
    dy = 0.01
    type = ['a','b','c']
    Uplift_discrete = DB.Uplift_Te_discrete
    for ip in range(3):
        fg = figure(figsize=(18*cm, 9*cm)) 
        ax0 = fg.add_axes([bx, by, sx, sy])
        ax1 = fg.add_axes([bx+sx+dx,by,sx,sy])

        fn = os.path.join(path_save,'%s%s.png'%(figure_name,type[ip]))

        uplift = DB.uplift[:,ip]
        colors = ['royalblue','goldenrod','tomato']#,'orange','grey','pink']
        label_fig = [r'$v_s = 10$ [$\frac{\mathrm{cm}}{\mathrm{yr}}$]',r'$v_s = 5.0$ [$\frac{\mathrm{cm}}{\mathrm{yr}}$]',r'$v_s = 2.5$ [$\frac{\mathrm{cm}}{\mathrm{yr}}$]']

        T_u = np.sort(np.unique(T))
        T_u = T_u[T_u>0.0]
        for i in range(len(T_u)):
            ax0.scatter(depth_tearing[(T == T_u[i]),0],uplift[(T == T_u[i])],c=colors[i],s=50,edgecolor = 'k',label=label_fig[i])
            ax1.scatter(depth_tearing[(T == T_u[i]),0],Uplift_discrete[(T == T_u[i])],c=colors[i],s=50,edgecolor = 'k',label=label_fig[i])
        #ax0.axvline(2,linewidth=0.8,color='k',alpha=0.5)
        #ax0.axvline(94,linewidth=0.8,color='k',alpha=0.5)
        #ax1.axvline(2,linewidth=0.8,color='k',alpha=0.5)
        #ax1.axvline(94,linewidth=0.8,color='k',alpha=0.5)


        ax0.legend(loc='upper center', bbox_to_anchor=(1.1, 1.30),ncol=3, columnspacing=0.02,handletextpad=0.005, shadow=True,fontsize=8)
        ax0.tick_params(axis="y",direction="in")
        ax1.tick_params(axis="y",direction="in")
        ax0.tick_params(axis="x",direction="in")
        ax1.tick_params(axis="x",direction="in")

        ax0.tick_params(left=True,right=True,labelbottom=True) 
        ax0.set_ylim(0.01,100)
        ax1.set_ylim(0.01,100)


        plt.setp(ax0.spines.values(), linewidth=1.4)
        plt.setp(ax1.spines.values(), linewidth=1.4)
        ax1.tick_params(left=True,right=True,labelbottom=True,labelleft = False) 

    

        ax0.tick_params(width=1.2)
    

        ax0.set_xlabel(r'$d_{\mathrm{tearing}}$ [$\frac{\mathrm{cm}}{\mathrm{yr}}$]',fontsize=fnt_g.label_)
        ax0.set_ylabel(r'$\dot{H}_{\mathrm{mean}}$ [$\frac{\mathrm{mm}}{\mathrm{yr}}$]',fontsize=fnt_g.label_)
        ax1.set_xlabel(r'$d_{\mathrm{tearing}}$ [$\frac{\mathrm{cm}}{\mathrm{yr}}$]',fontsize=fnt_g.label_)


        ax0.xaxis.set_tick_params(labelsize=fnt_g.axis_)
        ax0.yaxis.set_tick_params(labelsize=fnt_g.axis_)
        
        
        ax1.xaxis.set_tick_params(labelsize=fnt_g.axis_)
        ax1.yaxis.set_tick_params(labelsize=fnt_g.axis_)

        #ax0.set_xscale('log')
        #ax1.set_xscale('log')
        props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)
        ax0.text(0.05, 0.96, '$[a]$', transform=ax0.transAxes, fontsize=fnt_g.label_,
            verticalalignment='top', bbox=props,color='white')
        ax1.text(0.05, 0.96, '$[b]$', transform=ax1.transAxes, fontsize=fnt_g.label_,
            verticalalignment='top', bbox=props,color='white')
    

        if ip == 0:
            ax0.set_yscale('log')
            ax1.set_yscale('log')

            print('Hell ya')
        

        fg.savefig(fn,dpi=600)
        plt.close()
        ax0 = []


def figure_S8(DB,path_save,figure_name):
    """
    Major update: better to put everything in a row otherwise it is horrible in a paper
    This figure represents the summary of the results as a function of vs,Va 
    """
    # figure name
    fn = os.path.join(path_save,'%s.png'%(figure_name))
    
    # Prepare variables
    vel_tearing = (DB.detachment_velocity)
    AVol        = DB.Avolume 
    T           = DB.Temp 
    SLim        = DB.StressLimit/1e6
    
    # Prepare axis of the figures 
    
    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(18*cm, 9*cm))  
    bx = 0.12
    by = 0.15
    sx = 0.27
    dx = 0.02
    sy = 0.50
    dy = 0.02
    ax0 = fg.add_axes([bx, by, sx, sy])
    ax1 = fg.add_axes([bx+dx+sx, by, sx, sy])
    ax2 = fg.add_axes([bx+2*dx+2*sx, by, sx, sy])
    ax3 = fg.add_axes([bx+0.5*sx, by+sy+dy, 2.0*sx+2*dx, 0.28])
    
    
    colors = ['royalblue','goldenrod','tomato']

    T_u = np.sort(np.unique(T))
    T_u = T_u[T_u>0.0]
    
    # find the global limit of the axis:
    #min_ax = 1.0
    #max_ax = np.nanmax(vel_tearing[(T > 0)])
    #max_ax = max_ax+max_ax*0.2

    cf1=ax0.scatter(AVol[(SLim==200.0)],-T[(SLim==200.0)],s=50,c=np.log10(vel_tearing[(SLim==200.0)]),cmap = 'cmc.lipari',edgecolor = 'k')
    cf2=ax1.scatter(AVol[(SLim==400.0)],-T[(SLim==400.0)],s=50,c=np.log10(vel_tearing[(SLim==400.0)]),cmap = 'cmc.lipari',edgecolor = 'k')
    cf3=ax2.scatter(AVol[(SLim==600.0)],-T[(SLim==600.0)],s=50,c=np.log10(vel_tearing[(SLim==600.0)]),cmap = 'cmc.lipari',edgecolor = 'k')
    
    
    
    #ax0.set_ytick(0.1,0.5,0.9)
    
    
    ax0.tick_params(axis="y",direction="in")
    ax0.tick_params(axis="x",direction="in")
    ax0.tick_params(left=True,right=True,labelbottom=True,labelleft = True) 
    ax1.tick_params(left=True,right=True,labelbottom=True,labelleft = False) 
    ax2.tick_params(left=True,right=True,labelbottom=True,labelleft = False) 

    ax1.tick_params(axis="y",direction="in")
    ax1.tick_params(axis="x",direction="in")
    ax2.tick_params(axis="y",direction="in")
    ax2.tick_params(axis="x",direction="in")
    
    plt.setp(ax0.spines.values(), linewidth=1.4)
    plt.setp(ax1.spines.values(), linewidth=1.4)
    plt.setp(ax2.spines.values(), linewidth=1.4)


    ax0.tick_params(width=1.2)
    ax1.tick_params(width=1.2)
    ax2.tick_params(width=1.2)
    
    ax0.set_ylabel(r'$v_{s}$ [$\frac{\mathrm{cm}}{\mathrm{yr}}$]',fontsize=fnt_g.label_)

    ax0.set_xlabel(r'$V_{a,dis}$ [$\frac{\mathrm{cm}^3}{\mathrm{mol}}$]',fontsize=fnt_g.label_)
    ax1.set_xlabel(r'$V_{a,dis}$ [$ \frac{\mathrm{cm}^3}{\mathrm{mol}}$]',fontsize=fnt_g.label_)
    ax2.set_xlabel(r'$V_{a,dis}$ [$ \frac{\mathrm{cm}^3}{\mathrm{mol}}$]',fontsize=fnt_g.label_)

    
    ax0.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax0.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax2.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax2.yaxis.set_tick_params(labelsize=fnt_g.axis_)




    lim_def_0 = np.round(np.log10([np.min(vel_tearing[T>0]),np.max(vel_tearing[T>0])]))
    lim_  = np.linspace(lim_def_0[0],lim_def_0[1],num=4)
    cbaxes,cbar = define_colorbar(cf2,ax3,lim_def_0,(lim_),r'$\mathrm{log}_{10}\left(v_{\mathrm{tearing}}\right)$')
    ax3.axis('off')
    v1 = np.max(T[T>0])
    v3 = np.min(T[T>0])
    v2 = np.unique(T[(T!=v1)&(T!=v3)&(T>0)])
    v2 = v2[0]
    ax0.set_yticks([-v1,-v2,-v3])
    ax1.set_yticks([-v1,-v2,-v3])
    ax2.set_yticks([-v1,-v2,-v3])

    ax0.set_ylim(-v1-100,-v2+100)
    ax1.set_ylim(-v1-100,-v2+100)
    ax2.set_ylim(-v1-100,-v2+100)


    l_abelr = ['$10.0$','$5.0$','$2.5$']
    ax0.set_yticklabels(l_abelr)
    
    props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)
    ax0.text(0.87, 0.94, '$[a]$', transform=ax0.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.87, 0.94, '$[b]$', transform=ax1.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    ax2.text(0.87, 0.94, '$[c]$', transform=ax2.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    
    ax0.text(0.05, 0.13, r'$\tau_{lim} = 200 [MPa]$', transform=ax0.transAxes, fontsize=fnt_g.axis_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.05, 0.13, r'$\tau_{lim} = 400 [MPa]$', transform=ax1.transAxes, fontsize=fnt_g.axis_,
        verticalalignment='top', bbox=props,color='white')
    ax2.text(0.05, 0.13, r'$\tau_{lim} = 600 [MPa]$', transform=ax2.transAxes, fontsize=fnt_g.axis_,
        verticalalignment='top', bbox=props,color='white')
    

    
    fg.savefig(fn,dpi=600)


def figure_S11(T1:Test,T2:Test,T3:Test,T4:Test,T5:Test,T6:Test,path_save:str,figure_name:str):
    
    def _extracte_average_uplift(time,t,Starting_tearing,Ending_tearing,T1):
        uplift = np.zeros([3,len(t)],dtype=float)

        ip = 0
        for i in range(len(time)):
            if (time[i]>= Starting_tearing) & (time[i]<=Ending_tearing): 
                dH_trench = _interpolate_2D(T1.FS.dH[:,:,i],T1.C.xg,T1.C.yg,T1.C.x_trench_p,T1.C.y_trench_p)
                uplift[0,ip] = np.nanmax(dH_trench)
                uplift[1,ip] = np.nanmean(dH_trench)
                uplift[2,ip] = np.nanmin(dH_trench)
                ip +=1 
    
        return uplift 

    
    def _create_plot_(T,ax,flag_axis,letter):
        
        ax.spines['bottom'].set_color('k')
        ax.spines['top'].set_color('w') 
        ax.spines['right'].set_color('k')
        ax.spines['left'].set_color('k')

        ax.tick_params(axis="y",direction="in")
        ax.tick_params(axis="x",direction="in")
        ax.tick_params(width=1.2)
        ax.xaxis.set_tick_params(labelsize=fnt_g.axis_)
        ax.yaxis.set_tick_params(labelsize=fnt_g.axis_)
        det_vec = T.Det.det_vec
                    # Filter the necking along the width

        det_vec[T.Det.depth_vec>=-60.0]=np.nan


        #Data
        time = T.time 
        tear_vel = T.Det.vel_tear[0,:]
        Starting_tearing  = np.nanmin(T.Det.det_vec[(T.C.x_sp>=100) & (T.C.x_sp<=1100) ])
        Ending_tearing    = np.nanmax(T.Det.det_vec[(T.C.x_sp>=100) & (T.C.x_sp<=1100) ])
        t = time[(time>=Starting_tearing)&(time<=Ending_tearing)]
        uplift = _extracte_average_uplift(time,t,Starting_tearing,Ending_tearing,T)
        mean_v = np.mean(tear_vel[(time>=Starting_tearing)&(time<=Ending_tearing)])

        
        # Main plot
        ax.plot(t,tear_vel[(time>=Starting_tearing)&(time<=Ending_tearing)],linewidth=1.2,color='firebrick')
        ax.axhline(mean_v,linestyle = '-.',linewidth = 0.8, color='firebrick')
        print(mean_v)
        axb = ax.twinx()
        axb.fill_between(t,uplift[2,:],uplift[0,:],color='blue',alpha=0.1,linewidth=1.2)
       
       
        ax.set_ylim([0,np.round(np.nanmax(tear_vel))])
        ax.set_yticks([0.0,np.round(np.nanmax(tear_vel))])
  
        ax.set_ylabel(r'$v_{\mathrm{tearing}}$ $[\frac{\mathrm{cm}}{\mathrm{yr}}]$',fontsize=fnt_g.label_)
        axb.set_ylabel(r'$\dot{H}$ $[\frac{\mathrm{mm}}{\mathrm{yr}}]$',fontsize=fnt_g.label_)

        ax.set_xlim([min(t),max(t)])
        ax.set_xticks([min(t),max(t)])
        if flag_axis==True:
            ax.set_xlabel('Tearing')
            ax.set_xticklabels(['Start','End'])
        else: 
            ax.set_xticklabels([])


        props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)
        ax.text(0.94, 0.94, letter, transform=ax.transAxes, fontsize=fnt_g.label_,
            verticalalignment='top', bbox=props,color='white')
    
        ax.text(0.12, 1.30, r'$time = %.2f-%.2f$ [Myr]'%(Starting_tearing,Ending_tearing), transform=ax.transAxes, fontsize=fnt_g.axis_,verticalalignment='top', bbox=props,color='white')
        
        return ax,axb

    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(19*cm, 18*cm))  
    bx = 0.15
    by = 0.10
    sx = 0.70
    dx = 0.02
    sy = 0.10
    dy = 0.05

    fn = os.path.join(path_save,'%s.png'%(figure_name))

    # Prepare axis 
    ax0 = fg.add_axes([bx, by, sx, sy])
    ax1 = fg.add_axes([bx, by+dy+sy, sx, sy])
    ax2 = fg.add_axes([bx, by+2*(dy+sy), sx, sy]) 
    ax3 = fg.add_axes([bx, by+3*(dy+sy), sx, sy])
    ax4 = fg.add_axes([bx, by+4*(dy+sy), sx, sy])
    ax5 = fg.add_axes([bx, by+5*(dy+sy), sx, sy])
    # Create Axis 
    ax0,ax0tw = _create_plot_(T6,ax0,True,'$[f]$')
    ax1,ax1tw = _create_plot_(T5,ax1,False,'$[e]$')
    ax2,ax2tw = _create_plot_(T4,ax2,False,'$[d]$')
    ax3,ax3tw = _create_plot_(T3,ax3,False,'$[c]$')
    ax4,ax4tw = _create_plot_(T2,ax4,False,'$[b]$')
    ax5,ax5tw = _create_plot_(T1,ax5,False,'$[a]$')

    fg.savefig(fn,dpi=600)


def _gif_topography(A,path_save):
    
    def fmt(x):
        s = f"{x:.1f}"
        if s.endswith("0"):
            s = f"{x:.0f}"
        return rf"{s} km" if plt.rcParams["text.usetex"] else f"{s} km"

    def define_colorbar_S(cf,ax,lim:list,ticks:list,label:str):
        
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    
        cbaxes = inset_axes(ax, borderpad=1.3,  width="100%", height="15%", loc=3)   
        cbar=plt.colorbar(cf,cax=cbaxes, ticks=ticks, orientation='horizontal',extend="both")
        print(lim[0])
        print(lim[1])
        cbar.set_label(label=label,size=fnt_g.label_) 
        cbar.ax.xaxis.set_ticks_position('top')
        cbar.ax.xaxis.set_label_coords(0.5,-0.15)
        cbar.ax.xaxis.set_tick_params(pad=0.1)
        cbar.ax.xaxis.set_label_position('bottom')
        cbar.ax.tick_params(labelsize=fnt_g.axis_)
    
        return cbaxes,cbar
    
    def _compute_velocity_migration(A):

        vel    = np.zeros([len(A.time),1],dtype=float)
        x_depo = np.zeros([len(A.time),1],dtype=float)
        y_depo = np.zeros([len(A.time),1],dtype=float)
        dH     = np.zeros([len(A.time),3],dtype=float)


        for ipic in range(len(A.time)):
            if A.time[ipic]>1:
                H_trench = _interpolate_2D(A.FS.H[:,:,ipic],A.C.xg,A.C.yg,A.C.x_trench_p,A.C.y_trench_p)
                dH_trench = _interpolate_2D(A.FS.dH_fil[:,:,ipic],A.C.xg,A.C.yg,A.C.x_trench_p,A.C.y_trench_p)

                Depo = np.where((H_trench<=0) & (A.C.x_sp>=100))
                if len(Depo[0])>0:
                    ind = Depo[0][0]
                    tt = A.time[ipic]
                    x_depo[ipic] = A.C.x_trench_p[ind]
                    y_depo[ipic] = A.C.y_trench_p[ind]
                    d=((x_depo[ipic]-x_depo[ipic-1])**2+(y_depo[ipic]-y_depo[ipic-1])**2)**0.5
                    d *= 1e3*1e2
                    v_sig = np.sign(x_depo[ipic]-x_depo[ipic-1])
                    dt = (A.time[ipic]-A.time[ipic-1])*1e6
                    vel[ipic] = v_sig*d/dt
                    dH[ipic,0]   = np.nanmin(dH_trench)
                    dH[ipic,1] = np.nanmean(dH_trench)
                    dH[ipic,2] = np.nanmax(dH_trench)
                    print('v depocenter is %.2f'%vel[ipic])
        
        return x_depo, y_depo,vel,dH



    x_depo,y_depo,vel,dH = _compute_velocity_migration(A)


    path_saveb = os.path.join(path_save,'FM_def',A.Test_Name)
    if os.path.isdir(path_saveb)==False:
        os.mkdir(path_saveb)

 
    check = 0
    for ipic in range(len(A.time)):
        if A.time[ipic]>1.0:
            if check ==0:
                itf = ipic
                check = 1

            figure_name = 'Figure_TopoGif_%d'%ipic
            fg = figure()        
            cm = 1/2.54  # centimeters in inches
            fg = figure(figsize=(14*cm, 14*cm))  
            bx = 0.15
            by = 0.10
            sx = 0.8
            sx2 = 0.7
            sxc = 0.6
            sy1 = 0.20
            sy2 = 0.14
            sy3 = 0.38
            dx = 0.04
            dy = 0.08
            dy2 = 0.05
            fn = os.path.join(path_save,'%s.png'%(figure_name))
            # Prepare axis 
            ax0 = fg.add_axes([bx, by, sx2, sy1])
            ax1 = fg.add_axes([bx,by+sy1+dy,sxc,sy2])
            ax2 = fg.add_axes([bx, by+sy1+sy2+2*dy, sx, sy3])
            

            fn = os.path.join(path_saveb,figure_name)
            p1 = ax2.contourf(A.C.xg,A.C.yg,A.FS.H[:,:,ipic],levels = [-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5],cmap = 'cmc.oleron')#linewidths=0.5)
            ax2.plot(A.C.x_trench_p,A.C.y_trench_p,linewidth = 2.0,linestyle = 'dashdot',label = r'Slab position',color = 'firebrick')
            ax2.scatter(x_depo[ipic],y_depo[ipic],s=50,marker='D',c = 'firebrick')
            ax2.set_xlabel(r'$x$ [km]',fontsize=fnt_g.label_)
            ax2.set_ylabel(r'$y$ [km]',fontsize=fnt_g.label_)
            ax2.legend(loc='upper right')
            ax2.xaxis.set_tick_params(labelsize=fnt_g.axis_)
            ax2.yaxis.set_tick_params(labelsize=fnt_g.axis_)
            
            cbaxes,cbar =define_colorbar_S(p1,ax1,[-4.5,4.5],[-4.0,-2.0,0.0,2.0,4.0],r'${{H}}$ [$\mathrm{km}$]')
            ax1.axis('off')
            ax0.plot(range(itf,ipic),vel[itf:ipic],linewidth=1.5,color='forestgreen')
            ax0.set_xlabel(r'$time$ [Myr]',fontsize=fnt_g.label_)
            ax0.set_ylabel(r'$v_{mig}$ $[\frac{\mathrm{cm}}{\mathrm{yr}}]$',fontsize=fnt_g.label_)
            
            Starting_tearing  = np.nanmin(A.Det.det_vec[(A.C.x_sp>=100) & (A.C.x_sp<=1100) ])
            Ending_tearing  = np.nanmax(A.Det.det_vec[(A.C.x_sp>=100) & (A.C.x_sp<=1100) ])
            indTb = np.where(A.time==Starting_tearing)
            indTe = np.where(A.time==Ending_tearing)
            indTb = indTb[0][0]
            indTe = indTe[0][0]

            ax0.spines['bottom'].set_color('k')
            ax0.spines['top'].set_color('w') 
            ax0.spines['right'].set_color('k')
            ax0.spines['left'].set_color('k')
            ax0.set_xlim(np.round(indTb)-10,len(A.time))
            ax0.set_ylim([0,np.max(vel)])
            ax0.axvspan(indTb, indTe, alpha=0.1, color='red')
            axb = ax0.twinx()
            axb.plot(range(itf,ipic),(dH[itf:ipic,1]),linewidth=0.8,linestyle=':',color='blue')
            axb.set_ylim([np.nanmin(dH[:,0]),np.nanmax(dH[:,2])])
            axb.fill_between(range(itf,ipic),(dH[itf:ipic,0]),(dH[itf:ipic,2]),linewidth=0.2,color='blue', alpha = 0.3)

            axb.set_ylabel(r'$\dot{\bar{H}}$ $[\frac{\mathrm{mm}}{\mathrm{yr}}]$',fontsize=fnt_g.label_)
            ax0.set_xlim(np.round(indTb)-10,len(A.time))
            ax0.set_xticks([np.round(indTb)-10,np.round(indTb)-5,np.round(indTb),np.round(indTe),np.round(indTe)+5])
            ax0.set_xlabel([np.round(A.time[np.round(indTb)-10],3),np.round(A.time[np.round(indTb)-5],3),np.round(A.time[np.round(indTb)],3),np.round(A.time[np.round(indTe)],3)
                            ,np.round(A.time[np.round(indTe)],3)])

            axb.set_xticks([np.round(indTb)-10,np.round(indTb)-5,np.round(indTb),np.round(indTe),np.round(indTe)+5])

            axb.set_xlabel([np.round(A.time[np.round(indTb)-10],3),np.round(A.time[np.round(indTb)-5],3),np.round(A.time[np.round(indTb)],3),np.round(A.time[np.round(indTe)],3)
                            ,np.round(A.time[np.round(indTe)],3)])



            ax0.tick_params(axis="y",direction="in")
            ax0.tick_params(axis="x",direction="in")
            ax0.tick_params(width=1.2)


            ax2.tick_params(axis="y",direction="in")
            ax2.tick_params(axis="x",direction="in")
            ax2.tick_params(width=1.2)


            
            props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)

            ax2.text(0.05, 0.95, r'$t= %.2f Myr$ $v_{mig} = %.2f [\frac{\mathrm{cm}}{\mathrm{yr}}]$'%(A.time[ipic],vel[ipic]), transform=ax2.transAxes, fontsize=fnt_g.axis_,verticalalignment='top', bbox=props,color='white')

            fg.savefig(fn,dpi=600,transparent=False)



def figure_S12(T1:Test,T2:Test,T3:Test,T4:Test,T5:Test,T6:Test,path_save:str,figure_name:str):
    
    def _compute_velocity_migration(T):

        vel    = np.zeros([len(T.time),1],dtype=float)
        x_depo = np.zeros([len(T.time),1],dtype=float)
        y_depo = np.zeros([len(T.time),1],dtype=float)
        dH     = np.zeros([len(T.time),3],dtype=float)


        for ipic in range(len(T.time)):
            if T.time[ipic]>1:
                H_trench = _interpolate_2D(T.FS.H[:,:,ipic],T.C.xg,T.C.yg,T.C.x_trench_p,T.C.y_trench_p)
                dH_trench = _interpolate_2D(T.FS.dH_fil[:,:,ipic],T.C.xg,T.C.yg,T.C.x_trench_p,T.C.y_trench_p)

                Depo = np.where((H_trench<=0) & (T.C.x_sp>=100))
                if len(Depo[0])>0:
                    ind = Depo[0][0]
                    tt = T.time[ipic]
                    x_depo[ipic] = T.C.x_trench_p[ind]
                    y_depo[ipic] = T.C.y_trench_p[ind]
                    d=((x_depo[ipic]-x_depo[ipic-1])**2+(y_depo[ipic]-y_depo[ipic-1])**2)**0.5
                    d *= 1e3*1e2
                    v_sig = np.sign(x_depo[ipic]-x_depo[ipic-1])
                    dt = (T.time[ipic]-T.time[ipic-1])*1e6
                    vel[ipic] = v_sig*d/dt
                    dH[ipic,0]   = np.nanmin(dH_trench)
                    dH[ipic,1] = np.nanmean(dH_trench)
                    dH[ipic,2] = np.nanmax(dH_trench)
                    print('v depocenter is %.2f'%vel[ipic])
        
        return x_depo, y_depo,vel,dH





    
    def _create_plot_(T,ax,flag_axis,letter):

        x_depo,y_depo,vel,dH = _compute_velocity_migration(T)

        ipic = len(T.time)
        itf = np.where(T.time>=1.0)
        itf = itf[0][0]
        
        ax.spines['bottom'].set_color('k')
        ax.spines['top'].set_color('w') 
        ax.spines['right'].set_color('k')
        ax.spines['left'].set_color('k')

        ax.tick_params(axis="y",direction="in")
        ax.tick_params(axis="x",direction="in")
        ax.tick_params(width=1.2)
        ax.xaxis.set_tick_params(labelsize=fnt_g.axis_,pad=7)
        ax.yaxis.set_tick_params(labelsize=fnt_g.axis_)


        ax.plot(range(itf,ipic),vel[itf:ipic],linewidth=1.5,color='forestgreen')
        ax.set_ylabel(r'$v_{mig}$ $[\frac{\mathrm{cm}}{\mathrm{yr}}]$',fontsize=fnt_g.label_)
            
        Starting_tearing  = np.nanmin(T.Det.det_vec[(T.C.x_sp>=100) & (T.C.x_sp<=1100) ])
        Ending_tearing  = np.nanmax(T.Det.det_vec[(T.C.x_sp>=100) & (T.C.x_sp<=1100) ])
        indTb = np.where(T.time==Starting_tearing)
        indTe = np.where(T.time==Ending_tearing)
        indTb = indTb[0][0]
        indTe = indTe[0][0]

        ax.spines['bottom'].set_color('k')
        ax.spines['top'].set_color('w') 
        ax.spines['right'].set_color('k')
        ax.spines['left'].set_color('k')
        ax.set_xlim(np.round(indTb)-5,np.round(indTe)+1)

        ax.set_ylim([0,np.max(vel)])
        ax.axvspan(indTb, indTe, alpha=0.1, color='red')
        axb = ax.twinx()
        axb.plot(range(itf,ipic),(dH[itf:ipic,1]),linewidth=0.8,linestyle=':',color='blue')
        axb.set_ylim([np.nanmin(dH[:,0]),np.nanmax(dH[:,2])])
        axb.fill_between(range(itf,ipic),(dH[itf:ipic,0]),(dH[itf:ipic,2]),linewidth=0.2,color='blue', alpha = 0.3)

        axb.set_ylabel(r'$\dot{\bar{H}}$ $[\frac{\mathrm{mm}}{\mathrm{yr}}]$',fontsize=fnt_g.label_)
        ax.set_xlim(np.round(indTb)-10,len(T.time))
        ax.set_xticks([np.round(indTb)-10,np.round(indTb)-5,np.round(indTb),np.round(indTe)])
        
        labels_x = np.round([T.time[np.round(indTb)-10],T.time[np.round(indTb)-5],T.time[np.round(indTb)],T.time[np.round(indTe)]],3)
        axb.set_xticklabels([r'$%.1f$'%(labels_x[0]),r'$%.1f$'%(labels_x[1]),r'$%.1f$'%(labels_x[2]),r'$%.1f$'%(labels_x[3])])
        axb.set_xlim(np.round(indTb)-5,np.round(indTe)+1)

        axb.set_xlabel([])



        ax.tick_params(axis="y",direction="in")
        ax.tick_params(axis="x",direction="in")
        ax.tick_params(width=1.2)


        axb.xaxis.set_tick_params(labelsize=fnt_g.axis_)
        axb.yaxis.set_tick_params(labelsize=fnt_g.axis_)

        props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)

        ax.text(0.01, 0.95,letter, transform=ax.transAxes, fontsize=fnt_g.axis_,verticalalignment='top', bbox=props,color='white')


        if flag_axis == True: 
            ax.set_xlabel(r'$time$ [Myr]')


        
        return ax,axb

    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(19*cm, 18*cm))  
    bx = 0.15
    by = 0.10
    sx = 0.70
    dx = 0.02
    sy = 0.10
    dy = 0.05
    path_save = os.path.join(path_save,'FM_def')

    fn = os.path.join(path_save,'%s.png'%(figure_name))

    # Prepare axis 
    ax0 = fg.add_axes([bx, by, sx, sy])
    ax1 = fg.add_axes([bx, by+dy+sy, sx, sy])
    ax2 = fg.add_axes([bx, by+2*(dy+sy), sx, sy]) 
    ax3 = fg.add_axes([bx, by+3*(dy+sy), sx, sy])
    ax4 = fg.add_axes([bx, by+4*(dy+sy), sx, sy])
    ax5 = fg.add_axes([bx, by+5*(dy+sy), sx, sy])
    # Create Axis 
    ax0,ax0tw = _create_plot_(T6,ax0,True,'$[f]$')
    ax1,ax1tw = _create_plot_(T5,ax1,False,'$[e]$')
    ax2,ax2tw = _create_plot_(T4,ax2,False,'$[d]$')
    ax3,ax3tw = _create_plot_(T3,ax3,False,'$[c]$')
    ax4,ax4tw = _create_plot_(T2,ax4,False,'$[b]$')
    ax5,ax5tw = _create_plot_(T1,ax5,False,'$[a]$')

    fg.savefig(fn,dpi=600)



 



def figure_10(A,B,tA,tB,path_save,figure_name):
    
    def fmt(x):
        s = f"{x:.1f}"
        if s.endswith("0"):
            s = f"{x:.0f}"
        return rf"{s} km" if plt.rcParams["text.usetex"] else f"{s} km"

    def define_colorbar_S(cf,ax,lim:list,ticks:list,label:str):
        
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    
        cbaxes = inset_axes(ax, borderpad=1.3,  width="100%", height="15%", loc=3)   
        cbar=plt.colorbar(cf,cax=cbaxes, ticks=ticks, orientation='horizontal',extend="both")
        print(lim[0])
        print(lim[1])
        cbar.set_label(label=label,size=fnt_g.label_) 
        cbar.ax.xaxis.set_ticks_position('top')
        cbar.ax.xaxis.set_label_coords(0.5,-0.15)
        cbar.ax.xaxis.set_tick_params(pad=0.1)
        cbar.ax.xaxis.set_label_position('bottom')
        cbar.ax.tick_params(labelsize=fnt_g.axis_)
    
        return cbaxes,cbar
    
    def _compute_velocity_migration(T):

        vel    = np.zeros([len(T.time),1],dtype=float)
        x_depo = np.zeros([len(T.time),1],dtype=float)
        y_depo = np.zeros([len(T.time),1],dtype=float)
        dH     = np.zeros([len(T.time),3],dtype=float)


        for ipic in range(len(T.time)):
            if T.time[ipic]>1:
                H_trench = _interpolate_2D(T.FS.H[:,:,ipic],T.C.xg,T.C.yg,T.C.x_trench_p,T.C.y_trench_p)
                dH_trench = _interpolate_2D(T.FS.dH_fil[:,:,ipic],T.C.xg,T.C.yg,T.C.x_trench_p,T.C.y_trench_p)

                Depo = np.where((H_trench<=0) & (A.C.x_sp>=100))
                if len(Depo[0])>0:
                    ind = Depo[0][0]
                    tt = T.time[ipic]
                    x_depo[ipic] = T.C.x_trench_p[ind]
                    y_depo[ipic] = T.C.y_trench_p[ind]
                    d=((x_depo[ipic]-x_depo[ipic-1])**2+(y_depo[ipic]-y_depo[ipic-1])**2)**0.5
                    d *= 1e3*1e2
                    v_sig = np.sign(x_depo[ipic]-x_depo[ipic-1])
                    dt = (T.time[ipic]-T.time[ipic-1])*1e6
                    vel[ipic] = v_sig*d/dt
                    dH[ipic,0]   = np.nanmin(dH_trench)
                    dH[ipic,1] = np.nanmean(dH_trench)
                    dH[ipic,2] = np.nanmax(dH_trench)
        return x_depo, y_depo,vel,dH


    def _create_pic_typeA(ax,T,tv,letter):
        x_depo,y_depo,vel,dH = _compute_velocity_migration(T)

        ipic = np.where(T.time>=tv)[0][0]

        ax.spines['bottom'].set_color('k')
        ax.spines['top'].set_color('k') 
        ax.spines['right'].set_color('k')
        ax.spines['left'].set_color('k')

        ax.tick_params(axis="y",direction="in")
        ax.tick_params(axis="x",direction="in")
        ax.tick_params(width=1.2)
        ax.xaxis.set_tick_params(labelsize=fnt_g.axis_)
        ax.yaxis.set_tick_params(labelsize=fnt_g.axis_)
        
        pa = ax.contourf(T.C.xg,T.C.yg,T.FS.H[:,:,ipic],levels = [-4.5,-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5],cmap = 'cmc.oleron')#linewidths=0.5)
        pb = ax.contour(T.C.xg,T.C.yg,T.FS.H[:,:,ipic],levels = [0.0],linewidths=1.2,colors='k')
        ax.set_xlim([-650,650])
        ax.set_ylim([-250,100])
        ax.plot(T.C.x_trench_p,T.C.y_trench_p,linewidth = 2.0,linestyle = 'dashdot',label = r'Slab position',color = 'firebrick')
        ax.scatter(x_depo[ipic],y_depo[ipic],s=50,marker='X',c = 'k')
        if (letter == '[c]') or letter == '[g]':
            ax.set_xlabel(r'$x$ [km]',fontsize=fnt_g.label_)
        else: 
            ax.set_xticklabels([])
        
        if (letter == '[a]') or letter == '[b]' or letter == '[c]':
            ax.set_ylabel(r'$y$ [km]',fontsize=fnt_g.label_)
        else:
            ax.set_yticklabels([])

        


        props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)

        ax.text(0.01, 0.95,letter, transform=ax.transAxes, fontsize=fnt_g.axis_,verticalalignment='top', bbox=props,color='white')
        ax.text(0.12, 1.30, r'$time = %.2f$ [Myr]'%(tv), transform=ax.transAxes, fontsize=fnt_g.axis_,verticalalignment='top', bbox=props,color='white')


        return ax,pa

    def _create_pic_typeB(ax,T,letter):
        x_depo,y_depo,vel,dH = _compute_velocity_migration(T)

        ipic = len(T.time)
        itf = np.where(T.time>=1.0)
        itf = itf[0][0]
        
        ax.spines['bottom'].set_color('k')
        ax.spines['top'].set_color('w') 
        ax.spines['right'].set_color('k')
        ax.spines['left'].set_color('k')

        ax.tick_params(axis="y",direction="in")
        ax.tick_params(axis="x",direction="in")
        ax.tick_params(width=1.2)
        ax.xaxis.set_tick_params(labelsize=fnt_g.axis_,pad=7)
        ax.yaxis.set_tick_params(labelsize=fnt_g.axis_)


        ax.plot(range(itf,ipic),vel[itf:ipic],linewidth=1.5,color='forestgreen')
        ax.set_ylabel(r'$v_{mig}$ $[\frac{\mathrm{cm}}{\mathrm{yr}}]$',fontsize=fnt_g.label_)
            
        Starting_tearing  = np.nanmin(T.Det.det_vec[(T.C.x_sp>=100) & (T.C.x_sp<=1100) ])
        Ending_tearing  = np.nanmax(T.Det.det_vec[(T.C.x_sp>=100) & (T.C.x_sp<=1100) ])
        indTb = np.where(T.time==Starting_tearing)
        indTe = np.where(T.time==Ending_tearing)
        indTb = indTb[0][0]
        indTe = indTe[0][0]

        ax.spines['bottom'].set_color('k')
        ax.spines['top'].set_color('w') 
        ax.spines['right'].set_color('k')
        ax.spines['left'].set_color('k')
        ax.set_xlim(np.round(indTb)-5,np.round(indTe)+1)

        ax.set_ylim([0,np.max(vel)])
        ax.axvspan(indTb, indTe, alpha=0.1, color='red')
        ax.set_xlim(np.round(indTb)-10,indTe+2)
        ax.set_xticks([np.round(indTb)-10,np.round(indTb),np.round(indTe)])
        
        labels_x = np.round([T.time[np.round(indTb)-10],T.time[np.round(indTb)],T.time[np.round(indTe)]],3)
        ax.set_xticklabels([r'$%.1f$'%(labels_x[0]),r'$%.1f$'%(labels_x[1]),r'$%.1f$'%(labels_x[2])])
   

        ax.tick_params(axis="y",direction="in")
        ax.tick_params(axis="x",direction="in")
        ax.tick_params(width=1.2)


        props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)

        ax.text(0.01, 0.95,letter, transform=ax.transAxes, fontsize=fnt_g.axis_,verticalalignment='top', bbox=props,color='white')


       
        ax.set_xlabel(r'$time$ [Myr]')


        
        return ax




    fn = os.path.join(path_save,figure_name)
    fg = figure()        
    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(18*cm, 15*cm))  
    bx = 0.10
    by = 0.07
    sx = 0.38  
    sx2 = 0.42
    dx  = 0.05
    dx2 = 0.08
    dy = 0.04 
    sy = 0.15 
    sy2 = 0.13

    ax0 = fg.add_axes([bx, by, sx2, sy2])
    ax1 = fg.add_axes([bx, by+sy2+sy+2*dy, sx, sy2])
    ax2 = fg.add_axes([bx, by+3*dy+2*sy+sy2, sx, sy2])
    ax3 = fg.add_axes([bx, by+4*dy+3*sy+sy2, sx, sy2])
    ax4 = fg.add_axes([bx+sx2+dx2, by, sx, sy2])
    ax5 = fg.add_axes([bx+sx2+dx, by+sy2+sy+2*dy, sx, sy2])
    ax6 = fg.add_axes([bx+sx2+dx, by+3*dy+2*sy+sy2, sx, sy2])
    ax7 = fg.add_axes([bx+sx2+dx, by+4*dy+3*sy+sy2, sx, sy2])

    ax_cb = fg.add_axes([0.2, by+sy2+dy, 0.6, 0.15])

    ax3,p1 = _create_pic_typeA(ax3,A,tA[0],'[a]')
    ax2,p2 = _create_pic_typeA(ax2,A,tA[1],'[b]')
    ax1,p3 = _create_pic_typeA(ax1,A,tA[2],'[c]')
    ax0= _create_pic_typeB(ax0,A,'[d]')

    ax7,p4 = _create_pic_typeA(ax7,B,tB[0],'[e]')
    ax6,p5 = _create_pic_typeA(ax6,B,tB[1],'[f]')
    ax5,p6 = _create_pic_typeA(ax5,B,tB[2],'[g]')
    ax4 = _create_pic_typeB(ax4,B,'[h]')

    cbar_ax,cbar =define_colorbar_S(p1,ax_cb,[-4.5,4.5],[-4.0,-2.0,0.0,2.0,4.0],r'${{H}}$ [$\mathrm{km}$]')
    ax_cb.axis('off')

 
    fg.savefig(fn,dpi=600,transparent=False)


