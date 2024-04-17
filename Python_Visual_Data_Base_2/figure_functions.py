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
from Class_Data_Base import _merge_database
from time import perf_counter 
import cmcrameri as cmc
import scipy 
import os
import latex
os.environ["PATH"] += os.pathsep + '/usr/local/texlive/2015/bin/x86_64-darwin'
print(os.getenv("PATH"))

#matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#plt.rcParams.update(params)
matplotlib.rcParams['text.usetex'] = True

"""
Function that I found in stackoverflow: 
https://stackoverflow.com/questions/33737736/matplotlib-axis-arrow-tip
"""
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
    label_ = 0.5/point_to_cm 
    axis_  = 0.4/point_to_cm
    title_ = 0.5/point_to_cm
    legend_= 0.35/point_to_cm

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
def make_figure_3(A,path_figure,figure_name,lm1,lm2):
    # figure name
    fn = os.path.join(path_figure,'%s.png'%(figure_name))
    
    # Prepare variables
    b0 = A.Det.D_x_t_det_F
    b1 = A.Det.tau_x_t_det_F

    a0 = A.C.x_sp

    i10,i20 = np.where(A.time==np.nanmin(A.Det.det_vec)),np.where(A.time==np.nanmax(A.Det.det_vec))
    string_title0 = r'$t = %s$ \textendash \ $%s$, $[Myrs]$' %("{:10.2f}".format(A.time[i10[0][0]-10]),"{:10.2f}".format(A.time[i20[0][0]]))

    # Prepare figure layout 
    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(15*cm, 15*cm))  
    bx = 0.07
    by = 0.1
    sx = 0.35
    dx = 0.08
    sy = 0.8
    dy = []
    # Prepare axis of the two end member 
    ax0 = fg.add_axes([bx, by, sx, sy])
    ax1 = fg.add_axes([bx+sx+dx, by, sx, sy])
    #cbm = fg.add_axes([bx+2*sx, 0.25, 0.01, 0.8])
    # -> coloring the plot as a function of the time
    # Solution: https://gist.github.com/charlottenosam/c63b6caa68bd117e35a4c6a7ca98007d
    # Normalize the array vals so they can be mapped to a color
    c_norm = mpl.colors.Normalize(vmin=np.nanpercentile(A.Det.det_vec,5), vmax=np.nanpercentile(A.Det.det_vec,95)) #MidpointNormalize(vmin=-0.1, vmax=A.time[i20[0][0]]-A.time[i10[0][0]],midpoint=0)#

    # Pick a colormap
    c_map  = mpl.cm.turbo

    # Scalar mappable of normalized array to colormap
    s_map  = mpl.cm.ScalarMappable(cmap=c_map, norm=c_norm)
    s_map.set_array([])
    time_det = A.time[i10[0][0]:i20[0][0]+1]
    it = i10[0][0]
    for v in time_det:
        a00=ax0.plot(a0,b0[:,it]/(A.IC.D0[0][0]/1000),linewidth = 1.0,color = s_map.to_rgba(v))
        a10=ax1.plot(a0,b1[:,it]/(A.IC.tau_co/1e6),linewidth = 1.0,color=s_map.to_rgba(v))
        it +=1 

    a01=ax0.plot(a0,b0[:,i10[0][0]-10:i10[0][0]-1]/(A.IC.D0[0][0]/1000),linewidth = 0.3,color = 'k',alpha = 0.2)
    a11=ax1.plot(a0,b1[:,i10[0][0]-10:i10[0][0]-1]/(A.IC.tau_co/1e6),linewidth = 0.3,color='k',alpha = 0.2)


    cb_ax = fg.add_axes([.88,0.25,.04,0.5])
   # cbar0 = fg.colorbar(s_map,ax=ax1,orientation='vertical',extend="both",label=r't, [Myr]',location = 'right',shrink=0.5)

   # cax = ax1.inset_axes([1.2, 0.25, 0.01, 1.0])
    #cbm.set_axis_off()
    cbar0 = fg.colorbar(s_map,cax=cb_ax,orientation='vertical',extend="both",location = 'right')
    cbar0.set_label(r't, [Myr]', labelpad=-20, y=1.1, rotation=0)

    ax0.set_ylim(lm1[0],lm1[1])
    #ax0.set_ytick(0.1,0.5,0.9)
    ax0.tick_params(axis="y",direction="in")
    ax0.tick_params(axis="x",direction="in")
    ax1.set_ylim(lm2[0],lm2[1])
    ax1.tick_params(left=True,right=True,labelleft=True) 
    ax1.tick_params(axis="y",direction="in")
    ax1.tick_params(axis="x",direction="in")
    
    plt.setp(ax0.spines.values(), linewidth=1.4)
    plt.setp(ax1.spines.values(), linewidth=1.4)

    ax0.tick_params(width=1.2)
    ax1.tick_params(width=1.2)
    
    ax0.set_xlabel(r'$x_s, [km]$',fontsize=fnt_g.label_)
    ax1.set_xlabel(r'$x_s, [km]$',fontsize=fnt_g.label_)
    ax0.set_ylabel(r'$D^{\dagger}, []$',loc='bottom',fontsize=18)
    ax1.set_ylabel(r'$\frac{\tau_{max}}{\tau_{lim}}, []$',loc='bottom',fontsize=18)

    ax0.set_xticks([200,600,1000])
    ax1.set_xticks([200,600,1000])
    ax1.set_yticks([lm2[0],(lm2[0]+lm2[1])/2,lm2[1]])
    ax0.set_yticks([lm1[0],(lm1[0]+lm1[1])/2,lm1[1]])
    ax0.yaxis.set_label_coords(-0.01,0.3)
    ax1.yaxis.set_label_coords(-0.01,0.3)

    ax0.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax0.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1.yaxis.set_tick_params(labelsize=fnt_g.axis_)

    props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)
    ax0.text(0.90, 0.98, '$[a]$', transform=ax0.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.90, 0.98, '$[b]$', transform=ax1.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')






    #ax1.set_ytick(0.1,0.5,0.9)
     

    fg.savefig(fn,dpi=600)

def make_figure_4(A,B,path_figure,figure_name):
    # figure name
    fn = os.path.join(path_figure,'%s.png'%(figure_name))
    
    # Prepare variables
    c0 = A.FS.Uplift_det[:,:,0,0]
    c1 = B.FS.Uplift_det[:,:,0,0]
    c0[c0<=np.nanmean(c0)]=np.nan
    c1[c1<=np.nanmean(c1)]=np.nan
    

    lim_0 = [np.nanmin(c0),np.nanmean(c0),np.nanmax(c0)]
    lim_1 = [np.nanmin(c1),np.nanmean(c1),np.nanmax(c1)]

    
    #c0 = (c0-np.nanmin(c0))/(np.nanmax(c0)-np.nanmin(c0))
    #c1 = (c1-np.nanmin(c1))/(np.nanmax(c1)-np.nanmin(c1))

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
    fg = figure(figsize=(15*cm, 8*cm))  
    bx = 0.07
    by = 0.1
    sx = 0.35
    dx = 0.03
    sy = 0.7
    dy = []
    ax0 = fg.add_axes([bx, by, sx, sy])
    ax1 = fg.add_axes([bx+sx+dx, by, sx, sy])
    #norm0 = MidpointNormalize(vmin=lim_0[0], vmax=lim_0[2], midpoint=lim_0[1])
    cf0=ax0.pcolormesh(a,b,c0,cmap = 'cmc.lipari')
    cbar0 = fg.colorbar(cf0, ax=ax0,orientation='horizontal',extend="both")
    cbar0.ax.tick_params(labelsize=fnt_g.legend_)
    cbar0.set_label(label=r'${\dot{H}_{mean}}$, $[mm/yr]$',size=fnt_g.label_) 
    pl0 = ax0.plot(t_x,t_y,linewidth=1.3, color = 'red')
    
       
    #norm1 = MidpointNormalize(vmin=lim_1[0], vmax=lim_1[2], midpoint=lim_1[1])
    cf1=ax1.pcolormesh(a,b,c1,cmap = 'cmc.lipari')
    cbar1 = fg.colorbar(cf1, ax=ax1,orientation='horizontal',extend="both")
    cbar1.set_label(label=r'${\dot{H}_{mean}}$, $[mm/yr]$',size=fnt_g.label_) 
    cbar1.ax.tick_params(labelsize=fnt_g.legend_)
    pl1 = ax1.plot(t_x,t_y,linewidth=1.3, color = 'red')
    ax0.set_ylim(-500,500)
    #ax0.set_ytick(0.1,0.5,0.9)
    ax0.tick_params(axis="y",direction="in")
    ax0.tick_params(axis="x",direction="in")
    ax1.set_ylim(-500,500.0)
    ax1.tick_params(left=True,right=True,labelleft=False) 
    ax1.tick_params(axis="y",direction="in")
    ax1.tick_params(axis="x",direction="in")
    
    plt.setp(ax0.spines.values(), linewidth=1.4)
    plt.setp(ax1.spines.values(), linewidth=1.4)

    ax0.tick_params(width=1.2)
    ax1.tick_params(width=1.2)
    
    ax0.set_xlabel(r'$x, [km]$',fontsize=fnt_g.label_)
    ax1.set_xlabel(r'$x, [km]$',fontsize=fnt_g.label_)
    ax0.set_ylabel(r'$y, [km]$',fontsize=fnt_g.label_)
    ax0.yaxis.set_label_coords(-0.01,0.5)
    ax0.xaxis.set_label_coords(0.5,-0.01)
    ax1.xaxis.set_label_coords(0.5,-0.01)
    ax0.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax0.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax0.set_title(string_title0,fontsize=fnt_g.title_)    #ax0.yaxis.set_label
    ax1.set_title(string_title1,fontsize=fnt_g.title_)    #ax0.yaxis.set_label

    ax0.set_xticks([-500,500])
    ax1.set_xticks([-500,500])
    ax1.set_yticks([-500,500])
    ax0.set_yticks([-500,500])


    props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)
    ax0.text(0.87, 0.95, '$[a]$', transform=ax0.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.87, 0.94, '$[b]$', transform=ax1.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')




    fg.savefig(fn,dpi=600)

    
    
def make_figure5(DB,path_save,figure_name):
    
    # figure name
    fn = os.path.join(path_save,'%s.png'%(figure_name))
    
    # Prepare variables
    vel_tearing = (DB.detachment_velocity)
    AVol        = DB.Avolume 
    T           = DB.Temp 
    SLim        = DB.StressLimit/1e6
    
    # Prepare axis of the figures 
    
    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(9*cm, 18*cm))  
    bx = 0.15
    by = 0.08
    sx = 0.8
    dx = 0.00
    sy = 0.27
    dy = 0.01
    ax0 = fg.add_axes([bx, by+2*dy+2*sy, sx, sy])
    ax1 = fg.add_axes([bx, by+1*dy+1*sy, sx, sy])
    ax2 = fg.add_axes([bx, by, sx, sy]) 
    
    colors = ['royalblue','goldenrod','tomato']
    label_fig = [r'$v_c = 10, [cm/yr$]','$v_c = 5.0, [cm/yr$]','$v_c = 2.5, [cm/yr$]']

    T_u = np.sort(np.unique(T))

    for i in range(len(T_u)):
        ax0.scatter(AVol[(SLim==200.0) & (T == T_u[i])],vel_tearing[(SLim==200.0) & (T == T_u[i])],s=50,c=colors[i],edgecolor = 'k',label=label_fig[i])
    
    ax0.legend(loc='upper center', bbox_to_anchor=(0.45, 1.20),ncol=3, columnspacing=0.05,handletextpad=0.01, shadow=True,fontsize=8)
    for i in range(len(T_u)):
        ax1.scatter(AVol[(SLim==400.0) & (T == T_u[i])],vel_tearing[(SLim==400.0) & (T == T_u[i])],s=50,c=colors[i],edgecolor = 'k',label=label_fig[i])
    
    for i in range(len(T_u)):
        ax2.scatter(AVol[(SLim==600.0) & (T == T_u[i])],vel_tearing[(SLim==600.0) & (T == T_u[i])],s=50,c=colors[i],edgecolor = 'k',label=label_fig[i])
    
    
    #ax0.set_ytick(0.1,0.5,0.9)
    ax0.tick_params(axis="y",direction="in")
    ax0.tick_params(axis="x",direction="in")
    ax1.tick_params(left=True,right=True,labelbottom=False) 
    ax0.tick_params(left=True,right=True,labelbottom=False) 

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

    
    ax2.set_xlabel(r'$V_{a,dis}$, $[\mu m^3/Pa]$',fontsize=fnt_g.label_)
    ax0.set_ylabel(r'$v_{tearing}$, $[cm/yr]$',fontsize=fnt_g.label_)
    ax1.set_ylabel(r'$v_{tearing}$, $[cm/yr]$',fontsize=fnt_g.label_)
    ax2.set_ylabel(r'$v_{tearing}$, $[cm/yr]$',fontsize=fnt_g.label_)

    
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
    ax0.text(0.87, 0.95, '$[a]$', transform=ax0.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.87, 0.94, '$[b]$', transform=ax1.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    ax2.text(0.87, 0.94, '$[c]$', transform=ax2.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    
    ax0.text(0.05, 0.15, r'$\tau_{lim} = 200 [MPa]$', transform=ax0.transAxes, fontsize=fnt_g.axis_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.05, 0.15, r'$\tau_{lim} = 400 [MPa]$', transform=ax1.transAxes, fontsize=fnt_g.axis_,
        verticalalignment='top', bbox=props,color='white')
    ax2.text(0.05, 0.15, r'$\tau_{lim} = 600 [MPa]$', transform=ax2.transAxes, fontsize=fnt_g.axis_,
        verticalalignment='top', bbox=props,color='white')
    

    
    fg.savefig(fn,dpi=600)

def make_figure6(DB,path_save,figure_name):
    
    # figure name
    fn = os.path.join(path_save,'%s.png'%(figure_name))
    
    # Prepare variables
    vel_tearing = (DB.detachment_velocity)
    AVol        = DB.Avolume 
    T           = DB.Temp 
    uplift = DB.uplift[:,0]
    SLim        = DB.StressLimit/1e6
    
    # Prepare axis of the figures 
    
    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(9*cm, 9*cm))  
    bx = 0.2
    by = 0.2
    sx = 0.6
    dx = 0.00
    sy = 0.6
    dy = 0.01
    ax0 = fg.add_axes([bx, by, sx, sy])
    
    colors = ['royalblue','goldenrod','tomato']
    label_fig = [r'$v_c = 10 [\mathrm{cm/yr}$]','$v_c = 5.0 [\mathrm{cm/yr}$]','$v_c = 2.5 [\mathrm{cm/yr}]$']

    T_u = np.sort(np.unique(T))

    for i in range(len(T_u)):
        ax0.scatter(vel_tearing[(T == T_u[i])],uplift[(T == T_u[i])],s=50,c=colors[i],edgecolor = 'k',label=label_fig[i])
    
    ax0.legend(loc='upper center', bbox_to_anchor=(0.45, 1.15),ncol=3, columnspacing=0.02,handletextpad=0.005, shadow=True,fontsize=8)
    
    ax0.tick_params(axis="y",direction="in")
    ax0.tick_params(axis="x",direction="in")
    ax0.tick_params(left=True,right=True,labelbottom=True) 

    
    plt.setp(ax0.spines.values(), linewidth=1.4)
   

    ax0.tick_params(width=1.2)
   
    
    ax0.set_xlabel(r'$v_{\mathrm{tearing}}$, $[\mathrm{cm/yr}]$',fontsize=fnt_g.label_)
    ax0.set_ylabel(r'$\dot{H}_{\mathrm{mean}}$, $[\mathrm{mm/yr}]$',fontsize=fnt_g.label_)

    
    ax0.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax0.yaxis.set_tick_params(labelsize=fnt_g.axis_)

    ax0.set_xscale('log')
    ax0.set_yscale('log')

    
    fg.savefig(fn,dpi=600)
    
def initial_geotherm(path_save):
    z = np.linspace(start=0.0,stop=150.0,num=100)
    d = np.linspace(start=0.0,stop=100.0,num=100)
    L0 = 400.0 
    T_continent = compute_geoterhm(35.0,20,1350,600,z)
    T_Oplate    = compute_OC_geoterhm(20,1350,z,1e-6,30e6*60*60*365.25*24)
    T_Oplate    = (T_Oplate[1:]+T_Oplate[0:-1])/2
    T_slab_v0,l,lc   = compute_McKenzie_Field(20,1350,2.5,d,L0,-100)
    T_slab_v1,l,lc   = compute_McKenzie_Field(20,1350,5.0,d,L0,-100)
    T_slab_v2,l,lc   = compute_McKenzie_Field(20,1350,10.0,d,L0,-100)
    
    fn = os.path.join(path_save,'Initial_Temperature.svg')
    
    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(10*cm, 12*cm))  
    bx = 0.15
    by = 0.15
    sx1 = 0.35
    sx2 = 0.20
    dx = 0.06
    dx1 = 0.08
    sy1 = 0.25
    sy2 = 0.45
    dy = 0.05
    ax0 = fg.add_axes([bx, by+dy+sy2, sx1, sy1])
    ax1 = fg.add_axes([bx+dx+sx1, by+dy+sy2, sx1, sy1])
    ax2 = fg.add_axes([bx, by, sx2, sy2])
    ax3 = fg.add_axes([bx+dx1+sx2, by, sx2, sy2])
    ax4 = fg.add_axes([bx+2*dx1+2*sx2, by, sx2, sy2])
    plt.setp(ax0.spines.values(), linewidth=1.4)
    plt.setp(ax1.spines.values(), linewidth=1.4)
    plt.setp(ax2.spines.values(), linewidth=1.4)
    plt.setp(ax3.spines.values(), linewidth=1.4)
    plt.setp(ax4.spines.values(), linewidth=1.4)
    
    ax0.tick_params(axis="y",direction="in")
    ax0.tick_params(axis="x",direction="in")
    
    ax1.tick_params(axis="y",direction="in")
    ax1.tick_params(axis="x",direction="in")
    
    ax2.tick_params(axis="y",direction="in")
    ax2.tick_params(axis="x",direction="in")
    
    ax3.tick_params(axis="y",direction="in")
    ax3.tick_params(axis="x",direction="in")
    
    ax4.tick_params(axis="y",direction="in")
    ax4.tick_params(axis="x",direction="in")




    ax0.tick_params(width=1.2)
    ax1.tick_params(width=1.2)
    ax2.tick_params(width=1.2)
    ax3.tick_params(width=1.2)
    ax4.tick_params(width=1.2)

    
    ax0.spines['bottom'].set_color('white')
    ax0.spines['top'].set_color('black') 
    ax0.spines['right'].set_color('white')
    ax0.spines['left'].set_color('black')
    ax0.xaxis.tick_top()
    ax0.tick_params(labelbottom=False,labeltop=True)
    
    
    ax0.plot(T_continent[z<100],-z[z<100],linewidth = 2.0, color = 'red')
    ax0.set_ylim(-120,0)
    ax0.set_xlim(20,1350)
    
    

    
    ax1.spines['bottom'].set_color('white')
    ax1.spines['top'].set_color('black') 
    ax1.spines['right'].set_color('white')
    ax1.spines['left'].set_color('black')
    ax1.xaxis.tick_top()
    z_alt = (z[1:]+z[:-1])/2
    ax1.plot(T_Oplate[z_alt<100],-z_alt[z_alt<100],linewidth = 2.0, color = 'blue')
    ax1.set_ylim(-120,0)
    ax1.set_xlim(20,1350)
    ax1.tick_params(labelbottom=False,labeltop=True)

    
    
    t_field = np.linspace(start=0,stop=1200,num=12)
    ax2.contourf(-d,-l,np.transpose(T_slab_v0),cmap = 'cmc.lipari')
    ax2.axhline(y=-lc,color='lightgreen',linewidth=1.2)
    ax3.contourf(-d,-l,np.transpose(T_slab_v1),cmap = 'cmc.lipari')
    ax3.axhline(y=-lc,color='lightgreen',linewidth=1.2)
    ax4.contourf(-d,-l,np.transpose(T_slab_v2),cmap = 'cmc.lipari')
    ax4.axhline(y=-lc,color='lightgreen',linewidth=1.2)


    ax1.tick_params(left=True,right=False,labelleft=False) 
    ax3.tick_params(left=True,right=True,labelleft=False) 
    ax4.tick_params(left=True,right=True,labelleft=False) 

    ax0.set_xlabel(r'$T$, $[^{\circ}C]$',fontsize=fnt_g.label_)
    ax0.set_ylabel(r'$z$, $[km]$',fontsize=fnt_g.label_)
    ax1.set_xlabel(r'$T$, $[^{\circ}C]$',fontsize=fnt_g.label_)

    ax2.set_xlabel(r'$d$, $[km]$',fontsize=fnt_g.label_)
    ax2.set_ylabel(r'$\ell$, $[km]$',fontsize=fnt_g.label_)
    ax3.set_xlabel(r'$d$, $[km]$',fontsize=fnt_g.label_)
    ax4.set_xlabel(r'$d$, $[km]$',fontsize=fnt_g.label_)

    ax0.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax0.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax2.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax2.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax3.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax3.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax4.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax4.yaxis.set_tick_params(labelsize=fnt_g.axis_)

    ax0.set_xticks([0.0,600.0,1200.0])
    ax1.set_xticks([0.0,600.0,1200.0])
    ax0.set_xticklabels([r'$0$','','$1200$'])
    ax1.set_xticklabels([r'$0$','','$1200$'])




    ax0.set_yticks([-100,-35,0])
    ax1.set_yticks([-100,-35,0])

    ax2.set_yticks([-400,-200,0])
    ax3.set_yticks([-400,-200,0])
    ax4.set_yticks([-400,-200,0])
    
    l_abel = ['$400$','','$0$']
    ax2.set_yticklabels(l_abel)
    ax2.set_xticks([-100.0,-50,0.0])
    ax3.set_xticks([-100.0,-50,0.0])
    ax4.set_xticks([-100.0,-50,0.0])
    ax2.set_xticklabels([r'$100$','','$0$'])
    ax3.set_xticklabels([r'$100$','','$0$'])
    ax4.set_xticklabels([r'$100$','','$0$'])
    ax0.yaxis.set_label_coords(-0.03,0.4)
    ax2.yaxis.set_label_coords(-0.03,0.5)
    
    ax0.xaxis.set_label_coords(0.5,1.2)
    ax1.xaxis.set_label_coords(0.5,1.2)
    ax2.xaxis.set_label_coords(0.5,-0.03)
    ax3.xaxis.set_label_coords(0.5,-0.03)

    ax4.xaxis.set_label_coords(0.5,-0.03)
    
    
    props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)
    ax0.text(0.87, 0.90, '$[b]$', transform=ax0.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.87, 0.90, '$[c]$', transform=ax1.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    ax2.text(0.08, 0.96, '$[d]$', transform=ax2.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    ax3.text(0.08, 0.96, '$[e]$', transform=ax3.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    ax4.text(0.08, 0.96, '$[f]$', transform=ax4.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    
    arrowed_spines(fg, ax0)
    arrowed_spines(fg, ax1)



    fg.savefig(fn,dpi=600,transparent=True)

    



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
    """
    
     a = (-1).^(i)./(i.*pi);
    b = (Re-(Re.^2+i^2*pi^2).^(0.5)).*obj.l_slab(~isnan(obj.d_slab)).*sc;
    c = sin(i.*pi.*(1-abs(obj.d_slab(~isnan(obj.d_slab)).*sc)));
    e = exp(b);
    """
    
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

    
    
    
    