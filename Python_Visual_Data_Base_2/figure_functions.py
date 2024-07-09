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

from time import perf_counter 
plugin='DICOM'
import cmcrameri as cmc
import scipy 
import os
import latex
import imageio

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
    label_ = 0.45/point_to_cm 
    axis_  = 0.4/point_to_cm
    title_ = 0.45/point_to_cm
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
    b0 = A.Det.D_x_t_det
    b1 = A.Det.tau_x_t_det

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
    c_map  = cmc.batlow

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
    a11=ax1.plot(a0,b1[:,10:i10[0][0]-1]/(A.IC.tau_co/1e6),linewidth = 0.3,color='k',alpha = 0.2)


    cb_ax = fg.add_axes([.88,0.25,.04,0.5])
   # cbar0 = fg.colorbar(s_map,ax=ax1,orientation='vertical',extend="both",label=r't, [Myr]',location = 'right',shrink=0.5)

   # cax = ax1.inset_axes([1.2, 0.25, 0.01, 1.0])
    #cbm.set_axis_off()
    cbar0 = fg.colorbar(s_map,cax=cb_ax,orientation='vertical',extend="both",location = 'right')
    cbar0.set_label(r'$t$ /Myr', labelpad=-20, y=1.1, rotation=0)

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
    
    ax0.set_xlabel(r'$x_{strike}$ /[km]',fontsize=fnt_g.label_)
    ax1.set_xlabel(r'$x_{strike}$ /[km]',fontsize=fnt_g.label_)
    ax0.set_ylabel(r'$D^{\dagger}$ /n.d.',loc='bottom',fontsize=18)
    ax1.set_ylabel(r'$\tau^{\dagger}_{max}$ /n.d.',loc='bottom',fontsize=18)

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

def make_figure_4(A,B,path_figure,figure_name,det):
    # figure name
    fn = os.path.join(path_figure,'%s.png'%(figure_name))
    
    # Prepare variables
    if det == 1: 
        c0 = A.FS.total_uplift_NT
        c1 = B.FS.total_uplift_NT
    elif det == 2:
        c0 = A.FS.total_uplift_Te
        c1 = B.FS.total_uplift_Te
    else: 
        c0 = A.FS.total_uplift_LT
        c1 = B.FS.total_uplift_LT
        
 
    #c0[c0<=np.nanmean(c0)]=np.nan
    #c1[c1<=np.nanmean(c1)]=np.nan
    

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
    sx = 0.45
    dx = 0.03
    sy = 0.7
    dy = []
    ax0 = fg.add_axes([bx, by, sx, sy])
    ax1 = fg.add_axes([bx+sx+dx, by, sx, sy])
    norm0 = MidpointNormalize(vmin=lim_0[0], vmax=lim_0[2], midpoint=0.0)
    levels_0 = np.linspace(np.round(lim_0[0]),np.round(lim_0[2]),20)
    cf0=ax0.contourf(a,b,c0,cmap = 'cmc.cork',norm=norm0,levels = levels_0)
    cbar0 = fg.colorbar(cf0, ax=ax0,orientation='horizontal',extend="both",norm=norm0,ticks=[np.round(lim_0[0]), 0.0, np.round(lim_0[2])],shrink = 0.8)
    cbar0.ax.tick_params(labelsize=fnt_g.legend_)
    cbar0.set_label(label=r'${{\Delta H}}$ /$\mathrm{m}$',size=fnt_g.label_) 
    pl0 = ax0.plot(t_x,t_y,linewidth=1.3, color = 'red')
    
       
    norm1 = MidpointNormalize(vmin=lim_1[0], vmax=lim_1[2], midpoint=0.0)
    levels_1 = np.linspace(np.round(lim_1[0]),np.round(lim_1[2]),20)
    cf1=ax1.contourf(a,b,c1,cmap = 'cmc.cork',norm=norm1,levels=levels_1)
    cbar1 = fg.colorbar(cf1, ax=ax1,orientation='horizontal',extend="both",norm=norm1,ticks=[np.round(lim_1[0]), 0.0, np.round(lim_1[2])],shrink = 0.8)
    cbar1.set_label(label=r'${{\Delta H}}$ /$\mathrm{m}$',size=fnt_g.label_) 
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
    
    ax0.set_xlabel(r'$x$, /km',fontsize=fnt_g.label_)
    ax1.set_xlabel(r'$x$, /km',fontsize=fnt_g.label_)
    ax0.set_ylabel(r'$y$, /km',fontsize=fnt_g.label_)
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
    label_fig = [r'$v_c = 10$ /$\frac{\mathrm{cm}}{\mathrm{yr}}$',r'$v_c = 5.0$ /$\frac{\mathrm{cm}}{\mathrm{yr}}$',r'$v_c = 2.5$ /$\frac{\mathrm{cm}}{\mathrm{yr}}$']

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
    
    ax0.set_ylabel(r'$v_{tearing}$ /$\frac{\mathrm{cm}}{\mathrm{yr}}$',fontsize=fnt_g.label_)

    ax0.set_xlabel(r'$V_{a,dis}$, /$\mu \frac{\mathrm{m}^3}{\mathrm{Pa}}$',fontsize=fnt_g.label_)
    ax1.set_xlabel(r'$V_{a,dis}$, /$\mu \frac{\mathrm{m}^3}{\mathrm{Pa}}$',fontsize=fnt_g.label_)
    ax2.set_xlabel(r'$V_{a,dis}$, /$\mu \frac{\mathrm{m}^3}{\mathrm{Pa}}$',fontsize=fnt_g.label_)

    
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

def make_figure6(DB,path_save,figure_name):
    
    # figure name
    
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
    Uplift_discrete = DB.Uplift_max_discrete
    for ip in range(3):
        fg = figure(figsize=(18*cm, 9*cm)) 
        ax0 = fg.add_axes([bx, by, sx, sy])
        ax1 = fg.add_axes([bx+sx+dx,by,sx,sy])

        fn = os.path.join(path_save,'%s%s.png'%(figure_name,type[ip]))

        uplift = DB.uplift[:,ip]
        colors = ['royalblue','goldenrod','tomato','orange','grey','pink']
        label_fig = [r'$V_{a,dis} = 8$ /$\frac{\mathrm{m^3}}{\mathrm{Pa}}$',r'$V_{a,dis} = 10$ /$\frac{\mathrm{m^3}}{\mathrm{Pa}}$',r'$V_{a,dis} = 11$ /$\frac{\mathrm{m^3}}{\mathrm{Pa}}$',r'$V_{a,dis} = 12$ /$\frac{\mathrm{m^3}}{\mathrm{Pa}}$',r'$V_{a,dis} = 13 $ /$\frac{\mathrm{m^3}}{\mathrm{Pa}}$',r'$V_{a,dis} = 15$ /$\frac{\mathrm{m^3}}{\mathrm{Pa}}$']

        AVol_u = np.sort(np.unique(AVol))
        AVol_u = AVol_u[AVol_u>0.0]
        for i in range(len(AVol_u)):
            ax0.scatter(vel_tearing[(AVol == AVol_u[i])],uplift[(AVol == AVol_u[i])],c=colors[i],s=50,edgecolor = 'k',label=label_fig[i])
            ax1.scatter(vel_tearing[(AVol == AVol_u[i])],Uplift_discrete[(AVol == AVol_u[i])],c=colors[i],s=50,edgecolor = 'k',label=label_fig[i])
        ax0.axvline(2,linewidth=0.8,color='k',alpha=0.5)
        ax0.axvline(94,linewidth=0.8,color='k',alpha=0.5)
        ax1.axvline(2,linewidth=0.8,color='k',alpha=0.5)
        ax1.axvline(94,linewidth=0.8,color='k',alpha=0.5)


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
    

        ax0.set_xlabel(r'$v_{\mathrm{tearing}}$ /$\frac{\mathrm{cm}}{\mathrm{yr}}$',fontsize=fnt_g.label_)
        ax0.set_ylabel(r'$\dot{H}_{\mathrm{mean}}$ /$\frac{\mathrm{mm}}{\mathrm{yr}}$',fontsize=fnt_g.label_)
        ax1.set_xlabel(r'$v_{\mathrm{tearing}}$ /$\frac{\mathrm{cm}}{\mathrm{yr}}$',fontsize=fnt_g.label_)


        ax0.xaxis.set_tick_params(labelsize=fnt_g.axis_)
        ax0.yaxis.set_tick_params(labelsize=fnt_g.axis_)
        
        
        ax1.xaxis.set_tick_params(labelsize=fnt_g.axis_)
        ax1.yaxis.set_tick_params(labelsize=fnt_g.axis_)

        ax0.set_xscale('log')
        ax1.set_xscale('log')
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

    ax0.set_xlabel(r'$T$ /$^{\circ}\mathrm{C}$',fontsize=fnt_g.label_)
    ax0.set_ylabel(r'$z$ /km',fontsize=fnt_g.label_)
    ax1.set_xlabel(r'$T$ /$^{\circ}\mathrm{C}$',fontsize=fnt_g.label_)

    ax2.set_xlabel(r'$d$ /km',fontsize=fnt_g.label_)
    ax2.set_ylabel(r'$\ell$ /km',fontsize=fnt_g.label_)
    ax3.set_xlabel(r'$d$ /km',fontsize=fnt_g.label_)
    ax4.set_xlabel(r'$d$ /km',fontsize=fnt_g.label_)

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

def make_plot_dD(path_save,A,B):
    
    fn = os.path.join(path_save,'Test.png')
    D0 = 100 
    
    a0 = A.time
    a1 = A.Det.deltaD
    a2 = A.Det.minD
    a3 = A.Det.maxD 
    a4 = A.Det.meanD

    b0 = B.time
    b1 = B.Det.deltaD
    b2 = B.Det.minD
    b3 = B.Det.maxD 
    b4 = B.Det.meanD  
  

    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(12*cm, 6*cm))  
    bx = 0.15
    by = 0.2
    sx = 0.35
    dx = 0.10
    sy = 0.7
    
    ax0 = fg.add_axes([bx,by,sx,sy])
    ax0b = ax0.twinx()
    ax1 = fg.add_axes([bx+sx+dx,by,sx,sy])
    ax1b = ax1.twinx()

    ax0.plot(a0,a1,linewidth = 1.2, color = 'k')
    ax0b.fill_between(a0,a2,a3,alpha=0.4,color='firebrick')
    ax0b.plot(a0,a4,linewidth = 0.5,alpha=1.0,color='firebrick')
    

    ax1.plot(b0,b1,linewidth = 1.2, color = 'k')
    ax1b.fill_between(b0,b2,b3,alpha=0.4,color='firebrick')
    ax1b.plot(b0,b4,linewidth = 0.5,alpha=1.0,color='firebrick')


    plt.setp(ax0.spines.values(), linewidth=1.4)
    plt.setp(ax1.spines.values(), linewidth=1.4)
    
    plt.setp(ax0b.spines.values(), linewidth=1.4)
    plt.setp(ax1b.spines.values(), linewidth=1.4)
    

    ax0.tick_params(axis="y",direction="in")
    ax0.tick_params(axis="x",direction="in")
    
    ax1.tick_params(axis="y",direction="in")
    ax1.tick_params(axis="x",direction="in")


    ax0.set_ylabel(r'$\Delta D^{\dagger}$ /n.d.',fontsize=fnt_g.label_)
    ax1b.set_ylabel(r'$D^{\dagger}$ /n.d.',fontsize=fnt_g.label_)
    ax0.set_xlabel(r'$t$ /Myr',fontsize=fnt_g.label_)
    ax1.set_xlabel(r'$t$ /Myr',fontsize=fnt_g.label_)



    
    ax0.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax0.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax0b.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax0b.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1b.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1b.yaxis.set_tick_params(labelsize=fnt_g.axis_)


    ax0.set_yticks([0.0,round(np.nanmax(a1),1)])
    ax1.set_yticks([0.0,round(np.nanmax(b1),1)])
    ax0b.set_yticks([0.0,1.0])
    ax1b.set_yticks([0.0,1.0])


    ax0.yaxis.set_label_coords(-0.01,0.5)
    ax1b.yaxis.set_label_coords(1.01,0.5)


    ax0.spines['bottom'].set_color('k')
    ax0.spines['top'].set_color('w') 
    ax0.spines['right'].set_color('firebrick')
    ax0.spines['left'].set_color('k')

    ax1.spines['bottom'].set_color('k')
    ax1.spines['top'].set_color('w') 
    ax1.spines['right'].set_color('firebrick')
    ax1.spines['left'].set_color('k')

    ax0b.spines['bottom'].set_color('k')
    ax0b.spines['top'].set_color('w') 
    ax0b.spines['right'].set_color('firebrick')
    ax0b.spines['left'].set_color('k')

    ax1b.spines['bottom'].set_color('k')
    ax1b.spines['top'].set_color('w') 
    ax1b.spines['right'].set_color('firebrick')
    ax1b.spines['left'].set_color('k')


    ax0.spines[['right', 'top']].set_visible(False)
    ax1.spines[['right', 'top']].set_visible(False)
    ax0b.spines[[ 'top']].set_visible(False)
    ax1b.spines[[ 'top']].set_visible(False)


    ax1b.yaxis.label.set_color('firebrick')
    ax1b.tick_params(axis='y', colors='firebrick')
    ax0b.tick_params(axis='y', colors='firebrick')



    props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.5)
    ax0.text(0.82, 1.03, '$[a]$', transform=ax0.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.82, 1.03, '$[b]$', transform=ax1.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')

    
    fg.savefig(fn,dpi=600,transparent=False)
 
def make_figure7(DB,path_save,figure_name):
    
    # figure name
    fn = os.path.join(path_save,'%s.png'%(figure_name))
    
    # Prepare variables
    vel_tearing = (DB.detachment_velocity)
    AVol        = DB.Avolume 
    T           = DB.Temp 
    SLim        = DB.tau_max
    type = ['det','neck','LT']

    for i in range(3):

        fn = os.path.join(path_save,'%s%s.png'%(figure_name,type[i]))

        uplift = DB.uplift[:,0]
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
        label_fig = [r'$v_c = 10$ /$\frac{\mathrm{cm}}{\mathrm{yr}}$',r'$v_c = 5.0$ /$\frac{\mathrm{cm}}{\mathrm{yr}}$',r'$v_c = 2.5$ /$\frac{\mathrm{cm}}{\mathrm{yr}}$']

        T_u = np.sort(np.unique(T))
        T_u = T_u[T_u>0.0]


        for i in range(len(T_u)):
            ax0.scatter(AVol[(T == T_u[i])],SLim[(T == T_u[i])],s=50,c=colors[i],edgecolor = 'k',label=label_fig[i])

        ax0.legend(loc='upper center', bbox_to_anchor=(0.45, 1.15),ncol=3, columnspacing=0.02,handletextpad=0.005, shadow=True,fontsize=8)

        ax0.tick_params(axis="y",direction="in")
        ax0.tick_params(axis="x",direction="in")
        ax0.tick_params(left=True,right=True,labelbottom=True) 


        plt.setp(ax0.spines.values(), linewidth=1.4)
    

        ax0.tick_params(width=1.2)
    

        ax0.set_xlabel(r'$V_{a,\mathrm{dis}}$, $[\mathrm{\mu m^3/Pa}]$',fontsize=fnt_g.label_)
        ax0.set_ylabel(r'$\tau_{\mathrm{max}}/\tau_{\mathrm{lim}}$, $[]$',fontsize=fnt_g.label_)


        ax0.xaxis.set_tick_params(labelsize=fnt_g.axis_)
        ax0.yaxis.set_tick_params(labelsize=fnt_g.axis_)

        ax0.set_xscale('linear')
        ax0.set_yscale('linear')


        fg.savefig(fn,dpi=600)
       
    
def make_figure8(DB,path_save,figure_name):
    
    # figure name
    fn = os.path.join(path_save,'%s.png'%(figure_name))
    
    # Prepare variables
    vel_tearing = (DB.detachment_velocity)
    AVol        = DB.Avolume 
    T           = DB.Temp 
    SLim        = DB.StressLimit/1e6
    tau_M       = DB.tau_max
    
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
    label_fig = [r'$v_c = 10$ /$\frac{\mathrm{cm}}{\mathrm{yr}}$',r'$v_c = 5.0$ /$\frac{\mathrm{cm}}{\mathrm{yr}}$',r'$v_c = 2.5$ /$\frac{\mathrm{cm}}{\mathrm{yr}}$']

    T_u = np.sort(np.unique(T))
    T_u = T_u[T_u>0.0]


    for i in range(len(T_u)):
        ax0.scatter(AVol[(SLim==200.0) & (T == T_u[i])],tau_M[(SLim==200.0) & (T == T_u[i])],s=50,c=colors[i],edgecolor = 'k',label=label_fig[i])
    
    ax0.legend(loc='upper center', bbox_to_anchor=(0.45, 1.20),ncol=3, columnspacing=0.05,handletextpad=0.01, shadow=True,fontsize=8)
    for i in range(len(T_u)):
        ax1.scatter(AVol[(SLim==400.0) & (T == T_u[i])],tau_M[(SLim==400.0) & (T == T_u[i])],s=50,c=colors[i],edgecolor = 'k',label=label_fig[i])
    
    for i in range(len(T_u)):
        ax2.scatter(AVol[(SLim==600.0) & (T == T_u[i])],tau_M[(SLim==600.0) & (T == T_u[i])],s=50,c=colors[i],edgecolor = 'k',label=label_fig[i])
    
    
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
    ax0.set_ylabel(r'$\tau_{\mathrm{max}}/\tau_{\mathrm{lim}}$',fontsize=fnt_g.label_)
    ax1.set_ylabel(r'$\tau_{\mathrm{max}}/\tau_{\mathrm{lim}}$',fontsize=fnt_g.label_)
    ax2.set_ylabel(r'$\tau_{\mathrm{max}}/\tau_{\mathrm{lim}}$',fontsize=fnt_g.label_)

    
    ax0.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax0.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax2.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax2.yaxis.set_tick_params(labelsize=fnt_g.axis_)

    ax0.set_yscale('linear')
    ax1.set_yscale('linear')
    ax2.set_yscale('linear')


    
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

@timer
def _make_gif(test,ptsave_b):

    test_gf = os.path.join(ptsave_b,test.Test_Name)
    if not os.path.isdir(test_gf): 
        os.mkdir(test_gf)    
    
    cm = 1/2.54  # centimeters in inches
    file_list = []

    for ipic in range(len(test.time)):
        if test.time[ipic]>1.0:
            dH_trench = np.zeros(len(test.C.x_trench_p),dtype=float)

            # interpolate dh along the trench (filtered)
            dH_trench = _interpolate_2D(test.FS.dH[:,:,ipic],test.C.xg,test.C.yg,test.C.x_trench_p,test.C.y_trench_p)


            fna='Fig'+str(ipic)+'.png'
            fg = figure(figsize=(14*cm, 10*cm))  
            tick=r'$Time = %.3f$ /Myr' %(test.time[ipic])
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
            if (np.isnan(np.nanmin(dH_trench)) == False) and  (np.isnan(np.nanmax(dH_trench)) == False):
    #secax_y0b1= ax0b.secondary_yaxis(1.2, functions=(celsius_to_anomaly, anomaly_to_celsius))
    #secax_y0b1set_ylabel(r'$T - \ overline{T}\ [^oC]$')
                ax1.plot(test.C.x_sp,dH_trench,color='k',linewidth=1.2)
                ax1.axhline(y=0.0, color = 'k', linestyle=':', linewidth = 0.4)
                ax1.set_yticks([round(np.nanmin(dH_trench),2),round(np.nanmax(dH_trench),2)])
            else:
                ax1.set_yticks([-1,1])

            ax1.set_ylabel(r'$\dot{H}$ /$\frac{\mathrm{mm}}{\mathrm{yr}}$',size=fnt_g.label_)
            ax1b.plot(test.C.x_sp,test.Det.tau_x_t_det[:,ipic]/(test.IC.tau_co/1e6),color='forestgreen',linewidth=1.2) 
            ax1b.plot(test.C.x_sp,test.Det.D_x_t_det[:,ipic]/100,color='firebrick',linewidth=1.2) 
            ax1b.set_ylabel(r'$\tau^{\dagger}_{max}$ /n.d.',size=fnt_g.label_)
            ax1b.set_yticks([0.0,1.0]) 
            secax_y0b = ax1b.secondary_yaxis(1.12)
            secax_y0b.set_ylabel(r'$D^{\dagger}$ /n.d.',size=fnt_g.label_)
            secax_y0b.set_ylim([0.0,1.0])
            secax_y0b.spines['right'].set_color('white')
    
            ax1b.yaxis.label.set_color('forestgreen')
            ax1b.tick_params(axis='y', colors='k')
            secax_y0b.yaxis.label.set_color('firebrick')
            secax_y0b.tick_params(axis='y',color='white')
    
    #ax0.xaxis.set_label_coords(0.5,-0.02)
            ax1.yaxis.set_label_coords(-0.1,0.5)

            ax1b.yaxis.set_label_coords(1.03,0.5)
            secax_y0b.yaxis.set_label_coords(1.4,0.5)

           
            secax_y0b.yaxis.set_tick_params(labelsize=fnt_g.axis_)
            secax_y0b.set_yticklabels(['',''])


            
            plt.setp(ax0.spines.values(), linewidth=1.4)
            plt.setp(ax1.spines.values(), linewidth=1.4)
            plt.setp(ax1b.spines.values(), linewidth=1.4)

            
            buf = np.log10(test.Det.Psi[:,:,ipic]) # Hard coded value i know. 
            lim_psi = [round(np.nanpercentile(np.log10(test.Det.Psi[:,:,ipic]),5),0),round(np.nanpercentile(np.log10(test.Det.Psi[:,:,ipic]),95),0)]
            levels = np.linspace(np.round(0.1), np.round(0.75), num=10, endpoint=True, retstep=False, dtype=float)
            cf=ax0.pcolormesh(test.C.x_sp,test.C.zp,np.transpose(buf),cmap='cmc.lipari',vmin = lim_psi[0], vmax=lim_psi[1])
            cbar = fg.colorbar(cf,ax=ax0,orientation='horizontal',label=r'$\dot{\Psi}$ /$\mathrm{\frac{W}{m^3}}$',extend="both",shrink=0.5)
            
            depth = test.Det.depth_vec
            depth[depth>=-80]=np.nan
            
            condition = (test.C.x_sp > 100) & (test.C.x_sp <1100)
            cf2=ax0.plot(test.C.x_sp[condition==1],depth[condition==1],color='forestgreen',linewidth=0.1)
            dummy_var = test.C.x_sp[condition==1]/test.C.x_sp[condition==1]
            cf2=ax0.fill_between(test.C.x_sp[condition==1],dummy_var*np.nanmin(depth[condition==1]),dummy_var*np.nanmax(depth[condition==1]),color='forestgreen',alpha=0.2,linewidth=1.2)
            p1 = ax0.contour(test.C.x_sp,test.C.zp,np.transpose(test.Det.T[:,:,ipic]),levels = [800,900,1000,1100,1200],colors = 'k',linewidths=0.5)
            ax0.clabel(p1, p1.levels, inline=True, fmt=fmt2, fontsize=6)
            
            ax1.set_title(tick)
            ax0.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
            ax0.set_ylabel(r'$z$ /km')
            ax0.set_xlabel(r'$x_{strike}$ /km')
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
    





def make_figure_3N(A,path_figure,figure_name,time):
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
    
    tick0 = r'$t = %.2f$ /Myr' %(t0)
    tick1 = r'$t = %.2f$ /Myr' %(t1)
    tick2 = r'$t = %.2f$ /Myr' %(t2)
    
    ind0 = np.where(A.time>=t0)
    ind1 = np.where(A.time>=t1)
    ind2 = np.where(A.time>=t2)
    
    ind0 = ind0[0][0]
    ind1 = ind1[0][0]
    ind2 = ind2[0][0]
    
    dH0 = _interpolate_2D(A.FS.dH_fil[:,:,ind0],A.C.xg,A.C.yg,A.C.x_trench_p,A.C.y_trench_p)
    dH1 = _interpolate_2D(A.FS.dH_fil[:,:,ind1],A.C.xg,A.C.yg,A.C.x_trench_p,A.C.y_trench_p)
    dH2 = _interpolate_2D(A.FS.dH_fil[:,:,ind2],A.C.xg,A.C.yg,A.C.x_trench_p,A.C.y_trench_p)

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
    
    #secax_y0b = ax0b.secondary_yaxis(1.2, functions=(celsius_to_anomaly, anomaly_to_celsius))
    #secax_y0b.set_ylabel(r'$T - \overline{T}\ [^oC]$')
    ax0.plot(A.C.x_sp,dH0,color='k',linewidth=1.2)
    ax0.axhline(y=0.0, color = 'k', linestyle=':', linewidth = 0.4)
    ax0.set_yticks([round(np.min(dH0),2),round(np.max(dH0),2)])
    ax0.set_ylabel(r'$\dot{H}$ /$\frac{\mathrm{mm}}{\mathrm{yr}}$',size=fnt_g.label_)
    ax0b.plot(A.C.x_sp,A.Det.tau_x_t_det[:,ind0]/(A.IC.tau_co/1e6),color='forestgreen',linewidth=1.2) 
    ax0b.plot(A.C.x_sp,A.Det.D_x_t_det[:,ind0]/100,color='firebrick',linewidth=1.2) 
    ax0b.set_ylabel(r'$\tau^{\dagger}_{max}$ /n.d.',size=fnt_g.label_)
    ax0b.set_yticks([0.0,1.0]) 
    secax_y0b = ax0b.secondary_yaxis(1.12)
    secax_y0b.set_ylabel(r'$D^{\dagger}$ /n.d.',size=fnt_g.label_)
    secax_y0b.set_ylim([0.0,1.0])
    secax_y0b.spines['right'].set_color('white')
    
    ax0b.yaxis.label.set_color('forestgreen')
    ax0b.tick_params(axis='y', colors='k')
    secax_y0b.yaxis.label.set_color('firebrick')
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
    ax1.set_yticks([round(np.min(dH1),2),round(np.max(dH1),2)])
    ax1.set_ylabel(r'$\dot{H}$ /$\frac{\mathrm{mm}}{\mathrm{yr}}$',size=fnt_g.label_)
    ax1b.plot(A.C.x_sp,A.Det.tau_x_t_det[:,ind1]/(A.IC.tau_co/1e6),color='forestgreen',linewidth=1.2) 
    ax1b.plot(A.C.x_sp,A.Det.D_x_t_det[:,ind1]/100,color='firebrick',linewidth=1.2) 
    ax1b.set_ylabel(r'$\tau^{\dagger}_{max}$ /n.d.',size=fnt_g.label_)
    ax1b.set_yticks([0.0,1.0]) 
    secax_y1b = ax1b.secondary_yaxis(1.12)
    secax_y1b.set_ylabel(r'$D^{\dagger}$ /n.d.',size=fnt_g.label_)
    secax_y1b.set_ylim([0.0,1.0])
    secax_y1b.spines['right'].set_color('white')
    
    ax1b.yaxis.label.set_color('forestgreen')
    ax1b.tick_params(axis='y', colors='k')
    secax_y1b.yaxis.label.set_color('firebrick')
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
    ax2.set_yticks([round(np.min(dH2),2),round(np.max(dH2),2)])
    ax2.set_ylabel(r'$\dot{H}$ /$\frac{\mathrm{mm}}{\mathrm{yr}}$',size=fnt_g.label_)
    ax2b.plot(A.C.x_sp,A.Det.tau_x_t_det[:,ind2]/(A.IC.tau_co/1e6),color='forestgreen',linewidth=1.2) 
    ax2b.plot(A.C.x_sp,A.Det.D_x_t_det[:,ind2]/100,color='firebrick',linewidth=1.2) 
    ax2b.set_ylabel(r'$\tau^{\dagger}_{max}$ /n.d.',size=fnt_g.label_)
    ax2b.set_yticks([0.0,1.0]) 
    secax_y2b = ax2b.secondary_yaxis(1.12)
    secax_y2b.set_ylabel(r'$D^{\dagger}$ /n.d.',size=fnt_g.label_)
    secax_y2b.set_ylim([0.0,1.0])
    secax_y2b.spines['right'].set_color('white')
    
    ax2b.yaxis.label.set_color('forestgreen')
    ax2b.tick_params(axis='y', colors='k')
    secax_y2b.yaxis.label.set_color('firebrick')
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

        


    ax2.set_xticklabels([r'$200$','','$1200$'])
    ax2.set_xlabel(r'$x_{trench}$ /km',size=fnt_g.label_)

    props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)
    ax0.text(0.025, 1.15, tick0, transform=ax0.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.025, 1.15, tick1, transform=ax1.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
    ax2.text(0.025, 1.15, tick2, transform=ax2.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
    
    ax0.spines[['top']].set_visible(False)
    ax0b.spines[['left','top']].set_visible(False)
    ax1.spines[['top']].set_visible(False)
    ax1b.spines[['left','top']].set_visible(False)
    ax2.spines[['top']].set_visible(False)
    ax2b.spines[['left','top']].set_visible(False)

    #ax2.yaxis.set_label_coords(1.02,0.5)
    fg.savefig(fn,dpi=600)


def initial_topography(A,path_save):

    def fmt(x):
        s = f"{x:.1f}"
        if s.endswith("0"):
            s = f"{x:.0f}"
        return rf"{s} km" if plt.rcParams["text.usetex"] else f"{s} km"
        
    fg = figure()
        
    fn = os.path.join(path_save,'Figure_S2')
    ax = fg.gca()
    p1 = ax.contourf(A.C.xg,A.C.yg,A.FS.H[:,:,10],levels = [-3.5,-3.0,-2.5,-2.0,-1.5,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5],cmap = 'cmc.oleron',linewidths=0.5)
    #ax.clabel(p1, p1.levels, inline=True, fmt=fmt, fontsize=fnt_g.axis_)
    ax.plot(A.C.x_trench_p,A.C.y_trench_p,linewidth = 2.0,linestyle = 'dashdot',label = r'Slab position',color = 'firebrick')
    ax.set_xlabel(r'$x$/[km]',fontsize=fnt_g.label_)
    ax.set_ylabel(r'$y$/[km]',fontsize=fnt_g.label_)
    ax.legend(loc='upper right')
    ax.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    cbar0 = fg.colorbar(p1, ax=ax,orientation='horizontal',extend="both")
    cbar0.ax.tick_params(labelsize=fnt_g.legend_)
    cbar0.set_label(label=r'${{H}}$ /$\mathrm{km}$',size=fnt_g.label_) 
        
    fg.savefig(fn,dpi=600,transparent=False)


def exp_ernn(A,B,path_figure,figure_name):
    # figure name
    fn = os.path.join(path_figure,'%s.png'%(figure_name))
    
    # Prepare variables
    c0  = A.FS.time_series 
    c1  = B.FS.time_series
        
 
    


    

    i10 = ((A.time>1) & (A.time<np.nanmax(A.Det.det_vec[(A.C.x_sp>100) & (A.C.x_sp<1100)])))
    i11 = ((B.time>1) & (B.time<np.nanmax(B.Det.det_vec[(B.C.x_sp>100) & (B.C.x_sp<1100)])))
    t0 = A.time[i10==1]
    t1 = B.time[i11==1]
    t0 = t0[:-1]

    t1 = t1[:-1]
    
    
    string_title0 = r'Fast Tearing'
    string_title1 = r'Slow Tearing'
    # Prepare figure layout 
    cm = 1/2.54  # centimeters in inches
    fg = figure(figsize=(15*cm, 8*cm))  
    bx = 0.07
    by = 0.1
    sx = 0.45
    dx = 0.03
    sy = 0.7
    dy = []
    ax0 = fg.add_axes([bx, by, sx, sy])
    ax1 = fg.add_axes([bx+sx+dx, by, sx, sy])
    pl0 =  ax0.plot(t0/np.nanmax(A.Det.det_vec),c0[:,0],linewidth=1.3, color = 'firebrick',linestyle='dashed')
    pl0b = ax0.plot(t0/np.nanmax(A.Det.det_vec),c0[:,1],linewidth=1.3, color = 'firebrick')
    
       
 
    pl1   = ax1.plot(t1/np.nanmax(B.Det.det_vec),c1[:,0],linewidth=1.3, color = 'forestgreen',linestyle='dashed')
    pl1b  = ax1.plot(t1/np.nanmax(B.Det.det_vec),c1[:,1],linewidth=1.3, color = 'forestgreen')

    #ax0.set_ytick(0.1,0.5,0.9)
    ax0.tick_params(axis="y",direction="in")
    ax0.tick_params(axis="x",direction="in")
    ax1.tick_params(left=True,right=True,labelleft=False) 
    ax1.tick_params(axis="y",direction="in")
    ax1.tick_params(axis="x",direction="in")
    
    plt.setp(ax0.spines.values(), linewidth=1.4)
    plt.setp(ax1.spines.values(), linewidth=1.4)

    ax0.tick_params(width=1.2)
    ax1.tick_params(width=1.2)
    
    
    ax0.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax0.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax1.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax0.set_title(string_title0,fontsize=fnt_g.title_)    #ax0.yaxis.set_label
    ax1.set_title(string_title1,fontsize=fnt_g.title_)    #ax0.yaxis.set_label


    props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)
    ax0.text(0.87, 0.95, '$[a]$', transform=ax0.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.87, 0.94, '$[b]$', transform=ax1.transAxes, fontsize=fnt_g.label_,
        verticalalignment='top', bbox=props,color='white')




    fg.savefig(fn,dpi=600)

def make_figure_Sup(A,path_figure,figure_name,time):
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
    
    

    ax0=setup_axix(ax0,r'$\Psi$ / $\mathrm{W}/\mathrm{m}^3$')
    ax1=setup_axix(ax1,r'$T$ / $^{\circ}$C')
    ax2=setup_axix(ax2,r'$\tau^{\dagger}_{\mathrm{max}}$ / n.d.')
    it = 0 
    for i in ind0[0]:
        alpha = 0.1+(0.6/(len(ind0[0])-1))*it
        lw = 0.1+(0.6/(len(ind0[0])-1))*it
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
    ax2.set_xlabel(r'$x_{\mathrm{trench}}$ / km',fontsize=fnt_g.label_)


    props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)
    ax0.text(0.025, 1.15, '[a]', transform=ax0.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.025, 1.15, '[b]', transform=ax1.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
    ax2.text(0.025, 1.15, '[c]', transform=ax2.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
    
   

    #ax2.yaxis.set_label_coords(1.02,0.5)
    fg.savefig(fn,dpi=600)
    
def figure_experimental_supplementary(A,path_figure,figure_name):
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

    for i in range(len(time)):
        if time[i]>1.0:
            ind = np.where(time_1D[:,i]==np.max(time_1D[:,i]))
            ind = ind[0][0]
            max_up_x[i] = A.C.x_sp[ind]
    
    for i in range(len(time_pv)-1):
       
       
    
        t0 = time_pv[i]
        t1 = time_pv[i+1]
        t_buf = (time>=t0) & (time<t1)
        
        for ix in range(len(A.C.x_trench_p)):
            time_1D_d[ix,t_buf==1] = np.mean(time_1D[ix,t_buf==1])
            time_0D_d[t_buf==1]    = np.mean(time_0D[t_buf==1])
        
        
    
    cm = 1/2.54  # centimeters in inches

    fn = os.path.join(path_figure,'%s.png'%(figure_name))
    
    fg = figure(figsize=(14*cm, 12*cm))
    
    bx = 0.085
    by = 0.12
    dx = 0.12 
    dy = 0.02 
    sx = 0.35 
    sy1 = 0.8
    sy2 = 0.25
    # First panel with the pcolormap
    ax0 = fg.add_axes([bx, by, sx, sy1]) # t1 
    
    # 0D timeseries discrete
    ax1 = fg.add_axes([bx+sx+dx, by+2*sy2+2*dy, sx, sy2])           # t3 
    ax2 = fg.add_axes([bx+sx+dx, by+dy+sy2, sx, sy2])     # t2 
    ax3 = fg.add_axes([bx+sx+dx, by, sx, sy2])     # t2 
    
    p1 = ax0.pcolormesh(A.C.x_sp,time,np.transpose(time_1D_d[:-1, :-1]),shading='auto',cmap = 'cmc.imola',vmin = np.nanpercentile(time_1D_d[time_1D_d != 0.0],10), vmax=np.nanpercentile(time_1D_d[time_1D_d != 0.0],90))
    ax0.plot(max_up_x,time,color='red',linewidth = 2.0)
    #ax.clabel(p1, p1.levels, inline=True, fmt=fmt, fontsize=fnt_g.axis_)
    ax0.axvline(100,linewidth = 2.0,linestyle = 'dashdot',label = r'[b]',color = 'blue')
    ax0.axvline(500,linewidth = 2.0,linestyle = 'dashdot',label = r'[c]',color = 'forestgreen')
    ax0.axvline(1000,linewidth = 2.0,linestyle = 'dashdot',label = r'[d]',color = 'k')
    ax0.legend(loc='upper center', bbox_to_anchor=(0.45, 1.12),ncol=3, columnspacing=0.05,handletextpad=0.1, shadow=True,fontsize=8)
    ax0.set_xticks([100,500,1000])

    ax0.set_xlabel(r'$x_{trench}$/ km',fontsize=fnt_g.label_)
    ax0.set_ylabel(r'$time$/ Myr',fontsize=fnt_g.label_)
    ax0.set_ylim(1,time_max)
    ax0.xaxis.set_tick_params(labelsize=fnt_g.axis_)
    ax0.yaxis.set_tick_params(labelsize=fnt_g.axis_)
    cbar0 = fg.colorbar(p1, ax=ax0,orientation='horizontal',extend="both")
    cbar0.ax.tick_params(labelsize=fnt_g.legend_)
    cbar0.set_label(label=r'${{\dot{H}}}$ / $\mathrm{mm/yr}$',size=fnt_g.label_) 
    
    ia = np.where(A.C.x_sp>=100)
    ia = ia[0][0]
    
    ib = np.where(A.C.x_sp>=500)
    ib = ib[0][0]
    
    ic = np.where(A.C.x_sp>=1000)
    ic = ic[0][0]
    
    pb=ax1.plot(time,time_1D_d[ia,:],color = 'salmon',linewidth=1.2,label='Filtered Data')
    pb1 = ax1.plot(time,time_1D[ia,:],color = 'indigo',linewidth=0.7,label='Raw Data')
    ax1.legend(loc='upper center', bbox_to_anchor=(0.45, 1.25),ncol=2, columnspacing=0.05,handletextpad=0.1, shadow=True,fontsize=8)

    #pb2 = ax1.plot(time,time_0D_d,color='gainsboro',linewidth=0.6,label='Average')
    ax1.set_xlim(1,time_max)
    ax1.set_ylim(np.nanmin(time_1D_d[time_1D_d != 0.0]),np.nanmax(time_1D_d[time_1D_d != 0.0]))
    
    pc=ax2.plot(time,time_1D_d[ib,:],color = 'salmon',linewidth=1.2,label='Filtered Data')
    pc1 = ax2.plot(time,time_1D[ib,:],color = 'indigo',linewidth=0.7,label='Raw Data')
    #pc2 = ax2.plot(time,time_0D_d,color='gainsboro',linewidth=0.6,label='Average')
    ax2.set_xlim(1,time_max)
    ax2.set_ylim(np.nanmin(time_1D_d[time_1D_d != 0.0]),np.nanmax(time_1D_d[time_1D_d != 0.0]))

    pd=ax3.plot(time,time_1D_d[ic,:],color = 'salmon',linewidth=1.2,label='Filtered Data')
    pd1 = ax3.plot(time,time_1D[ic,:],color = 'indigo',linewidth=0.7,label='Raw Data')
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
    
    ax3.set_xlabel(r'$time$/ Myr',fontsize=fnt_g.label_)
    ax1.set_ylabel(r'$\dot{H}$/ mm/yr',fontsize=fnt_g.label_)
    ax2.set_ylabel(r'$\dot{H}$/ mm/yr',fontsize=fnt_g.label_)
    ax3.set_ylabel(r'$\dot{H}$/ mm/yr',fontsize=fnt_g.label_)




    props = dict(boxstyle='round', facecolor='black',edgecolor='none', alpha=0.8)
    ax0.text(0.025, 0.95, '[a]', transform=ax0.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
    ax1.text(0.025, 0.95, '[b]', transform=ax1.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
    ax2.text(0.025, 0.95, '[c]', transform=ax2.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
    ax3.text(0.025, 0.95, '[d]', transform=ax3.transAxes, fontsize=fnt_g.legend_,
        verticalalignment='top', bbox=props,color='white')
   

    
    fg.savefig(fn,dpi=600)

