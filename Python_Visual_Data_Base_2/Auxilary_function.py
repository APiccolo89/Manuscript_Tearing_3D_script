
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
    plt.close()


def  _scatter_plot_(Data:Data_Base,path_save:str,label_scatter:label_scatter,fields:list,name_figure,stress_limit):
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
    
    if stress_limit == 'PR_r':
        P = 400e6
    elif stress_limit == 'PR':
        P =200e6 
    else:
        P = 600e6
    
    z =  eval(z_f,globals(),Data.__dict__)
    m =  eval(m_f,globals(),Data.__dict__)
    p = Data.StressLimit
    m_u = np.unique(m) 
    sp = plt.scatter(x[(m==P)],y[(m==P)],60,z[(m==P)],marker=label_scatter.markers[0],cmap = label_scatter.colormap,edgecolors='k')
    cbar = fg.colorbar(sp,orientation='horizontal',extend="both",label=label_scatter.cbar_label)
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
    name_fig = '%s%s.png' %(name_figure,stress_limit)
    fn = os.path.join(path_save,name_fig)  
    #plt.show()
            
    fg.savefig(fn,dpi=300)
    plt.close()
    
    
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
    plt.close()
    
    
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
