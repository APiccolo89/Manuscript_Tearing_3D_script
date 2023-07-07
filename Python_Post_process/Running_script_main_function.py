
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 09:07:30 2022

@author: Andrea Piccolo
Slab Avalanching small project

1) Function that reads the vtk file and produce picture
    a) Read vtk, plot phase, topography 

"""

import sys,os
import numpy as np
from time import perf_counter 
import getopt
import argparse
    

class SLAB():
    def __init__(self,C,nstep):
        tz = np.sum(C.ind_z)
        self.D     = np.ones((tz,nstep),dtype=float)*(-1.0)
        self.T     = np.zeros((tz,nstep),dtype=float)
        self.nu    = np.zeros((tz,nstep),dtype=float)
        self.eps   = np.zeros((tz,nstep),dtype=float)
        self.tau   = np.zeros((tz,nstep),dtype=float)
        self.vz    = np.zeros((tz,nstep),dtype=float)
        self.vx    = np.zeros((tz,nstep),dtype=float)
        self.vis   = np.zeros((tz,nstep),dtype=float)
        self.F_T   = np.zeros((tz,nstep),dtype=float)
        self.F_B   = np.zeros((tz,nstep),dtype=float)
        self.dRho  = np.zeros((tz,nstep),dtype=float)
        self.vis_A    = np.zeros((tz,nstep),dtype=float)
        self.nu_A    = np.zeros((tz,nstep),dtype=float)
        self.det_vec =np.zeros(nstep)
        self.tau_det =np.zeros(nstep)
        self.D_det   = np.zeros(nstep)
        self.vz_S    = np.zeros(nstep)
        self.x1    = np.zeros(tz)
        self.x2    = np.zeros(tz)   
        self.LGV         = ["D","T","eps","tau","vz","vx","vis","nu","F_T"]
        self.Label       = ["$\delta_{ap} [km]$",     
                            "T [$^{\circ}C$]",
                            "$log_{10}(\epsilon_{II}) [1/s]$",
                            "$\tau_{II} [MPa]$",
                            "$v_x [cm/yrs]$",
                            "$v_z [cm/yrs]$",
                            "$\eta_{vep} [Pa\cdot s]$",
                            "$\eta_{creep} [Pa\cdot s]$",
                            "\u0394 \u03C1 $[kg/m^3]$",
                            "$\dot{\epsilon_{zz}} [1/s]$",
                            r"$\tau_{zz} [MPa]$",
                            r"F($\tau_{zz}$) [N/m]"
                            ]
        self.Colormap    = ["cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.oleron","cmc.oleron","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao"]
        self.Val         = [("min","max"),
                            ("min","max"),
                            ("min","max"),
                            ("min","max"),
                            ("min","max"),
                            ("min","max"),
                            ("min","max"),
                            ("min","max"),
                            ("min","max"),
                            ("min","max"),
                            ("min","max"),
                            ("min","max"),
                            ("min","max"),
                            ("min","max"),
                            ("min","max"),
                            ("min","max")]
        self.CV = ["T","nu","vis","eps","tau","vz","vx"]        


    def _update_C(self,C,Values,ipic,time):         
        
        self = self._Find_Slab_C(C,Values,ipic)
       
                
    def _Find_Slab_C(self,C,Values,ipic):
        
            """
            Better writing a bit: Now, at the beginning I was tayloring too much the code 
            with the problem at hand, as a consequence, I was not considering too much 
            different scenarios:
                Area of investigation -100 km to Length slab - 10 km (arbitrary measure)
                Everything that is outside this area is arbitrary set to be equal to -1.0 and infinite
                The moment in which the slab is detached is signalled by nan
            1) The continental lithosphere is accounted, as it deforms during the evolution of the experiments
            which would highlighting a detachment that does not exist. 
            ALG : 
                -> Everyvalues of the class set to -inf if outside the area of research 
                -> Compute the lithospheric budged (buf = OP+C1+C2)
                -> Find where buf >0.9 
                -> Compute the average and the thickness of the slab
                -> if D < 5 km -> effectively detached (5/80 5/(2*(40)) -> 1/16 of the origianl thickness): 
                        Most of the time the slab has non zero thickness due to presence of small patches of lithosphere
                        apperently this 5 km values worked well allmost for all the application
                
                1) LOOP over z 
                    if inside 
                        compute buf 
                        if exist slab
                        compute the averages
                        compute the thickness
                        if the thickness is less than 5 km
                        set everything to nan; 
                    if outside 
                        set everything to - inf 
            """
        
        
            tz    = np.sum(C.ind_z)
        
            self.x1    = np.zeros(tz)
            self.x2    = np.zeros(tz)
           
            z     = C.z 
            x     = C.x

            for i in range(tz):
         
                if(z[i]>-100):
                    buf = -1.0
                else:
                    buf = Values.OP[i,:]+Values.C1[i,:]+Values.C2[i,:]
                    buf[buf<0.95] = -1.0
                    if len(buf[buf ==-1.0]) == len(buf):
                        buf = -1.0 
                ind = np.where(buf>0.95)
                if(np.size(ind)==0 & (np.size(buf)>1)):
                    self.x2[i]=np.nan
                    self.x1[i]=np.nan
                    """
                    Compute the mean and standard deviation of certain value
                    """                    
                    CV = ["T","nu","vis","eps","epsxx","epszz","epsxz","tau","tauxx","tauzz","tauxz","vz","vx"]
                    for iv in self.CV:
                        eval(iv,globals(),self.__dict__)[i,ipic] = np.nan 
                        
                    self.nu_A[i,ipic] = np.nan
                    self.vis_A[i,ipic] = np.nan
                    self.dRho[i,ipic] = np.nan
                    self.D[i,ipic] = np.nan
                    if ((np.size(buf) == 1)):
                        self.x2[i]=- float("inf")
                        self.x1[i]=- float("inf")
                        """
                        Compute the mean and standard deviation of certain value
                        """                    
                        CV = ["T","nu","vis","eps","epsxx","epszz","epsxz","tau","tauxx","tauzz","tauxz","vz","vx"]
                        for iv in self.CV:
                            eval(iv,globals(),self.__dict__)[i,ipic] = - float("inf") 
                            
                        self.nu_A[i,ipic] = - float("inf")
                        self.vis_A[i,ipic] = - float("inf")
                        self.dRho[i,ipic] = - float("inf")
                        self.D[i,ipic] = - float("inf")
                    
                else:
                    self.x1[i]= x[ind[0][0]]
                    self.x2[i] = x[ind[0][-1]]
                    self.D[i,ipic] = self.x2[i]-self.x1[i] 
                    if(self.D[i,ipic] <8):
                        for iv in self.CV:
                            eval(iv,globals(),self.__dict__)[i,ipic] = np.nan 
                        self.nu_A[i,ipic] = np.nan
                        self.vis_A[i,ipic] = np.nan
                        self.dRho[i,ipic] = np.nan
                        self.D[i,ipic] = np.nan
                    else:
                        for iv in self.CV:
                            eval(iv,globals(),self.__dict__)[i,ipic] = np.mean(eval(iv,globals(),Values.__dict__)[i,ind])
                        ind_AA = (x>=self.x1[i]-60) & (x<self.x1[i])
                        ind_BB = ((x>self.x2[i]) & (x<=self.x2[i]+60))
                        ind_A   = np.where((ind_AA == True)| ind_BB == True)
                        self.nu_A[i,ipic] = np.mean(Values.nu[i,ind_A])
                        self.vis_A[i,ipic] = np.mean(Values.vis[i,ind_A])
                        ind2 = np.where(buf==0.0)
                        Rho_AST = np.mean(Values.Rho[i,ind2])
                        Rho_Slab = np.mean(Values.Rho[i,ind])
                        self.dRho[i,ipic] = Rho_AST-Rho_Slab
                        
            L_id = np.where(np.isnan(self.x1)==False) 
            z_b  = z[L_id[0][0]]
            self.vz_S[ipic] = np.nanmean(Values.vz[Values.OP>0.9])

            dRho = -np.nanmean(self.dRho[:,ipic])
            self.F_T[:,ipic] = 2*self.tau[:,ipic]*self.D[:,ipic]*1e9
            self.F_B[:,ipic] = dRho*9.81*self.D[:,ipic]*(z-z_b)*1e6                
            return self 
    
    def  _plot_average_C(self,t_cur,z,ptsave,ipic): 

        
        time_sim = "{:.3f}".format(t_cur)
        
        ic  = 0
        
        fna='Fig'+str(ipic)+'.png'
        ptsave_b=os.path.join(ptsave,'Averages')
        if not os.path.isdir(ptsave_b):
            os.mkdir(ptsave_b)

       
        var =self.LGV
        index = self.Label 

        for values in var: 

       
            fg = figure()
    
            tick='%s Time = %s Myrs' %(index[ic],time_sim)

            ptsave_c=os.path.join(ptsave_b,values)
            if not os.path.isdir(ptsave_c):
                os.mkdir(ptsave_c)
        
            fn = os.path.join(ptsave_c,fna)

            ax0 = fg.gca()
            if(values == "dF"):
                ax0.plot(self.F_T[:,ipic]/self.F_B[:,ipic],z, lw=1.2,c='k',ls='-')
            else:
                ax0.plot(eval(values,globals(),self.__dict__)[:,ipic],z, lw=1.2,c='k',ls='-')
            plt.grid(True)
            if((values == "eps") | (values == "epsxx") | (values == "epszz")| (values == "epsxz")):
                ax0.set_xscale("log")
            ax0.set_title(tick)
            ax0.tick_params(axis='both', which='major', labelsize=5)
            ax0.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
            plt.draw()    # necessary to render figure before saving
            fg.savefig(fn,dpi=300)
            ax0.plot()
            ic += 1
            val = [] 
            
            plt.close()
        
    def  _plot_slab_time(self,time,z,ptsave,TestName,Data_BASE_F,t_lim): 

        import cmcrameri.cm as cmc

        
        ic  = 0
        
        ptsave_b=os.path.join(ptsave,'Averages')
        if not os.path.isdir(ptsave_b):
            os.mkdir(ptsave_b)

        var =self.LGV
        index = self.Label   
        cmaps = self.Colormap
        LIM = self.Val
        self.det_vec, ind_z,TT,ZZ = self._find_detachment_age(time,z)
        self.tau_vec = self.tau[ind_z,:]
        self.D_det   = self.D[ind_z,:]
#        if len(self.det_vec) >0:
#            self._save_txt_DATABASE(self.det_vec,TestName,Data_BASE_F)
#        

        for values in var: 
            print(values)
            fna='Fig%s.png' %values
            if(values == "dF"):
                val = self.F_T/self.F_B
            else:
                val = eval(values,globals(),self.__dict__)

            
            if values == "D":
                val = val/80 

            
            cmap2 = eval(cmaps[ic])
            lm    = LIM[ic]
            if len(self.det_vec)>0:
                buf = val[ind_z,:]
                buf[buf==0.0] = np.nan
                lim_m = np.nanmin(buf)
                lim_M = np.nanmax(buf)
                buf =[]
            elif((lm[0]=="min") & len(self.det_vec)==0):
                lim_m = np.nanmin(val)
                lim_M = np.nanmax(val)
            else:
                lim_m = lm[0]
                lim_M = lm[1]
            LCont = np.linspace(lim_m,lim_M,num=5)

            if((values == "eps")| (values == "Psi")):
                val = np.log10(val)
                lim_m = np.log10(lim_m)
                lim_M = np.log10(lim_M)
            fg = figure()
    
            tick='%s' %(index[ic])

            ptsave_c=os.path.join(ptsave_b,values)
            if not os.path.isdir(ptsave_c):
                os.mkdir(ptsave_c)
        
            fn = os.path.join(ptsave_c,fna)

            ax0 = fg.gca()

            levels = np.linspace(lim_m,lim_M,num=10)

            cf=ax0.contourf(time,z,val, cmap=cmap2,levels = levels)
            if len(self.det_vec)>0:
                ax0.scatter(self.det_vec[0],self.det_vec[1],s = 50,marker ="X",color = "r")
                ax0.plot(TT,ZZ,lw=1.2,c='k',ls='--')
           
            ax0.set_ylim(-600,-100)
            ax0.set_xlim(0,t_lim)
            plt.grid(True)
            cbar = fg.colorbar(cf, ax=ax0,orientation='horizontal')
            cbar.set_label(index[ic])
            plt.xlabel('$time, [Myrs]$')
            plt.ylabel('$Depth, [km]$')
 
            ax0.set_title(tick)
            ax0.tick_params(axis='both', which='major', labelsize=5)
            ax0.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
            plt.draw()    # necessary to render figure before saving
            fg.savefig(fn,dpi=300)
            ax0.plot()
            ic += 1
            val = [] 
            
            plt.close()
        return ind_z
    def _find_detachment_age(self,time,z):
        
        """
        Long story short: I start the loop from the top to the bottom, 
        the first occurence of nan is defined as effective detachment, 
        this value has a bit of error attacched, but i figure out, one day
        a better cut off. 
        I select all the point along the time vector that fullfill the detachment
        condition and i took the lowest one and the associated depth and i save
        the detachment depth, the detachment time, and the value that are associated 
        with the previous time step avaiable. Took allmost 6 hrs to find out, better to 
        associated this part of the script with a lot of shame and terrible comment on my
        self. 
        ALG: 
            start loop from the last element 
            find where D is NAN and save time and space coordinate per each 
            append the indexes 
            
        check if exist any element in the vector created before
        if yes, check what is the minimum time in which nan appear 
        then save time, save the depth 
        retrieve the 1D vector for each of the value that have been computed during the post processing
        """
            
        det_vec =[]                             #Vector that collects the the value at detachment time
        ind_z = []                              #empty vector for ind_z
        lz = len(z)                             #length vector lz
        chos_x = []                             #index of the chosen 
        chos_T = []
        for i in range(lz-1,0,-1):              #Loop starting from the top to the bottom 
            buf = self.D[i,:]                   #1D(z) vector of D 
            if (z[i] <-100) & (z[i]>-350):      #Constraining better the area of research
                ind_T = np.where(np.isnan(buf)) #find where in buf (1D vec D) is nan
                if np.size(ind_T)>0:            #if the vector has a size
                    chos_x.append(i)            #append the index i 
                    chos_T.append(ind_T[0][0])  #append the index associated with the first nan (left->right)
        if len(chos_T) > 0:                     #AFTER LOOP. check chos_T, if chos_T has a at least one element
            TD = np.min(time[chos_T])           #Time of the detachment is the min within the time in which is detected detachment
            buf2 = time[chos_T]                 #collect all the time points associated time 
            buf3 = z[chos_x]                    #collect the z vector to produce t-z path of the detachment 
            i_buf   = np.where(buf2==TD)        #find the point of the detachment
            iz      = chos_x[i_buf[0][0]]       #find the point z of the detachment
            depth   = z[iz]                     #save the depth
            idt     = np.where(time==TD)        #save the time
            idT     = idt[0][0]                 #save the the index
             
            time_ = time[idT-1]                 #save time step before of the detachment, and relative value
            
            tau_  = self.tau[iz,idT-1]
            eps_  = self.eps[iz,idT-1]
            T_    = self.T[iz,idT-1]
            vis_  = self.vis[iz,idT-1]
            nu_   = self.nu[iz,idT-1]
            F_B   = self.F_B[iz,idT-1]
            F_T   = self.F_T[iz,idT-1]
            nu_A  = self.nu_A[iz,idT-1]
            vis_A = self.vis_A[iz,idT-1]
            det_vec = [TD, depth,time_,vis_,nu_,tau_,eps_,F_T,F_B,T_,vis_A,nu_A]
        else: 
            det_vec=[]
            iz= []
            buf2 =[]
            buf3 = []
        
               
        return det_vec, iz,buf2,buf3 
    
    def _save_txt_DATABASE(self,det_vec,TestName,ptsave):
        
        import csv 
        
        file_name = "Test_Detachment_Data_BASE.csv"
        
        """
        Write csv data files with the detachment vector 
        """
       
    
        filename = os.path.join(ptsave,file_name)
        if(os.path.isfile(filename)==False):
            header = ["TestName", "timedet","Depth","timets","effectiveviscosity","creepviscosity","tauII","tauxx","tauzz","tauxz","epsII","epsxx","epszz","epsxz","F_T","F_B","Temp","viscoeffAst","viscocreepAst"]
            f = open(filename, 'a')
            writer = csv.writer(f)

            writer.writerow(header) 
            f.close()
        f = open(filename, 'a')
        writer = csv.writer(f)
        row    = [TestName]
        for i in range(len(det_vec)):
            row.append(det_vec[i])
            
        writer.writerow(row)

        # close the file
        f.close()




"""
Import useful function 
"""
from Read_VTK_files_LAMEM import * #Read_VTK_files_LAMEM

from Read_VTK_files_LAMEM import  _file_list

from Initial_Input import * 

def _run_script_visualization(ptsave,Folder,Test_Name,l_path):
    t_INIZIO = perf_counter()

    
    dic_val= {
    "OP"       : "_Oceanic_Lit_ [ ]",
    "CC1"       : "_Continental_Crust_ [ ]",
    "T"        : "temperature [C]",
    "velocity" : "velocity [cm/yr]",
    "disp"     : "tot_displ [km]",
    "vis"      : "visc_total [Pa*s]",
    "nu"       : "visc_creep [Pa*s]",
    "eps"      : "j2_strain_rate [1/s]",
    "tau"      : "j2_dev_stress [MPa]",
    "Rho"      : "density [kg/m^3]",
    "gamma"    : "plast_strain [ ]"}
    
    phase_dictionary={
     "0" : ["Air","white"],
     "1" : ["Upper Crust", "silver"],
     "2" : ["Lower Crust", "dimgray"],
     "3" : ["Continental Lithosperhic mantle","lightgoldenrodyellow"],
     "4" : ["Continental Lithosperhic mantle2","palegoldenrod"],
     "5" : ["Upper Mantle", "aliceblue"],
     "6" : ["Slab","sandybrown"],
     "7" : ["Oceanic Crust","palegreen"],
     "8" : ["Weak Zone Mantle","rosybrown"],
     "9" : ["Upper Crust2","linen"],
     "10" : ["Lower  Crust 2","gainsboro"],
     "11" : ["Pelagic Sediment","lavender"],
     "12" : ["Prism","burlywood"],
     "13" : ["Passive Margin","darkkhaki"],
     "14" : ["Mantle 410 km","palegreen"],
     "15" : ["Slab 410 km","peachpuff"],
     "16" : ["Eclogite","palevioletred"],
     "17" : ["Lower  Mantle","darkseagreen"],
     "18" : ["Lower Slab","tan"]
     }
    
    fname=("SDet")
    dyn             ='%s.pvtr'      %fname
    surf            ='%s_surf.pvts' %fname
    phase           ='%s_phase.pvtr'%fname
    #passive_tracers ='%s_passive_tracers.pvtu' %fname
    ###############################################
    fname2 = fname+".pvd"
    fname=os.path.join(Folder,Test_Name,fname2)
    time, Flist, n_tim_step =_file_list(fname)    # Retrive the file list and the time step information
    ##############################################
    ts0 = Flist[0]
    fn=ts0[1:ts0.find('/')]
    Filename_0 = [os.path.join(Folder,Test_Name,fn,dyn),os.path.join(Folder,Test_Name,fn,phase)]
    folder_ = os.path.join(Folder,Test_Name)
    C = Coordinate_System(Filename_0,ptsave,(-600.0,600.0),(-600.0,50.0))
    Phase_DB = Phase_Data_Base(folder_)
    Initial_Condition = Initial_condition(folder_,Test_Name)
    Initial_Condition.tc = Initial_Condition.tc/3.5
    ###############################################
    FSurf = FS(C,len(time)) # Create the instance of free surface class 
    DYN   = VAL(C,dic_val)  # Create the instance of the .pvd file 
    Ph    = Phase(C,phase_dictionary) # Create the instance of the Phase field
    Slab  = SLAB(C,len(time))        # Create the instance of the Slab. 
    ################################################
    ###############################################################################
    # Read Information related to the initial condition
    ###############################################################################
    ipic= 0 
    for istp in Flist:
    ########################### Files and time stepping information############
        t1_start = perf_counter()
        fn=istp[1:istp.find('/')]
        t_cur=time[ipic]        
        Filename_dyn=os.path.join(Folder,Test_Name,fn,dyn)
        Filename_ph=os.path.join(Folder,Test_Name,fn,phase)
        Filename_s=os.path.join(Folder,Test_Name,fn,surf)
    #    Filename_ptr=os.path.join(Folder,Test_Name,fn,passive_tracers)
        
        ###### Retrieve the field that is needed ##################################
        t1 = perf_counter()
        DYN._update_Val(Filename_dyn,C)
        t2 = perf_counter()
        Values_time = t2-t1
        print("Dynamic value took","{:02}".format(Values_time))
        ###########################################################################
        t1 = perf_counter()
        Ph._update_phase(Filename_ph,C)  
        t2 = perf_counter()    
        average_plot_time = t2-t1
        print("Phase_update","{:02}".format(average_plot_time))

        ###########################################################################
        t1= perf_counter()
        FSurf._Update_(Filename_s,C,ipic)
        FSurf.ASCI_FILE(ipic,t_cur,Test_Name,ptsave,C.x)
        t2 = perf_counter()
        print("Free surface ","{:02}".format(t2-t1))

        ###########################################################################
        t1 = perf_counter()
        Slab. _update_C(C,DYN,ipic,time)
        t2 = perf_counter()
        print("Slab routine ","{:02}".format(t2-t1))
        ###########################################################################
        t1 = perf_counter()   
        if (ipic % 10000 == 0):
            # Plot the field of interest of dynamic 
            DYN._plot_maps_V(t_cur,C.z,C.x,ptsave,ipic)
            # Plot the phase plot 
            Ph._plot_phase_field(C, DYN, ptsave, ipic, t_cur,FSurf)
            # Plot the Slab averages 
            Slab._plot_average_C(t_cur,C.z,ptsave,ipic)          
            t2 = perf_counter()
        plot_map_time = t2-t1
        ###########################################################################
        #  Free Surface 
        #######################################################################
          
        FS_time = t2-t1
        t1_end = perf_counter()
        tstep_time_pr = t1_end-t1_start
        minutes, seconds = divmod(t1_end-t1_start, 60)
          
        ipic +=1 
        
        
        
        
        t_FINE = perf_counter()      
        
        minutes, seconds = divmod(t_FINE-t_INIZIO, 60)
        print("===========================================")
        
        print("||Script took ","{:02}".format(int(minutes)),"minutes and","{:05.2f}".format(seconds),' seconds ||' )
        
        print("===========================================")
    
    ind_z_t = Slab._plot_slab_time(time,C.z,ptsave,Test_Name,ptsave,60)


