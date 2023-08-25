# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 16:32:06 2023

@author: Andrea Piccolo
"""

import sys,os
import numpy as np
from time import perf_counter 
import getopt
import argparse
import matplotlib
import numpy as np 
import h5py 
from numba import jit
from numba import jitclass, types, typed, prange
from Read_VTK_files_LAMEM import * #Read_VTK_files_LAMEM

from Read_VTK_files_LAMEM import  _file_list
from Read_VTK_files_LAMEM import findnodes
from Read_VTK_files_LAMEM import bilinearinterpolation
from Read_VTK_files_LAMEM import trilinearinterpolation
from Read_VTK_files_LAMEM import linearinterpolation


from Parser_File import * 

    

class SLAB():
    def __init__(self,C,IG,nstep):
        Slab_Geometry = IG.Slab
        self.P1 = Slab_Geometry.Boundary_B[0][0]
        self.P2 = Slab_Geometry.Boundary_B[0][2]
        t= (C.xp>=self.P1)&(C.xp<=self.P2)
        self.ind_boundary = np.where(t==True)
        tx  = len(self.ind_boundary[0])

        
        tz = len(C.zp)
        self.D     = np.ones((tx,tz,nstep),dtype=float)*(-1.0)
        self.T     = np.zeros((tx,tz,nstep),dtype=float)
        self.nu    = np.zeros((tx,tz,nstep),dtype=float)
        self.eps   = np.zeros((tx,tz,nstep),dtype=float)
        self.tau   = np.zeros((tx,tz,nstep),dtype=float)
        self.vis   = np.zeros((tx,tz,nstep),dtype=float)
        self.F_T   = np.zeros((tx,tz,nstep),dtype=float)
        self.F_B   = np.zeros((tx,tz,nstep),dtype=float)
        self.Rho  = np.zeros((tx,tz,nstep),dtype=float)
        self.Psi = np.zeros((tx,tz,nstep),dtype=float)
        self.L   = np.zeros((tx,tz,nstep),dtype=float)
        self.det_vec =np.zeros((tx,nstep),dtype=float)
        self.x1    = np.zeros((tx,tz,nstep),dtype=float)
        self.x2    = np.zeros((tx,tz,nstep),dtype=float)   
        self.LGV         = ["D","T","tau","F_B","F_T",'eps','Psi']
        self.Label       = ["$D^{\dagger} []$",     
                            "T [$^{\circ}C$]",
                            "$\\tau^{\dagger}_{II} []$$",
                            "$F_B [\\frac{N}{m}]$",
                            "$2 \cdot D \\tau [\\frac{N}{m}]$",
                            "$\dot{\epsilon}^{\dagger}_{II} []$",
                            "$\Psi [\\frac{W}{m^3}]$"
                            ]
        self.Colormap    = ["cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.oleron","cmc.oleron","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao"]
        self.Val         = [(0.1,0.85),
                            (700,1200),
                            ("min","max"),
                            (5e12,1e13),
                            (5e12,1e13),
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
        self.CV = ["T","nu","vis","eps","tau",'Rho','Psi']       


    def _update_C(self,C,FS,Ph,IG,ipic,time):
        # Create an array and select the phase that belongs to slab
        lay_ph = [Ph.Phase][0]
        lay_ph = lay_ph.astype('float64') # convert the type (phase field is uint8)
        ip = len(IG.Slab.Phases)
        # For whatever reason the type has been change, and the 1000 is 232 in uint8
        # Loop over the phases that belongs to slab and change the value
        for ic in range(ip):
            lay_ph[lay_ph==IG.Slab.Phases[ic]]=(1000)
        # Create the array to identify the slab
        lay_ph[lay_ph<(1000)] = -1.0
        lay_ph[lay_ph==(1000)] = 1.0 
        t1 = perf_counter()
        # loop over the x nodes of the margin
        for i in range(len(self.ind_boundary[0])):
            # Select the the actual node within x array Ph[k,j,ix],where ix node that belongs to interval xa-xb
            ix = self.ind_boundary[0][i]
            # Loop over the variable that you want to see in 2D plot of the slab
            for iv in self.CV:
                # Prepare buffer variable to fill with the new one
                buf_var_ph = eval(iv,globals(),Ph.__dict__)[:,:,ix]
                D     = np.zeros((np.sum(C.ind_zp)),dtype=float)
                buf_var   = np.zeros(np.sum(C.ind_zp),dtype=float)
                z_bottom   = np.zeros(np.sum(C.ind_zp),dtype=float)
                x1    = np.zeros(np.sum(C.ind_zp))
                x2    = np.zeros(np.sum(C.ind_zp))
                # Run the function 
                D,L0,buf_var,x1,x2 = _Find_Slab_PERFORM_C(C.yp,C.zp,buf_var_ph,lay_ph[:,:,ix],np.sum(C.ind_zp),ix,D,buf_var,x1,x2,z_bottom)
                # update the class 
                self.D[i,:,ipic] = D
                self.L[i,:,ipic] =L0
                eval(iv,globals(),self.__dict__)[i,:,ipic] = buf_var
                self.x1[i,:,ipic] = x1
                self.x2[i,:,ipic] =x2
            # Compute the additional variable (i.e. F_T = 2*D*tau), F_B
            self.F_T[i,:,ipic] = 2*self.D[i,:,ipic]*1e3*self.tau[i,:,ipic]*1e6
            self.F_B[i,:,ipic] = self.D[i,:,ipic]*1e3*self.L[i,:,ipic]*(self.Rho[i,:,ipic]-3300*(1-3e-5*1325))*9.81
            
        t2 =perf_counter()    
             
        print('time find_slab',"{:02}".format(t2-t1), 's')

                
    def _Find_Slab_C(self,C,Values,ph,ipic,ix,ix1):
        
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
        
        
            tz    = np.sum(C.ind_zp)
        
            self.x1    = np.zeros(tz)
            self.x2    = np.zeros(tz)
           
            z     = C.zp 
            x     = C.yp

            for i in range(tz):
         
                if(z[i]>-80):
                    buf = -1.0
                else:
                    buf = ph[i,:,ix]
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
                    for iv in self.CV:
                        eval(iv,globals(),self.__dict__)[ix1,i,ipic] = np.nan 
                        
                    self.nu_A[ix1,i,ipic] = np.nan
                    self.vis_A[ix1,i,ipic] = np.nan
                    self.dRho[ix1,i,ipic] = np.nan
                    self.D[ix1,i,ipic] = np.nan
                    if ((np.size(buf) == 1)):
                        self.x2[i]=- float("inf")
                        self.x1[i]=- float("inf")
                        """
                        Compute the mean and standard deviation of certain value
                        """                    
                        for iv in self.CV:
                            eval(iv,globals(),self.__dict__)[ix1,i,ipic] = - float("inf") 
                            
                        self.nu_A[ix1,i,ipic] = - float("inf")
                        self.vis_A[ix1,i,ipic] = - float("inf")
                        self.dRho[ix1,i,ipic] = - float("inf")
                        self.D[ix1,i,ipic] = - float("inf")
                    
                else:
                    self.x1[i]= x[ind[0][0]]
                    self.x2[i] = x[ind[0][-1]]
                    self.D[ix1,i,ipic]= self.x2[i]-self.x1[i] 
                    if(self.D[ix1,i,ipic] <8):
                        for iv in self.CV:
                            if (iv != 'F_T') & (iv != 'F_B'):
                                eval(iv,globals(),self.__dict__)[ix1,i,ipic] = np.nan 
                        self.nu_A[ix1,i,ipic] = np.nan
                        self.vis_A[ix1,i,ipic]= np.nan
                        self.dRho[ix1,i,ipic] = np.nan
                        self.D[ix1,i,ipic] = np.nan
                    else:
                        for iv in self.CV:
                            if iv == "tau_max":
                                eval(iv,globals(),self.__dict__)[ix1,i,ipic] = np.max(Values.tau[i,ind,ix])
                            elif iv == "tau_min":
                                eval(iv,globals(),self.__dict__)[ix1,i,ipic] = np.min(Values.tau[i,ind,ix])
                            else: 
                                if (iv != 'F_T') & (iv != 'F_B'):
                                    eval(iv,globals(),self.__dict__)[ix1,i,ipic] = np.mean(eval(iv,globals(),Values.__dict__)[i,ind,ix])
                            
                        ind_AA = (x>=self.x1[i]-60) & (x<self.x1[i])
                        ind_BB = ((x>self.x2[i]) & (x<=self.x2[i]+60))
                        ind_A   = np.where((ind_AA == True)| ind_BB == True)
                        self.nu_A[ix1,i,ipic] = np.mean(Values.nu[i,ind_A])
                        self.vis_A[ix1,i,ipic] = np.mean(Values.vis[i,ind_A])
                        Rho_AST = np.mean(Values.Rho[i,ind_A,ix])
                        Rho_Slab = np.mean(Values.Rho[i,ind,ix])
                        self.dRho[ix1,i,ipic] = Rho_AST-Rho_Slab
                        
            L_id = np.where(np.isnan(self.x1)==False) 
            z_b  = z[L_id[0][0]]
            temporary_drho = [self.dRho[ix1,:,ipic]][0]
            temporary_drho[temporary_drho==-np.inf]=np.nan
            dRho = -np.nanmean(temporary_drho)
            self.F_T[ix1,:,ipic] = 2*self.tau[ix1,:,ipic]*self.D[ix1,:,ipic]*1e9
            self.F_B[ix1,:,ipic] = dRho*9.81*self.D[ix1,:,ipic]*(z-z_b)*1e6                
            return self 
    
    def  _plot_average_C(self,t_cur,x,z,ptsave,ipic,Slab_Geo,IC): 

        x = x[self.ind_boundary]
        # Compute the arclength as a function of x
        boundary_geometry = Slab_Geo.Boundary_B[2]
        x_a = Slab_Geo.Boundary_B[0][0]
        y_a = Slab_Geo.Boundary_B[0][1]
        c   = boundary_geometry[0]
        cy  = boundary_geometry[1]
        R   = boundary_geometry[2]
        c_  = (x_a-c)**2-R**2+y_a**2
        b_ = 2.*y_a
        delta = np.sqrt(b_**2-4*c_)
        center_y = [(b_-delta)/2,(b_+delta)/2]
        y_c = np.min(center_y)
        y = (R**2-(x-c)**2)**(0.5)+y_c
        d = np.sqrt((x-x_a)**2+(y-y_a)**2)
        theta = 2*np.arcsin(d/(2*R))
        x_s = R*theta; 

        # => compute y with the data of the boundary
        # => compute the arc-length 
        # => create new vector 
        
        t_dim_less =  t_cur/(IC.tc/3.5) 
        time_sim = "{:.3f}".format(t_cur)
        time_dimen = "{:.3f}".format(t_dim_less[0][0])
        
        
        fna='Fig'+str(ipic)+'.png'
        ptsave_b=os.path.join(ptsave,'Averages')
        if not os.path.isdir(ptsave_b):
            os.mkdir(ptsave_b)

       
        var =self.LGV
        index = self.Label 
        zz,xx = np.meshgrid(z,x_s)
        it = 0 
        for values in var: 

       
            fg = figure()
    
            tick=r'%s Time = %s Myrs, $t^{\dagger}$ = %s' %(index[it],time_sim,time_dimen)

            ptsave_c=os.path.join(ptsave_b,values)
            if not os.path.isdir(ptsave_c):
                os.mkdir(ptsave_c)
        
            fn = os.path.join(ptsave_c,fna)

            ax0 = fg.gca()
            if values == 'D':
                cor = 1/(IC.D0[0][0]/1e3)
            elif values =='tau':
                cor = 1/(IC.tau0[0][0]/1e6)
            elif values == 'eps':
                cor = 1/IC.epsc[0][0]
            else:
                cor = 1.0
            buf = eval(values,globals(),self.__dict__)[:,:,ipic]*cor
            buf[buf == -np.inf]=np.nan
            lm    = self.Val[it]
            if lm[0]=="min":
                lm1=np.nanmin(buf)
            else:
                lm1 = lm[0]
            if lm[1]=="max":
                lm2=np.nanmax(buf)
            else:
                lm2 = lm[1]
            plt.grid(True)

            if((values == "eps") | (values == "Psi")):
                cf=ax0.pcolormesh(xx,zz,np.log10(buf),cmap='inferno',vmin=np.log10(lm1),vmax=np.log10(lm2))
                cbar = fg.colorbar(cf,ax=ax0,orientation='horizontal')
            else:
                cf=ax0.pcolormesh(xx,zz,buf,cmap='inferno',vmin=lm1,vmax=lm2)
                cbar = fg.colorbar(cf,ax=ax0,orientation='horizontal')

            

            ax0.set_title(tick)
            ax0.tick_params(axis='both', which='major', labelsize=5)
            ax0.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
            plt.draw()    # necessary to render figure before saving
            fg.savefig(fn,dpi=300)
            ax0.plot()
            it += 1
            val = [] 
            
            plt.close()
        
    def  _plot_slab_time(self,time,z,ptsave,TestName,Data_BASE_F,t_lim,IC): 

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
        self.tau_d_min = self.tau_min[ind_z,:]
        self.tau_d_max = self.tau_max[ind_z,:]
        if len(self.det_vec) >0:
            self._save_txt_DATABASE(self.det_vec,TestName,Data_BASE_F)
        

        for values in var: 
            print(values)
            fna='Fig%s.png' %values
            if(values == "dF"):
                val = self.F_T/self.F_B
            else:
                val = eval(values,globals(),self.__dict__)

            
            if values == "D":
                val = val/80 
    #        if values == "T":
     #           val = (val-IC.T_av)/IC.Tc

            
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
            lim_m = 0.1
            lim_M = 1.0
            print(lim_m)
            print(lim_M)
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
        
class Initial_condition(): 
    def __init__(self,Phase_Slab,Phase_Mantle,vIC): 
        """
        Parameters
        ----------
        path :[string]  
            I assume that the path for each category of tests is the same
        TestName : [string] 
            Testname just for addressing to the right test 

        Fill up the initial condition in structures, as such in a self-consistent manner
        the rheology, density, and so forth are filled up for then being saved in the 
        data base. 
        -------
        ?) reading geometry as well? 
        ?) reading the python file to generate the initial setup (later on)
        """
        # Data that are input of the script 
        self.D0,self.L0 = vIC.Slab.D0*1e3,vIC.Slab.L0*1e3
        self.T_av       = vIC.Slab.avTS+273.15
        self.TP         = vIC.Slab.TP+273.15
        rho_slab         = Phase_Slab.Density.rho*(1-Phase_Slab.Density.alpha*(self.T_av-273.15))
        rho_mantle      = Phase_Mantle.Density.rho*(1-Phase_Slab.Density.alpha*(self.TP-273.15))
        # Reference Buoyancy force 
        self.F_B0       = 9.81*(rho_slab-rho_mantle)*self.L0*self.D0
        print("Initial Bouancy force F_B is %e"%(self.F_B0))
        ###############################################
        # Compute the characteristic time following Thielmann 
        ###############################################
        # P_ref assuming a thickness of 100 km with a material of average density of 
        # 2900 
        Pr  = 2900*100e3*9.81 
        tau0  = self.F_B0/(2*self.D0)
        Bn    = Phase_Slab.Rheology.Dislocation.B
        n     = Phase_Slab.Rheology.Dislocation.n 
        Bd    = Phase_Slab.Rheology.Diffusion.B
        Vd    = Phase_Slab.Rheology.Diffusion.V
        Vn    = Phase_Slab.Rheology.Dislocation.V
        En    = Phase_Slab.Rheology.Dislocation.E
        Ed    = Phase_Slab.Rheology.Diffusion.E
        R     = 8.31446261815324
        expD = np.exp(-(Ed+Pr*Vd)/(R*self.T_av))
        expN = np.exp(-(En+Pr*Vn)/(R*self.T_av))
        eta0D = 0.5*(1/Bd)*tau0**(1-1)*1/expD
        eta0N = 0.5*(1/Bn)*tau0**(1-n)*1/expN
        # Compute the average viscosity of the upper mantle along the slab using the analytical formula that I derived
        # The analytical solution is not accounting for the adiabatic heating. 	
        expDM0 = np.exp(-(Ed+Pr*Vd)/(R*self.TP)) # Exponential term for the um with reference condition
        expNM0 = np.exp(-(En+Pr*Vn)/(R*self.TP)) # 
        eta0DM0 = 0.5*(1/Bd)*1/expDM0 # Reference viscosity diffusion
        eta0NM0 = 0.5*(1/Bn)*tau0**(1-n)*1/expNM0 # Reference viscosity dislocation
        Cd = (Vd)/(R*self.TP) # Exponential correction 
        Cn = (Vn)/(R*self.TP)
        w_m = rho_mantle*9.81 # weight force
        expNM = np.exp(Cn*w_m*self.L0) # exponential term @ bottom of the slab
        expDM = np.exp(Cd*w_m*self.L0)
        eta0DM = (eta0DM0/(Cd*w_m*self.L0))*(expDM-1) # Average dislocation/diffusion creep viscosity
        eta0NM = (eta0NM0/(Cn*w_m*self.L0))*(expNM-1)
        eta_ref_UM = (1/eta0DM+1/eta0NM)**(-1) #Reference effective visocosity of the um 
        print("The average effective viscosity of the mantle at reference condition is %3f"%(np.log10((1/eta0DM+1/eta0NM)**(-1))))
        eta_S_ref = (1/eta0D+1/eta0N)**(-1)
        eta_S_co  = (1/eta0D+1/eta0N+1/1e24)**(-1)
        print("The effective viscosity at the reference condition of the slab is %3f and with cut off is %3f"%(np.log10(eta_S_ref), np.log10(eta_S_co)))
        print("The initial phi without cut off is %2f and with cutoff is %2f and if the viscosity of the mantle is low than cutoff %2f" %(np.log10(eta_ref_UM/eta_S_ref), np.log10(eta_ref_UM/eta_S_co), np.log10(1e18/eta_S_co)))
        Bd_cuoff = 1/(2*1e24) # Bd upper cut off
        tc    = (Bn*(tau0)**n*expN + Bd*tau0*expD+Bd_cuoff*tau0) # epsc of the slab
        print("Strain rate dislocation creep is %3f"%(np.log10(Bn*tau0**n*expN)))
        print("Strain rate diffusion creep is %3f"%(np.log10(Bd*expD)))
        print("Strain rate upper cut off creep is %3f"%(np.log10(Bd_cuoff*tau0)))
        print("Characteristic strain rate is %3f"%(np.log10(tc)))
        self.epsc = tc 
        tc = tc**(-1)
        tc    /= (365.25*60*60*24*1e6) 
        print("The characteristic time of the simulation is %1f" %tc)
        print("The analytical time  of detachment is %1f" %(tc/n))
        self.tc = tc        #Characteristic deformation time 
        self.td = tc/n      #Detachment time
        self.tau0 = tau0    #tau0 
        self.eta0DS = eta0D #Diffusion creep viscosity slab @ reference
        self.eta0DN = eta0N #Dislocation creep viscosity slab @ reference
        self.xiUS   = eta0D/eta0N # Characteristic viscosity contrast @ reference condition
        self.xiUM   = eta0DM/eta0NM # Characteristic viscosity contrast @ reference condition UM
        self.Psi_R  = eta_ref_UM/eta_S_ref # Real Psi (i.e., the Psi computed with the real data)
        self.Psi_R_Sco = eta_ref_UM/eta_S_co #  Psi accounting for the upper cut off
        self.Psi_co = 1e18/eta_S_co # Psi in case the effective viscosity of the um is below the cut off viscosity
        self.eta_ref_UM = eta_ref_UM # the effective viscosity of the upper mantle 
        self.Tc     = self.TP-self.T_av # Characteristic temperature
        self.T_av   = self.T_av-273.15
        self.TP         = self.TP-273.15
        self.eta_ref_S = eta_S_co 

        print("Reference viscosity diffusion creep of the slab is %f"%(np.log10(eta0D)))
        print("Reference viscosity dislocation creep of the slab is %f"%(np.log10(eta0N)))
        print("XiS %f"%(np.log10(self.xiUS)))


def _plot_D_D0(Slab,IC,ptsave,time,Test_Name,t_lim): 
    # 
    fg = figure()

    tick='%s' %(Test_Name)

    fna = "D_vs_Analytical.png"
    tt = np.linspace(0,1.0,num=10000)
    analytical_solution = (1-3.5*tt)**(1/3.5)

    ptsave_c=os.path.join(ptsave)
    if not os.path.isdir(ptsave_c):
        os.mkdir(ptsave_c)

    fn = os.path.join(ptsave_c,fna)

    ax0 = fg.gca()

    ax0.plot(time/IC.td,Slab.D_det/(IC.D0/1e3), lw=1.2,c='k',ls='--', label = r"Simulations")
    ax0.plot(tt/(1/3.5),analytical_solution, lw=1.5,c='k',ls='dashdot', label = r"Analytical Solution")
    ax0.scatter(time/IC.td,Slab.D_det/(IC.D0/1e3),c='r',s=10)
    plt.xlim(0, t_lim)
    plt.ylim(0.1, 1.0)
    ax2 = ax0.twinx()  
    ax2.set_ylim(0.5,5.0)
    ax2.plot(time/IC.td,Slab.tau_vec,lw=1.2,c='b',label=r"$\tau^{\dagger}_{SIM}$")
    ax2.plot(time/IC.td,Slab.tau_d_min,lw=0.8,c='b',label=r"$\tau^{\dagger}_{SIM}$")
    ax2.plot(time/IC.td,Slab.tau_d_max,lw=0.8,c='b',label=r"$\tau^{\dagger}_{SIM}$")
    ax2.fill_between(time/IC.td, Slab.tau_d_min, Slab.tau_d_max,alpha=0.2)
    ax2.plot(tt/(1/3.5),1/analytical_solution,lw=0.8,c='b',ls='dashdot')
    ax2.yaxis.set_ticks(np.arange(0.5,5.0,0.5))

    ax0.grid(True)
    #plt.yscale('log')
    ax0.set_title(tick)
    plt.xlabel(r'$t^{\dagger}, [n.d.]$')
    plt.ylabel(r'$D^{\dagger}, [n.d.]$')
    ax0.tick_params(axis='both', which='major', labelsize=5)
    ax0.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
    ax2.tick_params(axis='both', which='major', labelsize=5)
    ax2.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
   
    plt.draw()    # necessary to render figure before saving
    fg.savefig(fn,dpi=300,transparent=False)
    ax0.plot()
    val = [] 
    
    plt.close()
    


def _plot_F_T_F_B(Slab,IC,ptsave,time,Test_Name,t_lim,ind_z): 
    
    fg = figure()

    tick='%s' %(Test_Name)

    fna = "Buoyancy_FORCE.png"
    tt = np.linspace(0,1.0,num=10000)
    analytical_solution = (1-3.5*tt)**(1/3.5)

    ptsave_c=os.path.join(ptsave)
    if not os.path.isdir(ptsave_c):
        os.mkdir(ptsave_c)

    fn = os.path.join(ptsave_c,fna)

    ax0 = fg.gca()

    ax0.plot(time/IC.td,Slab.F_T[ind_z,:]/IC.F_B0, lw=1.2,c='k',ls='--', label = r"Simulations")
    
    
    leg = plt.legend()
    plt.xlim(0, t_lim)
    plt.ylim(0, 1.1)
    plt.grid(True)
    #plt.yscale('log')
    ax0.set_title(tick)
    plt.xlabel(r'$\frac{t}{t_d}, [n.d.]$')
    plt.ylabel(r'$\frac{F_T}{F_B0}, [n.d.]$')
    ax0.tick_params(axis='both', which='major', labelsize=5)
    ax0.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
    plt.draw()    # necessary to render figure before saving
    fg.savefig(fn,dpi=300,transparent=False)
    ax0.plot()
    val = [] 
    
    plt.close()
    

    

def _plot_vz(Slab,IC,ptsave,time,Test_Name,t_lim,ind_z): 
    
    fg = figure()

    tick='%s' %(Test_Name)

    fna = "velocity_slab_detached.png"
    tt = np.linspace(0,1.0,num=10000)
    analytical_solution = (1-3.5*tt)**(1/3.5)

    ptsave_c=os.path.join(ptsave)
    if not os.path.isdir(ptsave_c):
        os.mkdir(ptsave_c)

    fn = os.path.join(ptsave_c,fna)

    ax0 = fg.gca()

    ax0.plot(time/IC.td,Slab.vz_S[:,0], lw=1.2,c='b',ls='--', label = r"$v_z mean [cm/yrs]$")
    ax0.plot(time/IC.td,Slab.vz_S[:,1], lw=1.2,c='b',ls=':', label = r"$v_z median [cm/yrs]$")
    ax0.plot(time/IC.td,Slab.vz_S[:,3], lw=0.5,c='b',ls=':')
    ax0.plot(time/IC.td,Slab.vz_S[:,4], lw=0.5,c='b',ls=':')
    plt.axvline(x=Slab.SDet.t_det_td,c = 'r')
    
    ax0.fill_between(time/IC.td,Slab.vz_S[:,3],Slab.vz_S[:,4],color="b",alpha=0.1)

    plt.grid(True)
    
    leg = plt.legend()
    plt.xlim(0, t_lim)
    plt.grid(True)
    #plt.yscale('log')
    ax0.set_title(tick)
    plt.xlabel(r'$\frac{t}{t_d}, [n.d.]$')
    plt.ylabel(r'$v_z, [cm/yrs]$')
    ax0.tick_params(axis='both', which='major', labelsize=5)
    ax0.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
    plt.draw()    # necessary to render figure before saving
    fg.savefig(fn,dpi=300,transparent=False)
    ax0.plot()
    val = [] 
    
    plt.close()
    

def _plot_time_map_surface(x,time,buf,Field,Test_Name,cmap,ptsave,clim,t_lim):
    import cmcrameri.cm as cmc
    cmap2 =eval(cmap)
    ptsave_b=os.path.join(ptsave,'Maps')
    if not os.path.isdir(ptsave_b):
        os.mkdir(ptsave_b)
    
    fg = figure()
    ax0 = fg.gca()
    
    fna=('%s_%s.png' %(Test_Name,Field))
    cbar_title = ('%s, mm/yrs' %(Field))
    fn = os.path.join(ptsave_b,fna)
    cf=ax0.pcolormesh(x, time, np.transpose(buf),cmap = cmap2, shading='gouraud', vmin=clim[0], vmax=clim[1])
    cbar = fg.colorbar(cf, ax=ax0,orientation='horizontal')
    ax0.tick_params(axis='both', which='major', labelsize=5)
    ax0.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
    ax0.set_ylim(np.min(time),t_lim)
    ax0.set_xlim(np.min(x),np.max(x))
    plt.xlabel('$x, [km]$')
    plt.ylabel('$Time, [Myrs]$')
    cbar.set_label(cbar_title)
 
    plt.draw()    # necessary to render figure before saving
    fg.savefig(fn,dpi=300)
    fg.clear()
    plt.close()

class Initial_Geometry():
    def __init__(self, file_mat):
        self.Continent1 = Terrane_Geo(file_mat,'Continent1')
        self.Continent2 = Terrane_Geo(file_mat,'Continent2')
        self.Slab       = Trench(file_mat,'Trench')
        self.Ocean_BG   = Terrane_Geo(file_mat,'Ocean_BG')
        
                 

class Terrane_Geo():
      def __init__(self, file_mat,name):
          mat = h5py.File(file_mat)
          TB = mat['TB']
          buf = TB[name]
          TI = buf['Thermal_type']
          # Thermal Information 
          self.Age = np.array(TI['Age'])
          self.Moho_T = np.array(TI['Moho'])
          self.Moho_z = np.array(TI['Moho_d'])
          string_type = []
          aa = np.array(TI['Type'])
          self.Type   = ''.join((str(np.string_(i))[2]) for i in aa[:])
          try:
              self.vel    = np.concatenate(np.array(mat[TI['vel'][0][0]]))
          except:
              self.vel = np.array(TI['vel'])
          # Boundary information 
          Bounds = buf['Boundary']
          self.B_main_coordinate = np.array([np.array(Bounds['x1'])[0][0],np.array(Bounds['x2'])[0][0],np.array(Bounds['y1'])[0][0],np.array(Bounds['y2'])[0][0]])
          self.Boundary_A = self._extract_boundaries(mat,Bounds,'A')
          self.Boundary_B = self._extract_boundaries(mat,Bounds,'B')
          self.Boundary_C = self._extract_boundaries(mat,Bounds,'C')
          self.Boundary_D = self._extract_boundaries(mat,Bounds,'D')
          try:
              self.Phases     = np.concatenate(np.array(buf['Stratigraphy/phases']))
              self.Thickness  = np.concatenate(np.array(buf['Stratigraphy/Tk']))
          except: 
              self.Phases     = np.concatenate(np.array(buf['Stratigraphy_Oceanic/phases']))
              self.Thickness  = np.concatenate(np.array(buf['Stratigraphy_Oceanic/Tk']))
          
      def _extract_boundaries(self,mat,Bounds,Boundary):
         B_A = []
         A=Bounds[Boundary]
         ind = np.array(A)[0][0]
         B_A.append(np.array([np.array(mat[ind])[0][0], np.array(mat[ind])[1][0],np.array(mat[ind])[2][0],np.array(mat[ind])[3][0]]))
         ind =  np.array(A)[1][0]
         B_A.append(''.join((str(np.string_(i))[2]) for i in np.array(mat[ind])[:]))
         if B_A[1] == 'none':
             B_A.append([])
         else:
             ind = np.array(A)[2][0]
             B_A.append(np.array([np.array(mat[ind])[0][0], np.array(mat[ind])[1][0],np.array(mat[ind])[2][0]]))
            
         return B_A
    
class Trench(Terrane_Geo):
    def __init__(self, file_mat,name):
        super().__init__(file_mat, name)
        mat = h5py.File(file_mat)
        TB = mat['TB']
        buf = TB[name]
        self.D0 = np.array(buf['D0'])
        self.L0 = np.array(buf['L0'])
        self.R  = np.concatenate(buf['R'])
        self.BList = str(np.string_(np.array(mat[buf['Boundaries_list'][0][0]])))[2]
        self.theta = np.concatenate(np.array(mat[buf['theta'][0][0]]))
        self.continent_S = np.concatenate(np.array(mat[buf['length_continent'][0][0]]))
        self.avTS = np.array(buf['TS'])
        self.TP   = np.array(buf['Thermal_information/TP'])
        
class Free_S_Slab_break_off(FS):
    def __init__(self,C,nstep):
        super().__init__(C,nstep)
        tx  = np.sum(C.ind_x)
        ty  = np.sum(C.ind_y)
        self.dH          = np.zeros((ty,tx,nstep),dtype=float)
        self.vz_M        = np.zeros((ty,tx,nstep),dtype=float)
        self.mean_stress = np.zeros((ty,tx,nstep),dtype=float)
        self.mean_eps = np.zeros((ty,tx,nstep),dtype=float)
        self.mean_dx     = np.zeros((ty,tx,nstep),dtype=float)
        self.mean_dy     = np.zeros((ty,tx,nstep),dtype=float)
        self.mean_dz     = np.zeros((ty,tx,nstep),dtype=float)
        self.HB           = np.ones((tx,nstep),dtype=float)*np.nan # topography boundary plates
        self.HBy          = np.ones((tx,nstep),dtype=float)*np.nan # coordinate boundary plates
        self.HBx          = np.ones((tx,nstep),dtype=float)*np.nan # coordinate boundary plates
        self.HMax         =np.ones((tx,nstep),dtype=float)*np.nan # maximum topography within +/-200 km boundary
        self.HMaxy        = np.ones((tx,nstep),dtype=float)*np.nan # coordinate maximum topography
        self.Hmin       = np.ones((tx,nstep),dtype=float)*np.nan # minimum frontal topography
        self.Hminy      = np.ones((tx,nstep),dtype=float)*np.nan # coordinate
        self.v_z_M       = np.ones((tx,nstep),dtype=float)*np.nan # maximum velocity
        self.v_z_m       = np.ones((tx,nstep),dtype=float)*np.nan #minimum velocity
        self.v_z_mean    = np.ones((tx,nstep),dtype=float)*np.nan #minimum velocity
        self.HMean       = np.ones((tx,nstep),dtype=float)*np.nan #minimum velocity
        self.Hmeang       = np.ones((nstep),dtype=float)*np.nan #minimum velocity

        self.x_s         = np.ones(tx)*np.nan
        self.ind_boundary = np.zeros(tx)*np.nan
        
        
        
    def _update_extra_variables(self,V:VAL,C:Coordinate_System,dt,ipic):
        self.dH[:,:,ipic]          = self._update_FS(dt,ipic,'dH')
        self.vz_M[:,:,ipic]        = self._update_FS(dt,ipic,'vz_M')
        self.mean_stress[:,:,ipic] = self.update_Mean_CC(V,C,ipic,'tau')
        self.mean_eps[:,:,ipic]    = self.update_Mean_CC(V,C,ipic,'eps')

        self.mean_dx[:,:,ipic]     = self.update_Mean_CC(V,C,ipic,'dx')
        self.mean_dy[:,:,ipic]     = self.update_Mean_CC(V,C,ipic,'dy')
        self.mean_dz[:,:,ipic]     = self.update_Mean_CC(V,C,ipic,'dz')
        
        
        
        
        self.LGV = ['dH','vz_M','mean_stress','mean_eps','mean_dz','Amplitude']
        self.Label  = ['dH','vz','tau','eps','dz','H']
        self.Colormap = ["cmc.cork","cmc.cork","cmc.bilbao","cmc.devon","cmc.cork","cmc.oleron"]
        self.Val         = [(-3,3),
                            (-3,3),
                            (0,"max"),
                            (1e-16,"max"),
                            ("min","max"),
                            (-4.0,3.0)]
        
        return self
        
        
    def _update_FS(self,dt,ipic,name_field):
        if name_field == 'dH':
            if ipic > 0:
                buf = (self.Topo[:,:,ipic]-self.Topo[:,:,ipic-1])/dt
            else:
                buf = self.dH[:,:,ipic]
        else:
            if ipic > 0:
                buf = (self.vz[:,:,ipic]+self.vz[:,:,ipic-1])/2
            else: 
                buf = self.vz[:,:,ipic]
        return buf 
    
    def update_Mean_CC(self,V:VAL,C:Coordinate_System,ipic,field_name):
        # Compute the continental crust
        Continental_crust                           = V.Sed+V.CC1+V.CC2 
        Continental_crust[Continental_crust>=0.8]      = 1.0
        Continental_crust[Continental_crust<0.8]    = np.nan
        buf = 0.0*Continental_crust
        buf = buf[0,:,:]
        # loop over x,y. Select colum: compute mean of those that has high crustal fraction 
        for i in range(len(C.y)):
            for j in range(len(C.x)): 
                column = eval(field_name,globals(),V.__dict__)[:,i,j]
                crust  = Continental_crust[:,i,j]
                if len(crust[crust==1]) == 0:
                    buf[i,j]=np.nan 
                else:
                    buf[i,j] = np.mean(column[crust==1])
        return buf 
    def _plot_maps_FS(self,t_cur,y,x,ptsave,ipic):
        import cmcrameri.cm as cmc
        from matplotlib.colors import LogNorm
        ptsave_b=os.path.join(ptsave,"WhMaps_FS")
        
        if not os.path.isdir(ptsave_b):
            os.mkdir(ptsave_b)
        
        values = self.LGV 
        index = self.Label
        cmaps = self.Colormap 
        LIM   = self.Val
        
        
        time_sim = "{:.3f}".format(t_cur)
      
        ic = 0  
        val = np.zeros((len(y),len(x)),dtype=float)
        fg = figure()
        ax0 = fg.gca()
        fna='Fig'+"{:03d}".format(ipic)+'.png'       
        cf =ax0.pcolormesh(x, y, val, shading='gouraud')
        cbar = fg.colorbar(cf, ax=ax0,orientation='horizontal',extend="both")
        cf0 = ax0.contour(x,y,self.Amplitude[:,:,ipic],levels = [-1000,-500,0,500,1000],colors = "k",linewidths=0.75)
        
       
        for name in values:
        
            cmap2 = eval(cmaps[ic])

            val = eval(name,globals(),self.__dict__)
            val = val[:,:,ipic]
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
            ax0.tick_params(axis='both', which='major', labelsize=5)
            ax0.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
            ax0.set_ylim(np.min(y),np.max(y))
            ax0.set_xlim(np.min(x),np.max(x))
            cbar.set_label(index[ic])
            ax0.set_title(tick)

            #plt.draw()    # necessary to render figure before saving
            
            fg.savefig(fn,dpi=300,transparent=False)
            
            
            val = [] 

            ic +=1 
        fg.clear
        plt.close()

    def _compute_relevant_information_topography(self,C: Coordinate_System,Slab_GEO:Trench,ipic:int):
        # Compute the plate boundary using the data of the boundary 
            """
            Approach: compute x-y of the boundary of the terrane
            -> interpolate topography at the boundary of the terrane and other variables
            -> loop over the nodes of the point of the margin:
                -> select the area between p(B)-200,p(B)+200 
                    -> find the maximum topography (store its node)
                    -> find the minimum topography:
                        -> behind and in front the maximum topography 
                    -> compute the averages, the wavelength (i.e. distance between the two minima)
            """
            # Extract information boundary and select x belongs to xa-xb 
            boundary_geometry = Slab_GEO.Boundary_B
            xa = boundary_geometry[0][0]
            xb = boundary_geometry[0][2]
            ind_boundary = (C.x>=xa) & (C.x<xb)
            x_B = C.x[ind_boundary]
            # Compute the boundary 
            y_B,x_s =_compute_length_coordinatesCB(ind_boundary,boundary_geometry,C.x)
            # function to retrieve the data. 
            vz = self.vz_M[:,:,ipic]
            H  = self.Amplitude[:,:,ipic]
            self.HBx[ind_boundary==True,ipic]=x_B[:]
            self.HBy[ind_boundary,ipic]=y_B[:]
            self.x_s = x_s

            
            """
            HB     topography of plate boundary  
            HBc    coordinate of the plate boundary
            HMax   topography of the maximum   
            HMaxc  coordinate maximum
            HminF  topography frontal minumum
            HminFc topography frontal minimum coordinates
            HminB  topography back minimum
            Hminc  coordinate topography back mininum
            v_z_M  max velocity     [HBc+/-200]
            v_z_m  minumum velocity [HBc+/-200]
            x_s    length from left to right 
            """
            ix_b = np.sum(ind_boundary)
            ix_ch = np.where(ind_boundary==True)
            for i_x in range(ix_b):
                i = ix_ch[0][i_x]
                yy = y_B[i_x] # 
                iy = find1Dnodes(C.y,yy,len(C.y))
                y1 = C.y[iy]
                y2 = C.y[iy+1]
                intp1 = H[iy,i]
                intp2 = H[iy,i]
                self.HB[i,ipic] = linearinterpolation(yy,y1,y2,intp1,intp2)
                # find area of interest: 
                ind_area = np.where((C.y<=self.HBy[i,ipic]+200)&(C.y>=self.HBy[i,ipic]-200))
                # find_maximum
                self.HMax[i,ipic] = np.max(H[ind_area,i])
                ind_max   = np.where(H[:,i]==self.HMax[i,ipic])
                # find_coordinate
                if C.y[ind_max[0][0]]<self.HBy[i,ipic]-200:
                    self.HMaxy[i,ipic]=np.nan
                else:
                    self.HMaxy[i,ipic]=C.y[ind_max[0][0]]
                self.Hmeang[ipic] = np.mean(H[ind_area[0][0],ind_boundary==True])
                # find_minimum_front
                if self.HMaxy[:,ipic].all() == False:
                    self.Hmin[i,ipic]  = np.min(H[ind_area,i])
                    self.HMean[i,ipic]  = np.mean(H[ind_area,i])  
                    self.Hminy[i,ipic] = C.y[H[:,i]==np.min(H[ind_area,i])]
                    self.v_z_M[i,ipic] = np.max(vz[ind_area,i])
                    self.v_z_m[i,ipic] = np.min(vz[ind_area,i])
                    self.v_z_mean[i,ipic] = np.mean(vz[ind_area,i])

                else:
                    self.Hmin[i,ipic]   = np.nan
                    self.Hminy[i,ipic]  = np.nan
                    self.HMean[i,ipic]  = np.mean(H[ind_area,i])  
                    self.v_z_M[i,ipic]   = np.nan
                    self.v_z_m[i,ipic]   = np.nan
                    self.v_z_mean[i,ipic] = np.nan
            self.ind_boundary=ind_boundary
            return self 
    def _plot_1D_plots_Free_surface(self,ipic: int,ptsave):
        val = ['HMax','v_z','Maps']
        ptsave=os.path.join(ptsave,'1D_surfaces')
        if not os.path.isdir(ptsave):
            os.mkdir(ptsave)
        ib = self.ind_boundary
        for v in val:
            ptsave_b=os.path.join(ptsave,v)
            if not os.path.isdir(ptsave_b):
                os.mkdir(ptsave_b)
            fig = plt.figure()
            ax2 = fig.gca()
            fna='Fig'+str(ipic)+'.png'
            fn=os.path.join(ptsave_b,fna)
            if v=='HB':
                for ip in range(10):
                    it = ipic - ip
                    alpha_v= 0.8-ip*(1/12)
                    if ip == 0: 
                        cc = 'r'
                    else:
                        cc = 'b'
                    if (it == 0) & (ip == 0) :
                        ax2.plot(self.x_s, self.HB[ib==True,0]-self.Hmean[0],c = cc,alpha = alpha_v,linewidth=alpha_v)
                        break
                    if (ip >0 ) & (it == 0 ):
                        ax2.plot(self.x_s, self.HB[ib==True,it]-self.Hmean[it],c = cc,alpha = alpha_v,linewidth=alpha_v)
                        break 
                    else: 
                        ax2.plot(self.x_s, self.HB[ib==True,it]-self.Hmean[it],c = cc,alpha = alpha_v,linewidth=alpha_v)
                ax2.plot(self.x_s, self.HB[ib==True,ipic]-self.Hmean[ipic],c = 'r',alpha = 1.0,linewidth=1.2)
                plt.xlabel('x, [km]')
                plt.ylabel('H, [km]')
                plt.xlim(0,1200)

            elif v == 'v_z':
                p1=ax2.plot(self.x_s,self.v_z_M[ib==True,ipic],c = 'b',alpha = 1.0,linewidth=0.8)
                p2=ax2.plot(self.x_s,self.v_z_m[ib==True,ipic],c = 'b',alpha = 1.0,linewidth=0.8)
                p3 =ax2.fill_between(self.x_s, self.v_z_m[ib==True,ipic], self.v_z_m[ib==True,ipic],color='blue',alpha=0.4)
                p4=ax2.plot(self.x_s,self.v_z_mean[ib==True,ipic],c = 'b',alpha = 1.0,linewidth=0.8)
                plt.xlabel('x, [km]')
                plt.ylabel('$v_z$, $[\frac{cm}{yrs}]$')
                plt.xlim(0,1200)

            elif v == 'HMax':
                p1=ax2.plot(self.x_s,self.HMax[ib==True,ipic]-self.Hmeang[ipic],c = 'r',alpha = 1.0,linewidth=1.0,linestyle=':')
                p2=ax2.plot(self.x_s,self.Hmin[ib==True,ipic]-self.Hmeang[ipic],c = 'b',alpha = 1.0,linewidth=1.2)
                plt.xlabel('x_s, [km]')
                plt.ylabel('H, [km]')
                plt.xlim(0,1200)
            elif v == 'Maps':
                p1=ax2.fill_between(self.HBx[ib==True,ipic], self.Hminy[ib==True,ipic],self.HMaxy[ib==True,ipic],color='blue',alpha=0.4)
                p3=ax2.plot(self.HBx[ib==True,ipic],self.HBy[ib==True,ipic],c = 'k',alpha = 1.0,linewidth=0.8,linestyle=':')
                plt.xlabel('x, [km]')
                plt.ylabel('y, [km]')
                
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

class Phase_det(Phase):
    def __init__(self,C,Phase_dic):
        super().__init__(C,Phase_dic)
        tx    = len(C.xp)
        ty    = len(C.yp)
        tz    = len(C.zp)
        self.tau   = np.zeros([tz,ty,tx],dtype = float)
        self.eps   = np.zeros([tz,ty,tx],dtype = float)
        self.T     = np.zeros([tz,ty,tx],dtype = float)
        self.nu    = np.zeros([tz,ty,tx],dtype = float)
        self.vis   = np.zeros([tz,ty,tx],dtype = float)
        self.Rho   = np.zeros([tz,ty,tx],dtype = float)
        self.Psi   = np.zeros([tz,ty,tx],dtype = float)
    def _interpolate_dyn_phase(self,V,C):
        # prepare the variables
        xp = C.xp
        yp = C.yp
        zp = C.zp
        x  = C.x
        y  = C.y
        z  = C.z 
        t1= perf_counter()
        val_ = ['T','nu','vis','Rho','tau','eps','Psi']
        for v in val_: 
            buf = eval(v,globals(),V.__dict__)
            buf2 = np.zeros([len(zp),len(yp),len(xp)],dtype = float)
            buf2 = function_interpolate(xp,yp,zp,x,y,z,buf,buf2)
            eval(v,globals(),self.__dict__)[:,:,:]=buf2 
            
        t2= perf_counter()
        print('time interpolation',"{:02}".format(t2-t1), 's')
        return self 
        
@jit(nopython=True,parallel=True)
def function_interpolate(xp,yp,zp,x,y,z,buf,buf2):
    for k in prange(len(zp)):
            for j in prange(len(yp)):
                for i in prange(len(xp)):
                    xx = xp[i]
                    yy = yp[j]
                    zz = zp[k]
                    ix = find1Dnodes(x,xx,len(x))
                    iy = find1Dnodes(y,yy,len(y))
                    iz = find1Dnodes(z,zz,len(z))
                    x1 = x[ix]
                    x2 = x[ix+1]
                    y1 = y[iy]
                    y2 = y[iy+1]
                    z1 = z[iz]
                    z2 = z[iz+1]
                    intp1,intp2,intp3,intp4,intp5,intp6,intp7,intp8 = findnodes(buf,ix,iy,iz)
                    val_ = trilinearinterpolation(xx,yy,zz,x1,x2,y1,y2,z1,z2,intp1,intp2,intp3,intp4,intp5,intp6,intp7,intp8)
                    buf2[k,j,i]=val_
    return buf2

@jit(nopython=True,parallel=True)
def _Find_Slab_PERFORM_C(x,z,buf_var_ph,ph,tz,ix,D,buf_var,x1,x2,z_bottom):    #(self,C,Values,ph,ipic,ix,ix1):

    for i in prange(tz):
        buf = ph[i,:]
        buf2 = buf_var_ph[i,:]
        condition = 0.0
        z_bottom[i]=1e6
        if(z[i]>-80):
            condition = 1.0
        else:
            """
            This is simply an highlight: i´m idiot. i failed to see a fucking error in the index. 8 hours wasted, and I should be punished for my cursed mind
            """
            counter=0
            for i_a in range(len(buf)):
                if buf[i_a]==-1:
                    counter +=1
            
            if counter == len(buf):
                condition = 2.0        
        if(condition>0.0):
            x2[i]=np.nan
            x1[i]=np.nan
            D[i] = np.nan                
            buf_var[i] = np.nan
            if (condition==1.0):
                x2[i]=-np.inf
                x1[i]=-np.inf
                buf_var[i] = -np.inf
                D[i] = -np.inf
            
        else:
            z_bottom[i]=z[i]
            i1 = 0
            i2 = -np.inf
            i1,i2 = _find_index(buf,i1,i2)
            x1[i]= x[i1]
            x2[i] = x[np.int(i2)]
            D[i]= x2[i]-x1[i] 
            if(x2[i]-x1[i] <8.0):
                buf_var[i] = np.nan
                D[i]=np.nan               
            else:
                buf_var[i] = _mean_numba(buf,buf2)
    L0  =z-np.min(z_bottom)
    return D,L0,buf_var,x1,x2
    
@jit(nopython=True)  
def _mean_numba(buf,buf2):
    mean = 0.0
    len_ = np.int(0)
    for i in range(len(buf)):
        if buf[i]>0.95:
            mean+= buf2[i]
            len_ += np.int(1)
    mean = mean/len_
    
    return mean

@jit(nopython=True)  
def _find_index(buf,i1,i2):
    """_summary_

    Args:
        buf (float64): buffer of the phase layout
        i1 (_type_): index 1 
        i2 (_type_): index2

    Returns:
        _type_: index 1 and index2 update. 
        Short disclaimer: I was using break for interrupting the 
        search for the i2. However, it appears that time to time there might 
        be [not valid, valid, not valid]. One solution is to do a moving average 
        for smothening the phase boundaries. Or, to keep looking. This might create 
        artifacts. So, I will try to introduce the moving average approach
        1. Discrete field -> [not valid, valid, not valid, valid:10,not valid] (i2) is 1 
    """
    # Thanks Arne for the tip. 
    for ib in range(len(buf)):
        if (buf[ib]>0.95) & (i1 == 0):
            c1 = ib+1
            c2 = ib+2
            if c1 >len(buf)-1:
                c1 = len(buf)-1
                c2 = len(buf)-1
            if(buf[ib+1]==-1)|(buf[ib+2]==-1):
                i1 == 0 
            else:
                i1 = ib
                i2 = 0
                break 

    for ib in range(len(buf),-1,-1):
        if (buf[ib]>0.95) & (i2 == 0):
            c1 = ib-1
            c2 = ib-2
            if c1 <0:
                c1 = 0
                c2 = 0
            
            if(buf[c1]==-1)|(buf[c2]==-1):
                i2 == 0 
            else:
                i2 = ib
                break 
        
    return i1,i2

#@jit(nopython=True)  
def _compute_length_coordinatesCB(ind_boundary,boundary_geometry,x):
    x = x[ind_boundary]
    # Compute the arclength as a function of x
    x_a = boundary_geometry[0][0]
    y_a = boundary_geometry[0][1]
    c   = boundary_geometry[2][0]
    cy  = boundary_geometry[2][1]
    R   = boundary_geometry[2][2]
    c_  = (x_a-c)**2-R**2+y_a**2
    b_ = 2.*y_a
    delta = np.sqrt(b_**2-4*c_)
    center_y = [(b_-delta)/2,(b_+delta)/2]
    y_c = np.min(center_y)
    y = (R**2-(x-c)**2)**(0.5)+y_c
    d = np.sqrt((x-x_a)**2+(y-y_a)**2)
    theta = 2*np.arcsin(d/(2*R))
    x_s = R*theta; 
    return y,x_s 