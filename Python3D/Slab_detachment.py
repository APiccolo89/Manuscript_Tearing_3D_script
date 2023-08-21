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

        
        tz = np.sum(C.ind_zp)
        self.D     = np.ones((tx,tz,nstep),dtype=float)*(-1.0)
        self.T     = np.zeros((tx,tz,nstep),dtype=float)
        self.nu    = np.zeros((tx,tz,nstep),dtype=float)
        self.eps   = np.zeros((tx,tz,nstep),dtype=float)
        self.tau   = np.zeros((tx,tz,nstep),dtype=float)
        self.vis   = np.zeros((tx,tz,nstep),dtype=float)
        self.F_T   = np.zeros((tx,tz,nstep),dtype=float)
        self.F_B   = np.zeros((tx,tz,nstep),dtype=float)
        self.dRho  = np.zeros((tx,tz,nstep),dtype=float)
        self.vis_A    = np.zeros((tx,tz,nstep),dtype=float)
        self.nu_A    = np.zeros((tx,tz,nstep),dtype=float)
        self.det_vec =np.zeros((tx,nstep),dtype=float)
        self.x1    = np.zeros(tz)
        self.x2    = np.zeros(tz)   
        self.LGV         = ["D","T","tau","F_B","F_T"]
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
        self.Val         = [("0.1","1.0"),
                            ("400","1200"),
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
        self.CV = ["T","nu","vis","eps","tau"]       


    def _update_C(self,C,V,Ph,IG,ipic,time):
        lay_ph = [Ph.Phase][0]
        lay_ph = lay_ph.astype('float64')
        ip = len(IG.Slab.Phases)
        # For whatever reason the type has been change, and the 1000 is 232 in uint8
        for ic in range(ip):
            print(ic)
            print(IG.Slab.Phases[ic])
            lay_ph[lay_ph==IG.Slab.Phases[ic]]=(1000)
        lay_ph[lay_ph<(1000)] = 0.0
        lay_ph[lay_ph==(1000)] = 1.0 
        
        for i in range(len(self.ind_boundary[0])-1):
            ix = self.ind_boundary[0][i]
            self = self._Find_Slab_C(C,Ph,lay_ph,ipic,ix,i)
       
                
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
                        ind2 = np.where(buf==0.0)
                        Rho_AST = np.mean(Values.Rho[i,ind2,ix])
                        Rho_Slab = np.mean(Values.Rho[i,ind,ix])
                        self.dRho[ix1,i,ipic] = 0.0#Rho_AST-Rho_Slab
                        
            L_id = np.where(np.isnan(self.x1)==False) 
            z_b  = z[L_id[0][0]]
            dRho = -np.nanmean(self.dRho[ix1,:,ipic])
            self.F_T[ix1,:,ipic] = 2*self.tau[ix1,:,ipic]*self.D[ix1,:,ipic]*1e9
            self.F_B[ix1,:,ipic] = dRho*9.81*self.D[ix1,:,ipic]*(z-z_b)*1e6                
            return self 
    
    def  _plot_average_C(self,t_cur,x,z,ptsave,ipic): 

        x = x[self.ind_boundary]        
        time_sim = "{:.3f}".format(t_cur)
        
        ic  = 0
        
        fna='Fig'+str(ipic)+'.png'
        ptsave_b=os.path.join(ptsave,'Averages')
        if not os.path.isdir(ptsave_b):
            os.mkdir(ptsave_b)

       
        var =self.LGV
        index = self.Label 
        zz,xx = np.meshgrid(z,x)

        for values in var: 

       
            fg = figure()
    
            tick='%s Time = %s Myrs' %(index[ic],time_sim)

            ptsave_c=os.path.join(ptsave_b,values)
            if not os.path.isdir(ptsave_c):
                os.mkdir(ptsave_c)
        
            fn = os.path.join(ptsave_c,fna)

            ax0 = fg.gca()
            if values == 'D':
                cor = 1/80
                lim = (0.1,0.9)
            else:
                cor = 1.0
                lim = (np.min(eval(values,globals(),self.__dict__)[:,:,ipic]),np.max(eval(values,globals(),self.__dict__)[:,:,ipic]))
           
            plt.grid(True)
            if((values == "eps") | (values == "epsxx") | (values == "epszz")| (values == "epsxz")):
                cf=ax0.pcolormesh(xx,zz,eval(values,globals(),self.__dict__)[:,:,ipic]*cor,cmap='cmc.bilbao')
            else:
                cf=ax0.pcolormesh(xx,zz,eval(values,globals(),self.__dict__)[:,:,ipic]*cor,cmap='cmc.bilbao',vmin=lim[0],vmax=lim[1])
            
            cbar = fg.colorbar(cf,ax=ax0,orientation='horizontal')

            ax0.set_title(tick)
            ax0.tick_params(axis='both', which='major', labelsize=5)
            ax0.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
            plt.draw()    # necessary to render figure before saving
            fg.savefig(fn,dpi=300)
            ax0.plot()
            ic += 1
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
        for i in range(len(C.y)-1):
            for j in range(len(C.x)-1): 
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
         
class Phase_det(Phase):
    def __init__(self,C,Phase_dic):
        super().__init__(C,Phase_dic)
        ind_x = C.ind_xp
        ind_y = C.ind_yp
        ind_z = C.ind_zp
        tx    = np.sum(ind_x)
        ty    = np.sum(ind_y)
        tz    = np.sum(ind_z)
        self.tau   = np.zeros([tz,ty,tx],dtype = float)
        self.eps   = np.zeros([tz,ty,tx],dtype = float)
        self.T     = np.zeros([tz,ty,tx],dtype = float)
        self.nu    = np.zeros([tz,ty,tx],dtype = float)
        self.vis   = np.zeros([tz,ty,tx],dtype = float)
        self.Rho   = np.zeros([tz,ty,tx],dtype = float)
    def _interpolate_dyn_phase(self,V,C):
        # prepare the variables
        xp = C.xp
        yp = C.yp
        zp = C.zp
        x  = C.x
        y  = C.y
        z  = C.z 
        t1= perf_counter()
        val_ = ['tau','eps','T','nu','vis','Rho']
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
    for k in prange(len(zp)-1):
            for j in prange(len(yp)-1):
                for i in prange(len(xp)-1):
                    xx = xp[i]
                    yy = yp[j]
                    zz = zp[k]
                    ix = find1Dnodes(x,xx,len(x)-1)
                    iy = find1Dnodes(y,yy,len(y)-1)
                    iz = find1Dnodes(z,zz,len(z)-1)
                    x1 = x[ix]
                    x2 = x[ix+1]
                    y1 = y[iy]
                    y2 = y[iy+1]
                    z1 = z[iz]
                    z2 = z[iz+1]
                    intp1,intp2,intp3,intp4,intp5,intp6,intp7,intp8 = findnodes(buf,ix,iy,iz)
                    buf2[k,j,i]=trilinearinterpolation(xx,yy,zz,x1,x2,y1,y2,z1,z2,intp1,intp2,intp3,intp4,intp5,intp6,intp7,intp8)
    return buf2
    
    
    