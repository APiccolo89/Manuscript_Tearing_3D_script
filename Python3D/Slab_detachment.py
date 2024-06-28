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
from numba.experimental import  jitclass
from numba import jit, types, typed, prange
from Read_VTK_files_LAMEM import * #Read_VTK_files_LAMEM
import scipy.ndimage as ndi
from Read_VTK_files_LAMEM import  _file_list
from Read_VTK_files_LAMEM import findnodes
from Read_VTK_files_LAMEM import bilinearinterpolation
from Read_VTK_files_LAMEM import trilinearinterpolation
from Read_VTK_files_LAMEM import linearinterpolation
import copy 

from Parser_File import * 

    

class SLAB():
    def __init__(self,C,IG,nstep):
        Slab_Geometry = IG.Slab
        self.P1 = Slab_Geometry.Boundary_B[0][0]
        self.P2 = Slab_Geometry.Boundary_B[0][2]
        t= (C.xp>=self.P1)&(C.xp<=self.P2)
        self.ind_boundary = np.where(t==True)
        tx  = len(self.ind_boundary[0])
        self.y_b,self.x_s =_compute_length_coordinatesCB(self.ind_boundary,Slab_Geometry.Boundary_B,C.xp)
        tz = len(C.zp)
        self.D     = np.ones((tx,tz,nstep),dtype=float)*(-1.0)
        self.T     = np.zeros((tx,tz,nstep),dtype=float)
        self.dDdt  = np.zeros((tx,tz,nstep),dtype=float)
        self.eps   = np.zeros((tx,tz,nstep),dtype=float)
        self.tau   = np.zeros((tx,tz,nstep),dtype=float)
        self.tau_max   = np.zeros((tx,tz,nstep),dtype=float)
        self.vis   = np.zeros((tx,tz,nstep),dtype=float)
        self.F_T   = np.zeros((tx,tz,nstep),dtype=float)
        self.Rho  = np.zeros((tx,tz,nstep),dtype=float)
        self.Psi = np.zeros((tx,tz,nstep),dtype=float)
        self.L   = np.zeros((tx,tz,nstep),dtype=float)
        self.det_vec =np.ones((tx),dtype=float)*np.nan
        self.y_vec =np.ones((tx),dtype=float)*np.nan
        self.x_vec =np.ones((tx),dtype=float)*np.nan
        self.x1    = np.zeros((tx,tz,nstep),dtype=float)
        self.x2    = np.zeros((tx,tz,nstep),dtype=float)
        self.W    = np.zeros((tz,nstep),dtype = float)
        self.tau_vec = np.ones((tx),dtype=float)*np.nan
        self.T_vec = np.ones((tx),dtype=float)*np.nan
        self.depth_vec = np.zeros((tx),dtype=float)
        self.dDdt_vec =  np.zeros((tx,nstep),dtype=float)
        self.nodes_tearing_=np.ones(nstep,dtype=int)*np.nan
        self.average_tearing_velocity=np.ones((3,nstep),dtype=float)*np.nan
        self.LGV         = ["D","tau_max",'eps','vis','Psi','T']
        self.Label       = ["$D^{\dagger} []$",
                            r"$\tau_{II,max},[MPa]$",
                            r"$\dot{\varepsilon}_{II,mean} [1/s]$",
                            r"$\eta_{S}, [Pas]$",
                            r"$\Psi, [W/m^3]$",
                            r"$T, [^/circ C]$"
                            ]
        self.Colormap    = ["cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.oleron","cmc.oleron","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao","cmc.bilbao"]
        self.Val         = [(0.1,0.85),
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
        self.CV = ["T","eps","tau",'Rho','Psi',"dDdt","tau_max","vis"]       


    def _update_C(self,C,FS,Ph,IG,ipic,tcur,dt):
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
                switch = 0.0
                # Prepare buffer variable to fill with the new one
                if iv != 'dDdt':
                    if iv == 'tau_max':
                        switch = 1
                    if iv != 'tau_max':
                        buf_var_ph = eval(iv,globals(),Ph.__dict__)[:,:,ix]
                    else:
                        buf_var_ph = Ph.tau[:,:,ix]
                        
                    D     = np.zeros((np.sum(C.ind_zp)),dtype=float)
                    buf_var   = np.zeros(np.sum(C.ind_zp),dtype=float)
                    z_bottom   = np.zeros(np.sum(C.ind_zp),dtype=float)
                    x1    = np.zeros(np.sum(C.ind_zp))
                    x2    = np.zeros(np.sum(C.ind_zp))
                    # Run the function 
                    D,L0,buf_var,x1,x2 = _Find_Slab_PERFORM_C(C.yp,C.zp,buf_var_ph,lay_ph[:,:,ix],np.sum(C.ind_zp),ix,D,buf_var,x1,x2,z_bottom,switch)
                    # update the class 
                    self.D[i,:,ipic] = D
                    self.L[i,:,ipic] =L0
                if iv == 'dDdt':
                    if ipic == 0:
                        self.dDdt[i,:,ipic]=D*0.0
                    else:
                        buf = (((self.D[i,:,ipic]-self.D[i,:,ipic-1])*1e3*1e2)/(dt*1e6))
                        buf[buf>=0.0] = 1.0
                        buf[buf<0.0] = np.abs(buf[buf<0.0])
                        self.dDdt[i,:,ipic] = buf
                else:
                    eval(iv,globals(),self.__dict__)[i,:,ipic] = buf_var
                
                self.x1[i,:,ipic] = x1
                self.x2[i,:,ipic] = x2
            # Compute the additional variable (i.e. F_T = 2*D*tau), F_B
            self.F_T[i,:,ipic] = 2*self.D[i,:,ipic]*1e3*self.tau[i,:,ipic]*1e6
        # detect_slab_detachment
        ix1 = np.zeros((len(C.zp)),dtype=int)
        ix2 = np.zeros((len(C.zp)),dtype=int)
        
        ix1,ix2,self.W[:,ipic]=compute_W_slab(self.D[:,:,ipic],C.xp[self.ind_boundary[0]],C.zp,self.W[:,ipic],ix1,ix2)
        self.det_vec,self.depth_vec,self.T_vec,self.tau_vec, self.dDdt_vec[:,ipic] =detect_slab_detachment(self.D[:,:,ipic],
                                                                                    C.xp[self.ind_boundary[0]],
                                                                                    C.zp,
                                                                                    tcur,
                                                                                    ipic,
                                                                                    self.det_vec,
                                                                                    self.tau_vec,
                                                                                    self.depth_vec,
                                                                                    self.T_vec,
                                                                                    self.T,
                                                                                    self.tau,
                                                                                    ix1,ix2,
                                                                                    self.x1,
                                                                                    self.y_vec,
                                                                                    self.x_vec,
                                                                                    self.dDdt_vec[:,ipic],
                                                                                    self.D[:,:,ipic])
        # deep copy the detachment vector (variables space = pointer space, it means that if i do not deep copy a variable, it modifies it and screw up everything)
        det_buf = np.copy(self.det_vec)
        tearing_ = np.isnan(det_buf)==False
        nodes_   = np.sum(tearing_)
        self.nodes_tearing_[ipic] = nodes_-1
        if ipic == 0:
            self.average_tearing_velocity[0,ipic] = 0.0
            self.average_tearing_velocity[1,ipic] = 0.0
            self.average_tearing_velocity[2,ipic] = 0.0

        else:
            dx_s = np.diff(self.x_s)
            dx_s_m = np.min(dx_s)
            dx_s_av = np.mean(dx_s)
            dx_s_M = np.max(dx_s)
            self.average_tearing_velocity[0,ipic] = ((self.nodes_tearing_[ipic]-self.nodes_tearing_[ipic-1])/dt)*dx_s_m
            self.average_tearing_velocity[1,ipic] = ((self.nodes_tearing_[ipic]-self.nodes_tearing_[ipic-1])/dt)*dx_s_av
            self.average_tearing_velocity[2,ipic] = ((self.nodes_tearing_[ipic]-self.nodes_tearing_[ipic-1])/dt)*dx_s_M
            print('minumum tearing velocity is:', "{:03}".format(self.average_tearing_velocity[0,ipic]), 'km/Myr' )
            print('average tearing velocity is:', "{:03}".format(self.average_tearing_velocity[1,ipic]), 'km/Myr' )
            print('maximum tearing velocity is:', "{:03}".format(self.average_tearing_velocity[2,ipic]), 'km/Myr' )
            self.average_tearing_velocity[:,ipic] = self.average_tearing_velocity[:,ipic]*(1e5)/(1e6)
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
        self.RB         = vIC.Slab.Boundary_B[2][2]
        self.T_av       = vIC.Slab.avTS+273.15
        string = 'Average temperature of the slab is: %2f degC' %(vIC.Slab.avTS)
        print(string)
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
        B_A_Buff = (''.join((str(np.string_(i))[2]) for i in np.array(mat[ind])[:]))
        if B_A_Buff == 'none':
            B_A.append([])
            B_A.append(np.array([np.array(0),np.array(0),np.array(0)]))

        else:
            B_A.append([1])
            ind = np.array(A)[2][0]
            B_A.append(np.array([np.array(mat[ind])[0][0], np.array(mat[ind])[1][0],np.array(mat[ind])[2][0]]))
            # Correction with the most recent python version 
            B_A = np.array(B_A,dtype=object)
           
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
        self.vx_M        = np.zeros((ty,tx,nstep),dtype=float)
        self.vy_M        = np.zeros((ty,tx,nstep),dtype=float)
        self.vz_M        = np.zeros((ty,tx,nstep),dtype=float)
        self.mean_stress = np.zeros((ty,tx,nstep),dtype=float)
        self.thickness = np.zeros((ty,tx,nstep),dtype=float)
        self.mean_eps = np.zeros((ty,tx,nstep),dtype=float)
        self.x_s         = np.ones(tx)*np.nan
        self.x_sp         = []

        self.ind_boundary = np.zeros(tx)*np.nan
          
    def _update_extra_variables(self,V:VAL,C:Coordinate_System,dt,ipic):
        self.dH[:,:,ipic]          = self._update_FS(dt,ipic,'dH')
        self.vz_M[:,:,ipic]        = self._update_FS(dt,ipic,'vz_M')
        self.vx_M[:,:,ipic]        = self._update_FS(dt,ipic,'vx_M')
        self.vy_M[:,:,ipic]        = self._update_FS(dt,ipic,'vy_M')
        self.mean_stress[:,:,ipic] = self.update_Mean_CC(V,C,ipic,'tau')
        self.mean_eps[:,:,ipic]    = self.update_Mean_CC(V,C,ipic,'eps')
        self.thickness[:,:,ipic] = self.update_Mean_CC(V,C,ipic,'thickness')
        self.LGV = ['dH','vz_M','mean_stress','mean_eps','Amplitude','thickness','F_z','vy_M']
        self.Label  = [r'$\dot{H}, [\frac{mm}{yr}]$',r'$v_z[\frac{mm}{yr}]$',r'$\bar{\tau}_{II} [MPa]$',r'$\bar{\dot{varepsilon}}_{II} [\frac{1}{s}]$','Topography, $[km]$','Lithosphere $[km]$',r'$F_z, [\frac{N}{m}]$','$v_y$ [$cm/yrs$]']
        self.Colormap = ["cmc.cork","cmc.cork","cmc.bilbao","cmc.devon","cmc.oleron","cmc.bilbao","cmc.hawaii","cmc.hawaii"]
        self.Val         = [(-10,10),
                            ("min","max"),
                            ("min","max"),
                            (1e-16,1e-13),
                            ("min","max"),
                            ("min","max"),
                            ("min","max"),
                            (-10,10),
			    (-10,10)]
        
        return self
        
        
    def _update_FS(self,dt,ipic,name_field):
        if name_field == 'dH':
            if ipic > 0:
                buf = (self.Topo[:,:,ipic]-self.Topo[:,:,ipic-1])/dt
                # Convert into mm/yr [all surface velocity must be scaled using mm/yrs]
                buf = buf*(1e6/1e6)
            else:
                buf = self.dH[:,:,ipic]*0.0
        else:
            if name_field == 'vx_M':
                buf_pr = self.vx
            elif name_field == 'vy_M':
                buf_pr = self.vy 
            else:
                buf_pr = self.vz
            if ipic > 0:
                buf = (buf_pr[:,:,ipic]+buf_pr[:,:,ipic-1])/2 # cm/yr
                buf = buf*10 # mm/yr
            else: 
                buf = buf_pr[:,:,ipic]
        return buf 
    
    def update_Mean_CC(self,V:VAL,C:Coordinate_System,ipic,field_name):
        # Compute the continental crust
        Continental_crust                           = V.Lit+V.OP 
        Continental_crust[Continental_crust>=0.8]      = 1.0
        Continental_crust[Continental_crust<0.8]    = np.nan
        buf = 0.0*Continental_crust
        buf = buf[0,:,:]
        # loop over x,y. Select colum: compute mean of those that has high crustal fraction 
        for i in range(len(C.y)):
            for j in range(len(C.x)): 
                crust  = Continental_crust[:,i,j]
                if len(crust[crust==1]) == 0:
                    buf[i,j]=np.nan 
                else:
                    if field_name == 'thickness':
                        buf[i,j]=np.abs(np.nanmax(C.z[(crust==1) & (C.z>-100)])-np.nanmin(C.z[(crust==1) & (C.z>-100)]))
                    else:
                        column = eval(field_name,globals(),V.__dict__)[:,i,j]
                        buf[i,j] = np.mean(column[(crust==1) & (C.z>-100)])
                    
        return buf 

#    def _compute_relevant_information_topography(self,C: Coordinate_System,Slab_GEO:Trench,ipic:int):
#        # Compute the plate boundary using the data of the boundary 
#            """
#            Approach: compute x-y of the boundary of the terrane
#            -> interpolate topography at the boundary of the terrane and other variables
#            -> loop over the nodes of the point of the margin:
#                -> select the area between p(B)-200,p(B)+200 
#                    -> find the maximum topography (store its node)
#                    -> find the minimum topography:
#                        -> behind and in front the maximum topography 
#                    -> compute the averages, the wavelength (i.e. distance between the two minima)
#            """
#            # Extract information boundary and select x belongs to xa-xb 
#            boundary_geometry = Slab_GEO.Boundary_B
#            xa = boundary_geometry[0][0]
#            xb = boundary_geometry[0][2]
#            ind_boundary = (C.x>=xa) & (C.x<xb)
#            x_B = C.x[ind_boundary]
#            # Compute the boundary 
#            y_B,x_s =_compute_length_coordinatesCB(ind_boundary,boundary_geometry,C.x)
#            # function to retrieve the data. 
#            vz = self.vz_M[:,:,ipic]
#            tau_C = self.mean_stress[:,:,ipic] 
#            H  = self.Amplitude[:,:,ipic]
#            self.HBx[ind_boundary==True,ipic]=x_B[:]
#            self.HBy[ind_boundary,ipic]=y_B[:]
#            self.x_s = x_s
#            ind_boundary_p = (C.xp>=xa) & (C.xp<xb)
#            x_p = C.xp[ind_boundary_p]
#            y_Bp,x_sp =_compute_length_coordinatesCB(ind_boundary_p,boundary_geometry,C.xp)
#            self.x_sp = x_sp
#            
#            """
#            HB     topography of plate boundary  
#            HBc    coordinate of the plate boundary
#            HMax   topography of the maximum   
#            HMaxc  coordinate maximum
#            HminF  topography frontal minumum
#            HminFc topography frontal minimum coordinates
#            HminB  topography back minimum
#            Hminc  coordinate topography back mininum
#            v_z_M  max velocity     [HBc+/-200]
#            v_z_m  minumum velocity [HBc+/-200]
#            x_s    length from left to right 
#            """
#            ix_b = np.sum(ind_boundary)
#            ix_ch = np.where(ind_boundary==True)
#            for i_x in range(ix_b):
#                i = ix_ch[0][i_x]
#                yy = y_B[i_x] # 
#                iy = find1Dnodes(C.y,yy,len(C.y))
#                y1 = C.y[iy]
#                y2 = C.y[iy+1]
#                intp1 = H[iy,i]
#                intp2 = H[iy,i]
#                self.HB[i,ipic] = linearinterpolation(yy,y1,y2,intp1,intp2)
#                # find area of interest: 
#                ind_area = np.where((C.y<=self.HBy[i,ipic]+200)&(C.y>=self.HBy[i,ipic]-200))
#                # find_maximum
#                self.HMax[i,ipic] = np.max(H[ind_area,i])
#                ind_max   = np.where(H[:,i]==self.HMax[i,ipic])
#                # find_coordinate
#                if C.y[ind_max[0][0]]<self.HBy[i,ipic]-200:
#                    self.HMaxy[i,ipic]=np.nan
#                else:
#                    self.HMaxy[i,ipic]=C.y[ind_max[0][0]]
#                self.Hmeang[ipic] = np.mean(H[ind_area[0][0],ind_boundary==True])
#                # find_minimum_front
#                if self.HMaxy[:,ipic].all() == False:
#                    self.Hmin[i,ipic]  = np.min(H[ind_area,i])
#                    self.HMean[i,ipic]  = np.mean(H[ind_area,i])  
#                    self.Hminy[i,ipic] = C.y[H[:,i]==np.min(H[ind_area,i])]
#                else:
#                    self.Hmin[i,ipic]   = np.nan
#                    self.Hminy[i,ipic]  = np.nan
#                    self.HMean[i,ipic]  = np.mean(H[ind_area,i])  
#                
#                self.v_z_M[i,ipic] = np.max(vz[ind_area,i])
#                self.v_z_m[i,ipic] = np.min(vz[ind_area,i])
#                self.v_z_mean[i,ipic] = np.mean(vz[ind_area,i])
#                self.mean_stress_1D[i,ipic] = np.mean(tau_C[ind_area,i])
#            self.ind_boundary=ind_boundary
#            return self 
#    def _plot_1D_plots_Free_surface(self,ipic: int,ptsave,S:SLAB,t_cur,D0):
#        val = ['HMax','v_z']
#        ptsave=os.path.join(ptsave,'1D_surfaces')
#        if not os.path.isdir(ptsave):
#            os.mkdir(ptsave)
#        ib = self.ind_boundary
#        for v in val:
#            time_sim = "{:.3f}".format(t_cur)
#            tick = r"t = %s [Myrs]" %(time_sim)
#            ptsave_b=os.path.join(ptsave,v)
#            if not os.path.isdir(ptsave_b):
#                os.mkdir(ptsave_b)
#            fig = plt.figure()
#            ax2 = fig.gca()
#            fna='Fig'+str(ipic)+'.png'
#            fn=os.path.join(ptsave_b,fna)
#            if v=='HB':
#                for ip in range(10):
#                    it = ipic - ip
#                    alpha_v= 0.8-ip*(1/12)
#                    if ip == 0: 
#                        cc = 'r'
#                    else:
#                        cc = 'b'
#                    if (it == 0) & (ip == 0) :
#                        ax2.plot(self.x_s, self.HMax[ib==True,0]-self.Hmeang[0],c = cc,alpha = alpha_v,linewidth=alpha_v)
#                        break
#                    if (ip >0 ) & (it == 0 ):
#                        ax2.plot(self.x_s, self.HMax[ib==True,it]-self.Hmeang[it],c = cc,alpha = alpha_v,linewidth=alpha_v)
#                        break 
#                    else: 
#                        ax2.plot(self.x_s, self.HMax[ib==True,it]-self.Hmeang[it],c = cc,alpha = alpha_v,linewidth=alpha_v)
#                ax2.plot(self.x_s, self.HB[ib==True,ipic]-self.Hmean[ipic],c = 'r',alpha = 1.0,linewidth=1.2)
#                plt.xlabel('x, [km]')
#                plt.ylabel('H, [km]')
#                plt.xlim(0,1200)
#
#            elif v == 'v_z':
#                p1=ax2.plot(self.x_s,self.v_z_M[ib==True,ipic],c = 'b',alpha = 1.0,linewidth=0.6)
#                p2=ax2.plot(self.x_s,self.v_z_m[ib==True,ipic],c = 'b',alpha = 1.0,linewidth=0.6)
#                p3 =ax2.fill_between(self.x_s, self.v_z_m[ib==True,ipic], self.v_z_M[ib==True,ipic],color='blue',alpha=0.2)
#                p4=ax2.plot(self.x_s,self.v_z_mean[ib==True,ipic],c = 'r',alpha = 1.0,linewidth=1.2)
#                ax2.set_xlabel(r'$\ell_{trench}, [km]$')
#                ax2.set_ylabel(r'$\dot{H}$, $[\frac{cm}{yrs}]$')
#                ax3 = ax2.twinx()
#                D_ = ndi.uniform_filter1d(S.dDdt_vec[:,ipic], size=20)
#                p5=ax3.plot(self.x_sp,D_/D0[0][0],c = 'k',alpha = 1.0,linewidth=1.2,linestyle='dashed')
#                ax3.set_ylim(0.1,1.0)
#                ax3.set_xlim(0.0,1200.0)
#                ax3.set_ylabel(r'$D^{\dagger}$, $[]$')
#                ax3.set_ylim()
#                fig.tight_layout()  # otherwise the right y-label is slightly clipped
#            elif v == 'HMax':
#                p1=ax2.plot(self.x_s,self.HMax[ib==True,ipic]-self.Hmeang[ipic],c = 'r',alpha = 1.0,linewidth=1.0,linestyle=':')
#                p2=ax2.plot(self.x_s,self.Hmin[ib==True,ipic]-self.Hmeang[ipic],c = 'b',alpha = 1.0,linewidth=1.2)
#                ax2.set_xlabel(r'$\ell_{trench}, [km]$')
#                ax2.set_ylabel(r'${H}$, $[\frac{cm}{yrs}]$')
#                plt.xlabel('x_s, [km]')
#                plt.ylabel('H, [km]')
#                plt.xlim(0,1200)
#            elif v == 'Maps':
#                p1=ax2.fill_between(self.HBx[ib==True,ipic], self.Hminy[ib==True,ipic],self.HMaxy[ib==True,ipic],color='blue',alpha=0.4)
#                p3=ax2.plot(self.HBx[ib==True,ipic],self.HBy[ib==True,ipic],c = 'k',alpha = 1.0,linewidth=0.8,linestyle=':')
#                plt.xlabel('x, [km]')
#                plt.ylabel('y, [km]')
#            p5 = ax2.scatter(self.x_sp,np.zeros(len(self.x_sp)),10,S.det_vec,cmap='inferno')
#            try: 
#                lim1_m = np.nanmin(S.det_vec)
#                lim2_m = np.nanmax(S.det_vec)
#            except: 
#                lim1_m = 0
#                lim2_m = 1
#            cbar2 = fig.colorbar(p5, ax=ax2,orientation='horizontal',extend="both")
#            p5.set_clim([lim1_m,lim2_m])
#            cbar2.vmin = lim1_m 
#            cbar2.vmax = lim2_m
#            ax2.set_title(tick)
#
#            #plt.grid(True)
#            #plt.xlabel(r'\ell_{trench}, [km]')
#            #plt.ylabel('H, [km]')
#            ax2.tick_params(axis='both', which='major', labelsize=5)
#            ax2.tick_params(axis='both',bottom=True, top=True, left=True, right=True, direction='in', which='major')
#        ###############################################################       
#            plt.draw()    # necessary to render figure before saving
#            fig.savefig(fn,dpi=300,transparent=False)
#            fig.clear()
#            plt.close()

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
def _Find_Slab_PERFORM_C(x,z,buf_var_ph,ph,tz,ix,D,buf_var,x1,x2,z_bottom,switch):    #(self,C,Values,ph,ipic,ix,ix1):

    for i in range(tz):#prange(tz):
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
            i2 = 0
            i1,i2 = _find_index(buf,i1,i2)
            x1[i]= x[i1]
            x2[i] = x[np.int32(i2)]
            D[i]= x2[i]-x1[i] 
            if(x2[i]-x1[i] <8.0):
                buf_var[i] = np.nan
                D[i]=np.nan               
            else:
                if switch == 0.0:
                    buf_var[i] = _mean_numba(buf,buf2)
                else: 
                    buf_var[i] = _max_numba(buf,buf2)
    L0  =z-np.min(z_bottom)
    return D,L0,buf_var,x1,x2
    
@jit(nopython=True)  
def _mean_numba(buf,buf2):
    mean = 0.0
    len_ = np.int32(0)
    for i in range(len(buf)):
        if buf[i]>0.95:
            mean+= buf2[i]
            len_ += np.int32(1)
    mean = mean/len_
    
    return mean
@jit(nopython=True)  
def _max_numba(buf,buf2):
    max = 0.0
    for i in range(len(buf)):
        if (buf[i]>0.95) & (max<buf2[i]):
            max = buf2[i]    
    return max

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
    for ib in range(len(buf)-1):
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

    for ib in range(len(buf)-1,-1,-1):
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

@jit(nopython=True,parallel=True)
def detect_slab_detachment(D,x,z,tcur,ipic,det_vec,tau_vec,depth_vec,T_vec,T,tau,ind1,ind2,x1,yvec,xvec,dDdt_vec,dDdt):
        
        """
        Input:
        =============================================================================================
        D => D[ix,iz,ipic] the thickness of the slab of the current timestep.{float}
        x => x[] the x vector{float}
        z => z vector {float}
        t_cur => current time {float}
        ipic => actual timestep 
        det_vec => vector containing the position w.r.t the vector x, and time of the first detachment
        tau_vec => vector containing the stress associated with the previous timestep
        depth_vector => vector that contains the depth of the first occurence of detachment for a given x 
        T_vec => vector containing the temperature shortly before the detachment 
        T=> T[ix,iz,ipic] => temperature saved from the detection of the slab
        tau=> "" => tau saved from the detection of the slab
        ind1 => lateral boundary of the slab 
        ind2 => lateral boundary of the slab of the current timestep
        ==================================================================================================
        Output: (vector that are used to update the information within the class of slab detachment)
        det_vec
        taz_vec
        depth_vec
        T_vec
        """
        
        for ix in prange(len(x)):
            dDdt_vec[ix]=200
            for i in range(len(z)-1):              #Loop starting from the top to the bottom 
                buf = D[ix,i]                   #1D(z) vector of D 
                if (z[i] <-80) & (z[i]>-300):      #Constraining better the area of research
                    if (buf<dDdt_vec[ix]) & (np.isnan(buf)==0)& (np.isnan(det_vec[ix])):
                        dDdt_vec[ix]=buf
                    if (np.isnan(buf)) & (ix > ind1[i]) &  (ix < ind2[i]) & (np.isnan(det_vec[ix])):
                        det_vec[ix] = tcur
                        depth_vec[ix]=z[i]
                        T_vec[ix]    =T[ix,i,ipic-1]
                        tau_vec[ix]  = tau[ix,i,ipic-1]
                        xvec[ix]    = x[ix]
                        yvec[ix]    = x1[ix,i,ipic-1]
                        dDdt_vec[ix] = np.nan
        return det_vec,depth_vec,T_vec,tau_vec,dDdt_vec
    
    
@jit(nopython=True)  
def compute_W_slab(D,x,z,W,ix1,ix2):
    """
    _summary_: function that compute the Width of the slab, and collects the indeces where to look after the detachment
    D = current thicnkess x(y)-z
    x,z = coordinate
    W the width vector (initialize to be zero)
    ix1,ix2 indexes. 
    The algorithm has been co-ideated with Arne (i.e., he gave me the right tip to make it efficient)
    """
    for iz in range(len(z)):
        ix1[iz] = -1
        x1   = np.nan
        buf = D[:,iz]
        for ix in range(len(x)):
            if (np.isnan(buf[ix])== False) & (buf[ix] != -1.0) & (buf[ix]!=-np.inf) :
                ix1[iz] = ix
                x1  = x[ix]
                break
            if ix == (len(x)):
                ix1[iz]=-int(1)
                x1 = np.nan
        ix2[iz] = -1
        x2   = np.nan
        if x1 != np.nan:
            for ixb in range(len(x)-1,-1,-1):
                if (np.isnan(buf[ixb])== False) & (buf[ixb] != -1.0) & (buf[ixb]!=-np.inf) :
                    ix2[iz] = ixb
                    x2  = x[ixb]
                    break
        else:
            ix2[iz] = -int(1)
            x2   = np.nan
        W[iz] = x2-x1
    """
    Short comments: most of the time that I spent for writing this function, was literally
    to find idiotic mistake, that I could have spotted if I stopped to think about them. 
    I had to introduce several indeces, I discovered that x = np.nan is not always working, and that could
    create future issue. Moreover, I needed to try the effects of range (yet again) to avoid additional future problem.
    Just a paragraph to self-punish me, and avoid that people that read this code consider me a worthwile asset. 
    """
    return ix1,ix2,W

# Passive Tracer routine for the current project 
class Basement_Passive_Tracer():
    def __init__(self,P:Passive_Tracers,C: Coordinate_System,F:FS,list_phases,levels,ts,ipic):
        self.ID_chosen = self.select_chosen(P,F,C,list_phases,levels,ipic)
        n_chosen = len(self.ID_chosen)
        self.P = np.zeros([ts,n_chosen],dtype = float)
        self.T = np.zeros([ts,n_chosen],dtype = float)
        self.x = np.zeros([ts,n_chosen],dtype = float)
        self.y = np.zeros([ts,n_chosen],dtype = float)
        self.z = np.zeros([ts,n_chosen],dtype = float)
        self.x[ipic,:] = P.x[self.ID_chosen]
        self.y[ipic,:] = P.y[self.ID_chosen]
        self.z[ipic,:] = P.z[self.ID_chosen]
        self.P[ipic,:] = P.P[self.ID_chosen]
        self.T[ipic,:] = P.T[self.ID_chosen]
        
        
    def select_chosen(self,P:Passive_Tracers,F:FS,C: Coordinate_System,list_phases:int,levels:float,ipic:int): 
        # Select the marker that belongs to the chosen phases 
        # From them, select the one that are within level[0] and level[1]
        # Extract the ID and send out
        ID_chosen = []
        Ph_lay = copy.deepcopy(P.Ph) # Deep copy the variable, as I do not want to modify the phase layer
        
        for il in list_phases:
            Ph_lay[P.Ph==il] = -100
        Ph_lay[Ph_lay>0] = 0
        Ph_lay[Ph_lay<0] = 1 
        # interpolate topography
        topo_marker = _interpolate_topography(F.Topo[:,:,ipic],C.x,C.y,P.x,P.y)
        # Select the ID_Chosen
        Ph_lay[(Ph_lay == 1) & (P.z>=topo_marker+levels[0]) & (P.z<=topo_marker+levels[1]) & (P.y<50)]=2
        ID_Chosen = np.where(Ph_lay == 2)
        ID_Chosen = ID_Chosen[0]
        print(len(ID_Chosen),' marker choosen')
        return ID_Chosen
    def _update_PTDB(self,P: Passive_Tracers,ipic:int,time:float):
        self.T[ipic,:] = P.T[self.ID_chosen]
        self.P[ipic,:] = P.P[self.ID_chosen]
        self.x[ipic,:] = P.x[self.ID_chosen]
        self.y[ipic,:] = P.y[self.ID_chosen]
        self.z[ipic,:] = P.z[self.ID_chosen]

    
# Decorate with numba and parallel     
@jit(nopython=True,parallel=True)
def _interpolate_topography(Topo:float,xg:float,yg:float,xp:float,yp:float):
    topo_marker = xp*0.0 
    for i in prange(len(xp)):
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
