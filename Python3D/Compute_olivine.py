from Read_VTK_files_LAMEM import * #Read_VTK_files_LAMEM
from Read_VTK_files_LAMEM import  _file_list
from Parser_File import *
# Read the main file for processing the simulation
from Running_script_main_function import _run_script_visualization
from Slab_detachment import *


Depth = 400e3
g     = 9.81
rho   = 3300
Pr    = rho*Depth*g
T     = 1350+273.15
epsilon0 = 1e-15
RF = Rheological_data_Base()
dOl_dif=RF.Diffusion_DryOlivine
dOl_dis =RF.Dislocation_DryOlivine
dOl_dis.V =0.8e-05


B_d0 = dOl_dif.B*10e3**3

eta_D=0.5*(B_d0*10e3**(-3))**(-1)*np.exp((dOl_dif.E+Pr*dOl_dif.V)/(dOl_dif.R*T))

eta_dis=0.5*(dOl_dis.B)**(-1/dOl_dis.n)*np.exp((dOl_dis.E+Pr*dOl_dis.V)/(dOl_dis.n*dOl_dis.R*T))*epsilon0**((1-dOl_dis.n)/dOl_dis.n)

B_c = (0.5*(eta_dis)**(-1)*np.exp((dOl_dis.E+Pr*dOl_dis.V)/(dOl_dis.n*dOl_dis.R*T))*epsilon0**((1-dOl_dis.n)/dOl_dis.n))**dOl_dis.n

B_d = (0.5*(1e22)**(-1)*np.exp((dOl_dif.E+Pr*dOl_dif.V)/(dOl_dif.n*dOl_dis.R*T))*epsilon0**((1-dOl_dif.n)/dOl_dif.n))**dOl_dif.n

print('The value of computed depth is %2e km'%(Depth/1e3))
print('The value of DBase P_r is %2e MPa'%(Pr/1e6))
print('The value of DBase Vn is %2e'%dOl_dis.V)
print('reference diffusion viscosity is diffusion:%2e and dislocation: %2e and the effective viscosity is %2e'%(eta_D,eta_dis,(1/eta_D+1/eta_dis)**(-1)))



print('bla')






