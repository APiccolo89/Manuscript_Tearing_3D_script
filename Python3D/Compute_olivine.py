from Read_VTK_files_LAMEM import * #Read_VTK_files_LAMEM
from Read_VTK_files_LAMEM import  _file_list
from Parser_File import *
# Read the main file for processing the simulation
from Running_script_main_function import _run_script_visualization
from Slab_detachment import *


Depth = 190e3
g     = 9.81
rho   = 3300
Pr    = rho*Depth*g
T     = 1350+273.15
epsilon0 = 1e-15
RF = Rheological_data_Base()
dOl_dif=RF.Diffusion_DryOlivine
dOl_dis =RF.Dislocation_DryOlivine

B_d0 = dOl_dif.B*10e3**3

eta_D=0.5*(B_d0*10e3**(-3))**(-1)*np.exp((dOl_dif.E+Pr*dOl_dif.V)/(dOl_dif.R*T))

eta_dis=0.5*(dOl_dis.B)**(-1/dOl_dis.n)*np.exp((dOl_dis.E+Pr*dOl_dis.V)/(dOl_dis.n*dOl_dis.R*T))*epsilon0**((1-dOl_dis.n)/dOl_dis.n)

B_c = (0.5*(eta_dis)**(-1)*np.exp((dOl_dis.E+Pr*dOl_dis.V)/(dOl_dis.n*dOl_dis.R*T))*epsilon0**((1-dOl_dis.n)/dOl_dis.n))**dOl_dis.n

B_d = (0.5*(1e22)**(-1)*np.exp((dOl_dif.E+Pr*dOl_dif.V)/(dOl_dif.n*dOl_dis.R*T))*epsilon0**((1-dOl_dif.n)/dOl_dif.n))**dOl_dif.n

print('The value of computed B_d is %2e'%B_d)
print('The value of DBase B_d is %2e'%dOl_dif.B)
print('reference diffusion viscosity is %2e'%eta_D)



print('bla')






