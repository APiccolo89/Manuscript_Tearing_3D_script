
clear all;
close all;
addpath matlab   %Here -> Folder2LaMEM/matlab => all the function that handle matlab files are there
addpath(genpath('2D_numerical_setup'))
npart = [3,3,3];
% See model setup in Paraview 1-YES; 0-NO
Paraview_output        = 1;
% Output parallel files for LaMEM, using a processor distribution file (msetup = parallel)
LaMEM_Parallel_output  = 1;
Parallel_partition     = 'ProcessorPartitioning_8cpu_4.1.2.bin';
[A,Gr] = Parse_LaMEM_bin(Parallel_partition,Paraview_output,LaMEM_Parallel_output,npart);
%% Phase Lists
phases.Ph_Ar  = [0,10]  ;% Air
phases.Ph_UC  = [1,2700]  ;% Upper Crust
phases.Ph_LC  = [2,2750]  ;% Lower Crust
phases.Ph_Clt = [3,3300] ;% Continental Lithosphere
phases.Ph_Clt2 = [4,3300] ;
phases.Ph_UM  = [5,3300]  ;% Upper Mantle
phases.Ph_OLt = [6,3324] ;% Ocean Lithosphere
phases.Ph_OC  = [7,3324]   ;%place holder
phases.Ph_sed_oc = [8,2680];
phases.Ph_WZ  = [9,3300]  ;% Weak Zone
phases.Ph_LC2  = [10,2700]  ;
phases.Ph_UC2  = [11,2750] ;
phases.Ph_cont_pr = [12,2680];
phases.Ph_pas_m = [13,2680];
%% Stratigraphic structure list 
% Stratigraphies type
continental_stratigraphy.phases = [phases.Ph_UC(1),phases.Ph_LC(1),phases.Ph_Clt(1)];
continental_stratigraphy.Tk     = [0.0,-15.0,-30.0,-100.0];
oceanic_stratigraphy.phases = [phases.Ph_OLt(1)];
oceanic_stratigraphy.Tk     = [0.0, -80.0];
continental_stratigraphy2.phases = [phases.Ph_UC2(1),phases.Ph_LC2(1),phases.Ph_Clt2(1)];
continental_stratigraphy2.Tk     = [0.0,-15.0,-30.0,-100.0];
%% Thermal Type designed for the setup
Thermal_TypeTrench = Thermal_Type;
Thermal_TypeTrench.Age = convert_Age_velocity(80,1); 
Thermal_TypeTrench.vel = convert_Age_velocity(5,2); 
Thermal_TypeTrench.Type = 'McKenzie';
% 
Thermal_TypeContinent = Thermal_Type; 
Thermal_TypeContinent.Type = 'ContinentalGeotherm'; 
Thermal_TypeContinent.Moho = 600; 
Thermal_TypeContinent.Moho_d = -35; 
%% General Information for the 
Thermal_information = Thermal_Information;
Thermal_information.TP = 1350; 
Thermal_information.TS = 20; 
Thermal_information=Thermal_information.kappa_calculation; 
Thermal_information.Ph_Air = 0; 
Thermal_information.Ph_Ast = 5; 
% List Trech
Continent1     = Terrane;
Continent2     = Terrane;
passive_margin1 = Passive_Margin; 
% Special Boundary
T = Trench; 
% Fill the properties of the terrane template
Continent1.Boundary = [-1000,min(Gr.y_g),0,max(Gr.y_g)];
Continent1.Stratigraphy=continental_stratigraphy; 
Continent1.Thermal_type = Thermal_TypeContinent; 
Continent1.Thermal_information = Thermal_information; 
% Fill the properties of the terrane template 
Continent2.Boundary = [0,min(Gr.y_g),1000,max(Gr.y_g)];
Continent2.Stratigraphy=continental_stratigraphy2; 
Continent2.Thermal_type = Thermal_TypeContinent; 
Continent2.Thermal_information = Thermal_information; 
% Fill the properties of the subduction zone
T.Boundary = 0.0;
T.Stratigraphy_Continental = continental_stratigraphy;
T.Stratigraphy_Oceanic = oceanic_stratigraphy;
T.theta=30;
T.theta_c = 20;
T.D0   = 80; 
T.L0   = 300; 
T.Decoupling_depth = -100; 
T.Thermal_type = Thermal_TypeTrench;
T.Thermal_information = Thermal_information; 
T.R = [120,T.R-T.D0]; 
T.Type_Subduction = 'Mode_1'; 
T.tk_WZ = 0.0; 

Terranes = struct('Continent1',Continent1,'Continent2',Continent2,'T',T,'Passive_Margin',passive_margin1);
%% Generic information numerical domain:

Create_Setup(Terranes,phases,Thermal_information,A,npart,Gr,Parallel_partition);
TB = cputime;
disp('====================================================================')
disp(['The Setup is finished and took ', num2str(round((TB-TA)./60)), ' min'])
disp('====================================================================')



















function [converted] = convert_Age_velocity(value,type)
        Myrs_sec = 365.25*60*60*24*1e6;
        cm_to_meter = 100; 
        if type == 1
            converted = value*Myrs_sec; 
        elseif type ==2 
            converted = value/(Myrs_sec/1e6)/cm_to_meter; 
        end
end