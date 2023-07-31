
clear all;
close all;
addpath matlab   %Here -> Folder2LaMEM/matlab => all the function that handle matlab files are there
addpath(genpath('2D_numerical_setup'))
npart = [3,3,3];
% See model setup in Paraview 1-YES; 0-NO
Paraview_output        = 1;
% Output parallel files for LaMEM, using a processor distribution file (msetup = parallel)
LaMEM_Parallel_output  = 1;
% = = = = = = = = = = = = Script initial setup = = = = = = = = = = = = = =
Parallel_partition     = 'ProcessorPartitioning_8cpu_4.1.2.bin';
[A,Gr] = Parse_LaMEM_bin(Parallel_partition,Paraview_output,LaMEM_Parallel_output,npart);
%========================================================================== 
% Phase Lists
% Phase structure containing the phase ID number and reference density
%==========================================================================
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
%==========================================================================
% Stratigraphic structure list 
% Stratigraphies types ~ Structure containing the number of phases and
% layer to introduce in a specific geological unit
%==========================================================================
continental_stratigraphy.phases = [phases.Ph_UC(1),phases.Ph_LC(1),phases.Ph_Clt(1)];
continental_stratigraphy.Tk     = [0.0,-15.0,-30.0,-100.0];
oceanic_stratigraphy.phases = [phases.Ph_OLt(1)];
oceanic_stratigraphy.Tk     = [0.0, -80.0];
continental_stratigraphy2.phases = [phases.Ph_UC2(1),phases.Ph_LC2(1),phases.Ph_Clt2(1)];
continental_stratigraphy2.Tk     = [0.0,-15.0,-30.0,-100.0];
continental_stratigraphy_s.phases =[phases.Ph_pas_m(1),phases.Ph_UC(1),phases.Ph_LC(1),phases.Ph_Clt(1)];
continental_stratigraphy_s.Tk     = [0.0,-12.0,-15.0,-30.0,-100.0];

%==========================================================================
% Thermal types for the specific setup
%==========================================================================

Thermal_TypeTrench = Thermal_Type;
Thermal_TypeTrench.Age = convert_Age_velocity(30,1); 
Thermal_TypeTrench.vel = convert_Age_velocity(3,2); 
Thermal_TypeTrench.Type = 'McKenzie';
% 
Thermal_TypeContinent = Thermal_Type; 
Thermal_TypeContinent.Type = 'ContinentalGeotherm'; 
Thermal_TypeContinent.Moho = 600; 
Thermal_TypeContinent.Moho_d = -35; 
%==========================================================================
% General thermal/phase information (i.e. TS,TP, Phase astenosphere, Phase air)
%==========================================================================
Thermal_information = Thermal_Information;
Thermal_information.TP = 1350; 
Thermal_information.TS = 20; 
Thermal_information=Thermal_information.kappa_calculation; 
Thermal_information.Ph_Air = 0; 
Thermal_information.Ph_Ast = 5; 
%==========================================================================
% Geological object: Terrane, Trench {To Do Ridge}
% Terrane and Trench are the two major classes, Passive margin is a class
% that belongs to a terrane. 
% 
%==========================================================================
Continent1     = Terrane;
Continent2     = Terrane;
T = Trench; 
passive_margin1 = Passive_Margin; % To apply to continent 1
%==========================================================================
% Filling up the type with specific user command
%==========================================================================
passive_margin1.Direction = 'right';
passive_margin1.ph_pas_m  = phases.Ph_pas_m(1);
passive_margin1.shape     = 'trapezoidal'; 
passive_margin1.Thermal_type_O = Thermal_TypeTrench;
passive_margin1.Thermal_type_C = Thermal_TypeContinent; 
passive_margin1.Stratigraphy= continental_stratigraphy;
passive_margin1.Thermal_information = Thermal_information;
%==========================================================================
Continent1.name = 'Continent 1 ';
Continent1.Boundary = [-1000,min(Gr.y_g),0,max(Gr.y_g)];
Continent1.Stratigraphy=continental_stratigraphy; 
Continent1.Thermal_type = Thermal_TypeContinent; 
Continent1.Thermal_information = Thermal_information; 
Continent1.Passive_Margin = passive_margin1; 
%==========================================================================
Continent2.name = 'Continent 2 ';
Continent2.Boundary = [0,min(Gr.y_g),1000,max(Gr.y_g)];
Continent2.Stratigraphy=continental_stratigraphy2; 
Continent2.Thermal_type = Thermal_TypeContinent; 
Continent2.Thermal_information = Thermal_information; 
%==========================================================================
T.Boundary = 0.0;
T.Stratigraphy_Continental = continental_stratigraphy_s;
T.Stratigraphy_Oceanic = oceanic_stratigraphy;
T.theta = 50;
T.D0   = 80; 
T.L0   = 400; 
T.Decoupling_depth = -100; 
T.Thermal_type = Thermal_TypeTrench;
T.Thermal_information = Thermal_information; 
T.R = [120,T.R-T.D0]; 
T.Type_Subduction = 'Mode_1'; 
T.tk_WZ = 20.0; 
T.ph_WZ = phases.Ph_WZ(1);
T.Subducted_crust_L = [0.0,0.0,0.0,-80.0,50.0,0.0];
%                     x1   z1
T.phase_prism = phases.Ph_cont_pr(1);
T.position_end_prism = 200; 

%==========================================================================
% Filling the Terrane structure:=> The order is important as the phase will
% overlap each other, so trench must be later than two embedding terranes
% and so forth. In general trenches come after any terrane. Ridge should be
% a similar class of passive margin as it create a different thermal
% distribution for the plate. 
% Ongoing issue: How to compute the plate-model with a fixed plate
% thickness? 
% Proper continental lithosphere geotherm. 
%==========================================================================
Terranes = struct('Continent1',Continent1,'Continent2',Continent2,'T',T);
%==========================================================================
TA = cputime; 
Create_Setup(Terranes,phases,Thermal_information,A,npart,Gr,Parallel_partition);
TB = cputime;
disp('====================================================================')
disp(['The Setup is finished and took ', num2str(round((TB-TA)./60)), ' min'])
disp('====================================================================')
%========================================================================%
function [converted] = convert_Age_velocity(value,type)
        Myrs_sec = 365.25*60*60*24*1e6;
        cm_to_meter = 100; 
        if type == 1
            converted = value*Myrs_sec; 
        elseif type ==2 
            converted = value/(Myrs_sec/1e6)/cm_to_meter; 
        end
end