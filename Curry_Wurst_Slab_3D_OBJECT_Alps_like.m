
clear all;
close all;
addpath matlab   %Here -> Folder2LaMEM/matlab => all the function that handle matlab files are there
addpath(genpath('3D_numerical_setup'))
npart = [3,3,3];
% See model setup in Paraview 1-YES; 0-NO
Paraview_output        = 1;
% Output parallel files for LaMEM, using a processor distribution file (msetup = parallel)
LaMEM_Parallel_output  = 1;
% = = = = = = = = = = = = Script initial setup = = = = = = = = = = = = = =
Parallel_partition     = 'ProcessorPartitioning_128cpu_4.4.8.bin';
[A,Gr] = Parse_LaMEM_bin(Parallel_partition,Paraview_output,LaMEM_Parallel_output,npart);
%========================================================================== 
% Phase Lists
% Phase structure containing the phase ID number and reference density
%==========================================================================
phases.Ph_Ar  = [0,10]  ;% Air
phases.Ph_UC  = [1,2700]  ;% Upper Crust
phases.Ph_LC  = [2,2800]  ;% Lower Crust
phases.Ph_Clt = [3,3300] ;% Continental Lithosphere
phases.Ph_Clt2 = [4,3300] ;
phases.Ph_UM  = [5,3300]  ;% Upper Mantle
phases.Ph_OLt = [6,3300] ;% Ocean Lithosphere
phases.Ph_OC  = [7,3300]   ;%place holder
phases.Ph_sed_oc = [8,2700];
phases.Ph_WZ  = [9,3300]  ;% Weak Zone
phases.Ph_LC2  = [10,2700]  ;
phases.Ph_UC2  = [11,2800] ;
phases.Ph_cont_pr = [12,2700];
phases.Ph_pas_m = [13,2700];
phases.Ph_OLt2 = [14,3300] ;% Ocean Lithosphere
phases.Ph_OC2  = [15,3300]   ;%place holder
phases.Ph_cont_pr2  = [16,2800]   ;%place holder



%==========================================================================
% Stratigraphic structure list 
% Stratigraphies types ~ Structure containing the number of phases and
% layer to introduce in a specific geological unit
%==========================================================================
continental_stratigraphy.phases = [phases.Ph_UC(1),phases.Ph_LC(1),phases.Ph_Clt(1)];
continental_stratigraphy.Tk     = [0.0,-15.0,-30.0,-100.0];
oceanic_stratigraphy.phases = [phases.Ph_OC(1),phases.Ph_OLt(1)];
oceanic_stratigraphy.Tk     = [0.0,-7.0, -100.0];
oceanic_stratigraphy2.phases = [phases.Ph_OC2(1),phases.Ph_OLt2(1)];
oceanic_stratigraphy2.Tk     = [0.0,-7.0, -100.0];
continental_stratigraphy2.phases = [phases.Ph_UC2(1),phases.Ph_LC2(1),phases.Ph_Clt2(1)];
continental_stratigraphy2.Tk     = [0.0,-15.0,-30.0,-100.0];
continental_stratigraphy_s.phases =[phases.Ph_pas_m(1),phases.Ph_UC(1),phases.Ph_LC(1),phases.Ph_Clt(1)];
continental_stratigraphy_s.Tk     = [0.0,-12.0,-15.0,-30.0,-100.0];

% Boundaries
Lx = max(Gr.x_g)-min(Gr.x_g);
Wy = max(Gr.y_g)-min(Gr.y_g); 
BG = BoundsT;
BG.c = [0.0,0.0]; 
BG.W = Wy; 
BG.L = Lx; 
BG = BG.compute_coordinate_boundary_composite; 
%North Continent 

%South Continent 
C2 = BoundsT; 
C2.c = [0.0,-500];
C2.W = 1000.0;
C2.L = 1000.0;
C2=C2.Create_arc_circumference_margin(1000,'B',1200);

% = North Continent with composite D Boundary
C1 = BoundsT;
C1.c = [0.0, Wy/4];
C1.W = Wy/2;
C1.L = Lx;
C1=C1.compute_coordinate_boundary_composite;
% Trench => ~ night mare, but let's try
% Rift background oceanic terranes


%==========================================================================
% Thermal types for the specific setup
%==========================================================================

Thermal_TypeTrench = Thermal_Type;
Thermal_TypeTrench.Age = convert_Age_velocity(30,1); 
Thermal_TypeTrench.vel = {[convert_Age_velocity(3,2),convert_Age_velocity(10,2)],'linear'}; 
Thermal_TypeTrench.Type = 'McKenzie';
%
Thermal_TypeOcean = Thermal_Type;
Thermal_TypeOcean.Age = convert_Age_velocity(30,1); 
Thermal_TypeOcean.Type = 'HalfSpace';
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
Ocean_BG       = Terrane; 
Continent1     = Terrane;
Continent2     = Terrane;
T = Trench; 
T2 = Trench;
T3 = Trench;
passive_margin1 = Passive_Margin; % To apply to continent 1
passive_margin2 = Passive_Margin; % To apply to continent 1
%==========================================================================
% Filling up the type with specific user command
%==========================================================================
passive_margin1.Direction = 'right';
passive_margin1.ph_pas_m  = phases.Ph_pas_m(1);
passive_margin1.shape     = 'rectangular'; 
passive_margin1.Thermal_type_O = Thermal_TypeTrench;
passive_margin1.Thermal_type_C = Thermal_TypeContinent; 
passive_margin1.Stratigraphy= continental_stratigraphy;
passive_margin1.Thermal_information = Thermal_information;
passive_margin1.Boundary_terrane_list = {'D'}; 
% Passive margin type 2
passive_margin2.Direction = 'right';
passive_margin2.ph_pas_m  = phases.Ph_pas_m(1);
passive_margin2.shape     = 'trapezoidal'; 
passive_margin2.Thermal_type_O = Thermal_TypeTrench;
passive_margin2.Thermal_type_C = Thermal_TypeContinent; 
passive_margin2.Stratigraphy= continental_stratigraphy2;
passive_margin2.Thermal_information = Thermal_information;
passive_margin2.Boundary_terrane_list = {'B'}; 
%==========================================================================
Continent1.name = 'Continent_Main';
Continent1.Boundary = C1;
Continent1.Stratigraphy=continental_stratigraphy; 
Continent1.Thermal_type = Thermal_TypeContinent; 
Continent1.Thermal_information = Thermal_information; 
Continent1.Passive_Margin = passive_margin1; 
%==========================================================================
Continent2.name = 'Indenter';
Continent2.Boundary = C2;
Continent2.Stratigraphy=continental_stratigraphy2; 
Continent2.Thermal_type = Thermal_TypeContinent; 
Continent2.Thermal_information = Thermal_information; 
Continent2.Passive_Margin = passive_margin2; 
%==========================================================================
Ocean_BG.name = 'Ocean Background';
Ocean_BG.Boundary = BG; 
Ocean_BG.Stratigraphy = oceanic_stratigraphy2;
Ocean_BG.Thermal_information = Thermal_information;
Ocean_BG.Thermal_type = Thermal_TypeOcean;
%==========================================================================
T.Boundary = C2;
T.Boundaries_list = {'B'}; 
T.Stratigraphy_Continental = continental_stratigraphy_s;
T.Stratigraphy_Oceanic = oceanic_stratigraphy;
T.theta = {[-90,-90],'none'};
T.D0   = 100; 
T.L0   = 400; 
T.Decoupling_depth = -100; 
T.Thermal_type = Thermal_TypeTrench;
T.Thermal_information = Thermal_information; 
T.R = [120-T.D0, 120]; 
T.Type_Subduction = 'Mode_1'; 
T.tk_WZ = 20.0; 
T.ph_WZ = phases.Ph_WZ(1);
T.Subducted_crust_L = [0.0,0.0,0.0,-80.0,50.0,0.0];
%                     x1   z1
T.phase_prism = {phases.Ph_cont_pr(1),phases.Ph_cont_pr2(1)};
T.position_end_prism = 80; 
T.length_continent = {[150,20],'linear'};
T.prism_depth = -70; 

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
Terranes = struct('Ocean_BG',Ocean_BG,'Continent1',Continent1,'Continent2',Continent2,'Trench',T);%,'T',T);
%==========================================================================
TA = cputime; 
Create_Setup(Terranes,phases,Thermal_information,A,npart,Gr,Parallel_partition,1);
TB = cputime;
disp('====================================================================')
disp(['The Setup is finished and took ', num2str(round((TB-TA)./60)), ' min'])
disp('====================================================================')
%========================================================================%















