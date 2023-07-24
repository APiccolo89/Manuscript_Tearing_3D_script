
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


T.R        = 40;   %curvature radius
T.theta_c  = 30;   %curvature radius ingested continental crust
T.theta_dc = 20;   % additional curvature to emulate passive margin (optional)
T.theta    = 70;   % curvature slab
T.tk_WZ    = 0;   % thickness of the weak zone
Area       = 300.*80; 
L          = (300.*2)./(1+sin(T.theta*pi/180).^2);
T.L0       = L;  % length of the slab from the bottom of the lithosphere
T.D0       = 80;   % Thickness of the slab
T.Cont     = 30;   % Thickness continental crust
T.C  = [0.0 -T.D0-T.R];
T.r  = [T.R T.R+T.D0];
T.r_WZ = [T.r(2), T.r(2)+T.tk_WZ];
T.D_WZ = -100;
T.Depth_C = -40;
T.Length = 500;
T.Type = 'Mode_1'; % 'Ribe_Mode'
T.CCS.phases = [13,1,2];
T.CCS.Stratigraphy=[0,-12,-15,-30];
T.Temperature = 'McKenzie';
T.vl          = 1.0; 
T.Decoupling = -100; 



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


Buffer         = Terrane;
Continent1     = Terrane;
Oceanic_Plate  = Terrane;
Continent2     = Terrane;


% Organize
%% Buffer Terranes: Terranes at the left most area of the numerical domains:
%  it is the column of material used as reference for the lithostatic
%  calculation. And, it is a low viscosity terranes that allow convergence
%  velocity;
%==========================================================================
Buffer.order        = 1;
Buffer.Type         = 'Buffer';
Buffer.x_lim        = [-1500.0,-1300.0];
Buffer.Phases       = [phases.Ph_sed_oc(1),phases.Ph_OC(1),phases.Ph_OLt(1)];
Buffer.Stratigraphy = [0.0,-2.0,-7.0,-80.0];
Buffer.Age          = 5.0;

%% Continental Terranes 1
%Terranes.Continent1
Continent1.order = 2;
Continent1.Type  = 'Continent';
Continent1.x_lim = [-1000.0,0.0];
Continent1.Phases = [phases.Ph_UC(1),phases.Ph_LC(1),phases.Ph_Clt(1)];
Continent1.Stratigraphy = [0.0,-15.0,-30.0,-100.0];
Continent1.Age          = 100.0;
Continent1.Passive_Margin = {[],'right'};
Continent1.Passive_Margin_Age = [100.0,40.0]; 
Continent1.Passive_Margin_phase = phases.Ph_pas_m;
%% Continent Terranes 2
Continent2.order = 3;
Continent2.Type  = 'Continent';
Continent2.x_lim = [0.0,1000];
Continent2.Phases = [phases.Ph_UC2(1),phases.Ph_LC2(1),phases.Ph_Clt2(1)];
Continent2.Stratigraphy = [0.0,-15.0,-30.0,-100.0];
Continent2.Age          = 100.0;
Continent2.Accretion_prism = 'Prism';
Continent2.Trench_properties = T;
Continent2.prism_phase      = phases.Ph_cont_pr;
%% Oceanic plate
Oceanic_Plate.order = 4;
Oceanic_Plate.Type  = 'Ocean';
Oceanic_Plate.x_lim = [0.0,0.0];
Oceanic_Plate.Phases = [phases.Ph_OLt(1)];
Oceanic_Plate.Stratigraphy = [0.0,-80.0];
Oceanic_Plate.Age          = 40;

Oceanic_Plate.Trench        = 'Subduction';
Oceanic_Plate.Trench_properties = T;

Terranes = struct('Buffer',Buffer,'Continent1',Continent1,'Continent2',Continent2,'Oceanic_Plate',Oceanic_Plate);
%% Generic information numerical domain:
Gen.T_P = 1350;
Gen.T_S = 20;
Gen.Ph_Air   = phases.Ph_Ar(1);
Gen.Ph_UM    = phases.Ph_UM(1);
Gen.WZ       = phases.Ph_WZ(1);
Gen.PrismPh = phases.Ph_cont_pr(1);
Gen.AvTS    = 800;
TA = cputime;
Create_Setup(Terranes,phases,Gen,A,npart,Gr,Parallel_partition);
TB = cputime;
disp('====================================================================')
disp(['The Setup is finished and took ', num2str(round((TB-TA)./60)), ' min'])
disp('====================================================================')

