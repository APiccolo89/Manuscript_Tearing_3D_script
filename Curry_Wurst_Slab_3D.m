
clear all;
close all; 
addpath matlab   %Here -> Folder2LaMEM/matlab => all the function that handle matlab files are there


npart = [3,3,3];
% See model setup in Paraview 1-YES; 0-NO
Paraview_output        = 1;
% Output parallel files for LaMEM, using a processor distribution file (msetup = parallel)
LaMEM_Parallel_output  = 1;
Parallel_partition     = 'ProcessorPartitioning_128cpu_4.4.8.bin';
[A,Gr] = Parse_LaMEM_bin(Parallel_partition,Paraview_output,LaMEM_Parallel_output,npart,1);
Radius_terranes = 700; 
arc_length = 1200; 

%phase database, first element phase.ph_Air(1) = phase number, 2nd element
%density (useful and necessary for the isostasy)

phases.Ph_Ar  = [0,10]  ;% Air 
phases.Ph_UC  = [1,2700]  ;% Upper Crust 
phases.Ph_LC  = [2,2800]  ;% Lower Crust 
phases.Ph_Clt = [3,3300] ;% Continental Lithosphere
phases.Ph_Clt2 = [4,3300] ; 
phases.Ph_WZ  = [8,3300]  ;% Weak Zone 
phases.Ph_OLt = [6,3367.957866] ;% Ocean Lithosphere 
phases.Ph_UM  = [5,3300]  ;% Upper Mantle 
phases.Ph_OC  = [7,3367.957866]   ;%place holder
phases.Ph_LC2  = [9,2700]  ;
phases.Ph_UC2  = [10,2800] ; 
phases.Ph_sed_oc = [11,2680];
phases.Ph_cont_pr = [12,2680];
phases.Ph_pas_m = [13,2680];

phases.Ph_OLt2 = [19,3367.957866] ;% Ocean Lithosphere 
phases.Ph_OC2  = [20,3367.957866]   ;%place holder
% Structure concerning slab 
T.R        = 40;   %curvature radius 
T.theta_c  = 90;   %curvature radius ingested continental crust
T.theta    = 90;   % curvature slab
T.tk_WZ    = 20;   % thickness of the weak zone
T.L0       = 300;  % length of the slab from the bottom of the lithosphere
T.D0       = 80;   % Thickness of the slab
T.C  = [0.0 -T.D0-T.R];  %center of curvature 
T.r  = [T.R T.R+T.D0]; %radius of curvature
T.r_WZ = [T.r(2), T.r(2)+T.tk_WZ]; %weak zone data
T.D_WZ = -100; %ultimate weakzone depth 
T.Type = 'Mode1'; % 'Ribe_Mode' "Mode1" create the slab following radius of curvature and tagent to the curvature, Ribe_Mode uses equation of Ribe et al 2010
T.Continent_phases = [10,9];
T.Stratigraphy     = [-15,-30];
T.Radius_terranes = Radius_terranes;
T.arc_length      = arc_length;
% List of terranes 
Continent1                = Terrane; %continent 1 
Oceanic_Plate             = Terrane; %oceanic plate 
Continent2                = Terrane; %continent 2
Oceanic_Plate_Background  = Terrane; %oceanic plate 



% Organize 
%% Buffer Terranes: Terranes at the left most area of the numerical domains: 
%  it is the column of material used as reference for the lithostatic
%  calculation. And, it is a low viscosity terranes that allow convergence
%  velocity; 
%==========================================================================


%% Continental Terranes 1 
[xa,xb,xc,yc]=find_terrane_size(Radius_terranes,arc_length);
disp(['x coordinates', num2str(xa)])
disp(['center coordinates', num2str(xc),'' ,num2str(yc)])
f_c2_handle = @(x) circumference_margin(xa,xb,0.0,0.0,Radius_terranes,x);%parabolic_margin(-750.0,0.0,70.0,x)
f_continent_handle = @(x) continental_theta_angle(x,0.0,arc_length,T.theta_c,0);
f_theta_handle = @(x) continental_theta_angle(x,0.0,arc_length,T.theta,60);
T.function_continent = f_continent_handle; 
%T.function_bending_plate = f_theta_handle; 

%% Oceanic plate 
Oceanic_Plate_Background.order = 1; 
Oceanic_Plate_Background.Type  = 'Ocean';
Oceanic_Plate_Background.x_lim = [-1500.0, 1500.0]; 
Oceanic_Plate_Background.y_lim  = [-1500.0, 0.0]; 
Oceanic_Plate_Background.Phases = [phases.Ph_OC2(1),phases.Ph_OLt2(1)];
Oceanic_Plate_Background.Stratigraphy = [0.0,-7.0,-80.0];
Oceanic_Plate_Background.Age          = 30.0; 
%Oceanic_Plate_Background.Trench        = 'Subduction'; 
%Oceanic_Plate_Background.Trench_properties = T;  

%Terranes.Continent1
Continent1.order = 2; 
Continent1.Type  = 'Continent';
Continent1.x_lim = [-1500.0,1500.0];
Continent1.y_lim = [0.0,1500.0];
Continent1.Phases = [phases.Ph_UC(1),phases.Ph_LC(1),phases.Ph_Clt(1)];
Continent1.Stratigraphy = [0.0,-15.0,-30.0,-80.0];
Continent1.Age          = 100.0; 
Continent1.Passive_Margin = {'segment1','segment2'};
Continent1.Passive_Margin_segment = {'along_x',[-1500.0,xa],[xb,1500.0]}; 
Continent1.Passive_Margin_phase = phases.Ph_pas_m;
Continent2.Accretion_prism = 'Prism'; 
Continent2.Trench_properties = {T,f_c2_handle,[xa,xb]}; 
Continent2.prism_phase      = phases.Ph_cont_pr;
%% Continent Terranes 2 
Continent2.order = 3; 
Continent2.Type  = 'Continent';
Continent2.x_lim = [xa,xb];
%f_c2_handle = @ (x) linear_margin(Continent2.x_lim(1),30,0.0,x);
Continent2.y_lim = {-500,f_c2_handle};
Continent2.Phases = [phases.Ph_UC2(1),phases.Ph_LC2(1),phases.Ph_Clt2(1)];
Continent2.Stratigraphy = [0.0,-15.0,-30.0,-80];
Continent2.Age          = 80.0; 

%% Oceanic plate 
 Oceanic_Plate.order = 4; 
 Oceanic_Plate.Type  = 'Ocean';
 Oceanic_Plate.x_lim = [xa,xb]; 
 Oceanic_Plate.y_lim = {f_c2_handle,f_c2_handle};
 Oceanic_Plate.Phases = [phases.Ph_OC(1),phases.Ph_OLt(1)];
 Oceanic_Plate.Stratigraphy = [0.0,-7.0,-80.0];
 Oceanic_Plate.Age          = 30.0; 
 Oceanic_Plate.Trench        = 'Subduction'; 
 Oceanic_Plate.Trench_properties = T;  

Terranes = struct('Oceanic_Plate_Background',Oceanic_Plate_Background,'Continent1',Continent1,'Continent2',Continent2,'Oceanic_Plate',Oceanic_Plate); 
%% Generic information numerical domain: 
Gen.T_P = 1375; 
Gen.T_S = 20; 
Gen.Ph_Air   = phases.Ph_Ar(1);
Gen.Ph_UM    = phases.Ph_UM(1);
Gen.WZ       = phases.Ph_WZ(1); 
Gen.PrismPh = phases.Ph_cont_pr(1);
TA = cputime;
Create_Setup(Terranes,phases,Gen,A,npart,Gr,Parallel_partition);
TB = cputime; 
disp('====================================================================')
disp(['The Setup is finished and took ', num2str(round((TB-TA)./60)), ' min'])
disp('====================================================================')
%=========================================================================%
%%
%==========================================================================
% Function for constructing the initial setup [Still tuned for 2D and
% monodirectional]
%==========================================================================

function Create_Setup(Terranes,ph,Gen,A,npart,Gr,Parallel_partition)
    
    RandomNoise             =   logical(0);
    Is64BIT                 =   logical(0);
    Phase = 0.0.*((A.Xpart));
    Temp  = 0.0.*((A.Xpart));

    % Set Mantle Phase
    Phase(:,:,:)  = Gen.Ph_Air;
    Temp(:,:,:)   = Gen.T_P;

    % Set up continental terranes: 
    % Temperature is computed only for the layer considered, while the rest
    % is assumed to be not affected by the half space cooling - just to
    % simplify a bit 
    % To Do => Convert this shit in Python
    % To Do => Create ad hoc function
    
    %=====================================================================%
    terranes_list = fieldnames(Terranes); 
    for it =1:length(terranes_list)
            t = Terranes.(terranes_list{it});
            disp(['Filling the ',terranes_list{it} ,' terranes .... '])

            a = cputime;
            [Phase,Temp] =  Set_Phase_Temperature(A,Phase,Temp,t,Gen);
            b = cputime; 
            time = b-a; 
            disp(['took ', num2str(round(time)), ' s'])
            disp('=================================')


    end
    %===========================================================================%
    % Set Air Phase
    % Adiabatic gradient 
    
    ind = Phase == 0 & A.Zpart<0.0;
    Phase(ind)  = Gen.Ph_UM;
    Temp(ind)   = Gen.T_P;
    %Temp = Temp+abs(A.Zpart).*0.3;
    ind = A.Zpart>0.0;
    Temp(ind)   = Gen.T_S;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    A.Xpart  =  permute(A.Xpart,[2 1 3]);
    A.Ypart  =  permute(A.Ypart,[2 1 3]);
    A.Zpart  =  permute(A.Zpart,[2 1 3]);


    % We can still manually change all this & include, say, a lower mantle
    A.nump_x = npart(1);
    A.nump_y = npart(2);
    A.nump_z = npart(3);
    A.Phase  = double(Phase); clear Phase
    A.Temp   = double(Temp);  clear Temp
    A.Phase  = permute(A.Phase,[2 1 3]);
    A.Temp   = permute(A.Temp, [2 1 3]);

    x = squeeze(A.Xpart(:,1,1));
    y = squeeze(A.Ypart(1,:,1));   
    z = squeeze(A.Zpart(1,1,:));   
    A.x      =  double(x(:));
    A.y      =  double(y(:));
    A.z      =  double(z(:));
    
    
   %[A,surf] = displace_phase_isostasy(ph,A,Gr,Gen); 

   % plot_initial_setup2D(A,surf); 
   
    A.RandomNoise = logical(0);

    clear Temp Phase


    % PARAVIEW VISUALIZATION
    FDSTAGWriteMatlab2VTK(A,'BINARY'); % default option

    % SAVE PARALLEL DATA (parallel)
    FDSTAGSaveMarkersParallelMatlab(A,Parallel_partition, Is64BIT);
    
end


function [Phase,Temp] =  Set_Phase_Temperature(A,Phase,Temp,Terranes,Gen)
%=========================================================================%
% Function that fills the phases. 
%=========================================================================%
% Input Parameter: 
%=========================================================================%
% A = Data Structure with the grid/particles data. 
% Phase, Temp Layer of phase or initial temperature field. 
% Terranes: object containing information of the terranes. 
% Gen     : generic data structure containing the information that are
% general (i.e. temperature of air, potential temperature, phase of prism,
% or passive margin or weak zone). 
%=========================================================================%
% Output parameter: 
%=========================================================================%
% Phase,Temp update. 
%=========================================================================%
% Structure of the function: 
% Fill layer between fixed bounds. Add subduction zone, add accretionary
% prism or passive margin and transitional terranes. Each of the generation
% of terranes is handled with specific functions. 
%=========================================================================%
trench = Terranes.Trench; 
passive_margin = Terranes.Passive_Margin;
accretion_prism = Terranes.Accretion_prism; 
A_fill_layer = cputime;
[Phase,Temp] = fill_layer(A,Terranes,Phase,Temp,Gen);
B_fill_layer = cputime;
t = B_fill_layer-A_fill_layer; 
disp(['           1. Fill the layer and initialise Temperature field...', num2str(round(t)), 's']);
if strcmp(trench,'Subduction')
    A_fill_Sub = cputime;
    [Phase,Temp] = fill_subduction(A,Terranes,Phase,Temp,Gen);
    B_fill_Sub = cputime;
    t = B_fill_Sub-A_fill_Sub; 
    disp(['        2. Fill the Subduction and initialise its Temperature field...', num2str(round(t)), 's']);

end
if strcmp(accretion_prism,'Prism')
   A_fill_Prism = cputime;
   [Phase] = generate_accretion_prism(A,Terranes,Phase);
   B_fill_Prism = cputime;
   t = B_fill_Prism-A_fill_Prism; 
   disp(['         2. Fill the Prism ...', num2str(round(t)), 's']);
end
if strcmp(passive_margin,'none') == 0
   % generate left/right passive margin 
   A_fill_Passive = cputime;

   for i=1:length(passive_margin)
        direction = Terranes.Passive_Margin{i};
        limits=Terranes.Passive_Margin_segment{i+1}; 
        [Phase,Temp] = generate_passive_margin(A,Phase,Temp,Terranes,direction,Gen,limits);
   end
     B_fill_Passive = cputime;
     t = B_fill_Passive-A_fill_Passive; 
     disp(['       3. Fill the Passive Margins...', num2str(round(t)), 's']);
   
end

end



function [Temp]= compute_temperature_profile(A,Temp,Type,Gen,Terranes,indx)
T_tk = [0.0, Terranes.Stratigraphy(end)];
k = Terranes.K./(Terranes.Cp.*3300); 
T_Age = Terranes.Age.*Terranes.secMyrsyear;
T_P   = Gen.T_P;
T_S   = Gen.T_S; 

if Type == 1
    ind_z = find(A.Zpart<T_tk(1) & A.Zpart>=T_tk(2) & indx == 1);
    erf_function = erf(A.Zpart(ind_z).*1000/2/(k*T_Age)^0.5);
    Temp(ind_z) = T_S - (T_P-T_S).*erf_function;
else
    ind = A < T_tk(1) & A>=T_tk(2);
    erf_function = erf(A(ind).*1000/2/(k*T_Age)^0.5);
    Temp(ind) = T_S - (T_P-T_S).*erf_function;
    ind = A < T_tk(2); 
    Temp(ind) = T_P;

end

end




function [Layout,angle_c] = find_slab_(A,Slab,type,Terranes)
%=========================================================================
% Short explanation, i need to yet find the best solution that works in any
% case provided. In general Ribe mode require a bit of work as such angle
% that are close to 90 are considered. However, the formulation requires
% the massive usage of tan function, that is not continuous for angle close
% to 90. So, if someone want to have a vertical slab, is always convinient
% to use the mode1. The Ribe mode is more suitable for model featuring an
% initial convergence, as it should represents the geometry of a bended
% plate under its own weight, and should give more reliable flexural stress
%=========================================================================
% Really convoluted, but it works: 
% 1) create a polygon for the slab 
% 2) find the curvature 
% 3) save the relative distannce w.t.r the top
% How it works? 
% 1) Extract the information of the slab
% 2) Select a rectangular area in the surrounding of the center of
% curvature: 
% Mode_1: 
% select 2 circumference arc that represents the top of the slab and the
% bottom of the slab. 
% The distance from the top surface represents the depth of the slab and is
% computed along the radius of the circumference (i.e., we have a point
% within the area in which the slab existing). Associated to this point we
% have its angular distance w.r.t. the center of the curvature. If the
% point belongs to the set of angular distance that we are interested to,
% we computed the distance .w.r.t. to the center (which is, surprise
% surprise is the radius) then we substract the coordinate of the top
% surface. distance information is then used to fill up the stratigraphy
% and to compute the half-space cooling model. 
% Point (x,z) => compute the distance with C (center of curvature)  
%             => compute the angle between the vertical vector stemming
%             from the center and the actual position vector of the point 
% 
%               u = [x(i,lx,lz);z(i,lx,lz)]; position VECTOR current point
%               v = [C(1);z(i,lx,lz)];       vertical vector stemming from
%               the center 
%               a = sqrt(u(1).^2+u(2).^2);   %vector magnitude position
%               b = sqrt(v(1).^2+v(2).^2);   % vector magnitude vertical
%               vector
%               d(i,lx,lz) = sqrt((u(1)-C(1))^2+(u(2)- C(2))^2); %distance
%               betwee point and center
%               arc_angleS(i,lx,lz) =
%               acos((z(i,lx,lz)-C(2))/d(i,lx,lz)).*180/pi; % angle between
%               vertical vector and actual position vector 
%               => if arc_angleS belongs to 0-theta and has a distance
%               consitent within top and bottom of the slab 
%                if (d(i,lx,lz)>=r(1) && d(i,lx,lz)<=r(2))
%                    if(arc_angleS(i,lx,lz)>=0.0 && arc_angleS(i,lx,lz)<=theta)
%
%                        Layout(i,lx,lz) = d(i,lx,lz)-r(2);
%                    end
%                end
%   I'm sure that exists a more efficient way to handle the index, but i'm
%   really dump to vectorise matlab. If you find a decent way, please, feel
%   free to violate my script and give me the solutions :) 
% Ribe mode. I did not find a better name for this way to find a slab. The
% ribe mode is relying on the initial flexure of slab provided by Prof.Dr.
% Ribe in his paper 2010. Long story short, the derivation relises on the
% fiber stress. What is the fiber midplane ? The theoretical background is
% the thin sheet layer theory or boundary layer theory. Within an
% elastic/viscous sheet exist a surface in which the stresses are equal to
% 0.0 the mid plane. If you bend a plate, you have extensional stresses at
% the top, and compressive stress at the bottom of the bulked plate. 
% The algoritm is to find the distance from this mid surface using his
% equation, the main issue is related to the resolution and to the
% definition of the area where to look. One way is to allow the user to
% tune until he find the best solution, the other is provided a more robust
% search algorithm by default. This mode allow to not generate unrealistic
% stress during the first timestep, and to prevents odds curvature. The
% best way to have a subduction is the spontaneous, however, it is a tedius
% task and time to time we need to save computational time. 
%
%
% find the area belonging to the curvature slab: 
Xvec = A.Xpart(1,:,1);
Yvec = A.Ypart(:,1,1);
Zvec = A.Zpart(1,1,:);
C     = Slab.C;
D0    = abs(Slab.D0);


%% New Approach to improve
fhandle = Terranes.y_lim{2};
fhandle_c = Slab.function_continent;
fhandle_arc_length = @(x,y) arc_length_x(x,y,Terranes.x_lim(1),0.0,Slab.Radius_terranes);
if isfield(Slab,'function_bending_plate')
    f_handle_theta = Slab.function_bending_plate; 
end

theta = Slab.theta;
theta_c = Slab.theta_c; 
ang_  = (theta)*pi/180;
ang_2 = (theta+90)*pi/180;

r     = Slab.r;
L0    = Slab.L0 - Slab.R;
sl    = 1               ;
arc_angleS = []; 

ind = (A.Xpart(:)>=Terranes.x_lim(1)) & (A.Xpart(:)<=Terranes.x_lim(2)) & (A.Ypart(:)>=(fhandle(A.Xpart(:)))) & (A.Ypart(:)<=(fhandle(A.Xpart(:)))+2.*L0) & A.Zpart(:)>=-1.5*L0 & A.Zpart(:)<=1.0; 
index = find(ind==1); 
Layout = ind*nan;
d      = Layout; 
l_ind = length(ind); 
 
 A_time = cputime; 
 x     = A.Xpart(:);
 z     = A.Zpart(:);
 y     = A.Ypart(:);
 x     = x(ind==1);
 y     = y(ind==1);
 z     = z(ind==1); 
 % Vector version 
 arc_angleS = x.*NaN; 
 continent = arc_angleS;
 s = fhandle_arc_length(x,fhandle(x)); 

 Layout_buf = x.*NaN; 
 u = zeros(2,length(x));
 v = u; 
 u(1,:) = y; 
 u(2,:) = z;
 v(1,:) = fhandle(x);
 v(2,:) = z; 
 d = sqrt((u(1,:)-fhandle(x)').^2+(u(2,:)- Slab.C(2)).^2);
 arc_angleS = acos((z-Slab.C(2))./d').*180/pi;
 if isfield(Slab,'function_bending_plate')
     c_slab = d>=r(1) & d<=(r(2)) & arc_angleS'<=f_handle_theta(s)' & arc_angleS'>=0;
     continent(arc_angleS' <= fhandle_c(s)' & arc_angleS'<f_handle_theta(s)')=1.0;
 else
     continent(arc_angleS'<fhandle_c(s)')=1.0;
     c_slab = d>=r(1) & d<=(r(2)) & arc_angleS'<=theta & arc_angleS'>=0;
 end
 Layout_buf(c_slab==1)= d(c_slab==1)-r(2); 
 Layout_buf(d<r(1) & d>r(2)) = NaN; 
 Layout_buf(arc_angleS>theta & arc_angleS<0)=NaN; 
  
 x_buf = squeeze(A.Xpart(1,:,1));
 x_buf = x_buf(x_buf>=Terranes.x_lim(1) & x_buf<=Terranes.x_lim(2)); 
 if strcmp(type,'Weak')
     r(1) = r(2);
     r(2) = r(2)+Slab.tk_WZ;
     if Slab.theta == 90
         sl = 0;
     end
     %=============

 end


 
 for  i=1:length(x_buf) 
 sections = x==x_buf(i);
  if isfield(Slab,'function_bending_plate')
    x_s = (x_buf(i));
    y_s = fhandle(x_buf(i));
    s_s = fhandle_arc_length(x_s,y_s); 
    ang_  = (f_handle_theta(s_s))*pi/180;
    ang_2 = (f_handle_theta(s_s)+90)*pi/180;
  end


 yy = y(sections==1);
 zz = z(sections==1); 

 C(1)     = x_buf(i);
 C(2)     = fhandle(x_buf(i));
 C(3)     = Slab.C(2);
 p1(1)    = C(2)+r(1)*sin(ang_);
 p1(2)    = C(3)+r(1)*cos(ang_);
 p2(1) = C(2)+r(2)*sin(ang_);
 p2(2) = C(3)+r(2)*cos(ang_);
 p3(1)    = p1(1)+L0*sin(ang_2);
 p3(2)    = p1(2)+L0*cos(ang_2);
 p4(1)    = p2(1)+L0*sin(ang_2);
 p4(2)     = p2(2)+L0*cos(ang_2);
% 
 Py(1) = p1(1);
 Py(2) = p2(1);
 Py(3) = p4(1);
 Py(4) = p3(1);
 %Pz = [p1(2,:) p2(2,:) p4(2,:) p3(2,:)];
 Pz(1) = p1(2);
 Pz(2) = p2(2);
 Pz(3) = p4(2);
 Pz(4) = p3(2);
%       
  [in,~] = inpolygon(yy,zz,Py,Pz);
 P = zeros(2,length(zz)); 
 P(1,:) = yy;
 P(2,:) = zz; 
 Layout_sections = Layout_buf(sections==1);
 Layout_sections(in==1) = -find_distance_linear(P(:,in==1),p2,p4);
 Layout_buf(sections==1) = Layout_sections; 
 end

        
 B_time = cputime; 
 disp(['Time Loop = ', num2str(B_time-A_time)]);
angle_c = Layout*0.0; 
 Layout(index)=Layout_buf; 
angle_c(index)= continent; 
Layout = reshape(Layout,size(A.Xpart));
angle_c = reshape(angle_c,size(A.Xpart)); 
 bla = 0; 


if strcmp(Slab.Type,'Ribe_Mode')
    Layout = ones(size(x))*nan;
    ind_z = find(Zvec >= C(2)-Slab.L0 & Zvec<=0.0);
    ind_x = find(Xvec>=C(1) & Xvec <= C(1)+L0+L0*0.3);
    weak_flag = 0; 
    
    if strcmp(type,'Weak')
        weak_flag = 1; 
    end
    % Ribe equation can create artifacts if certain market that does not
    % belong to the slab are taken into account. Since most of the
    % algorithm lies on the distance estimation from a point and a generic
    % curve, it is not possible to predict before hand which marker is not
    % belonging to the slab. A quick hack is to select the markers that
    % feature a z position passing through the mid surface perpendicolarly.
    % 
    z_mid_f =-D0/2+tan(theta*pi/180).*L0 ;
    P_B  = [L0-sin(theta*pi/180), z_mid_f-cos(theta*pi/180)*D0/2];
    P_T  = [L0+sin(theta*pi/180), z_mid_f+cos(theta*pi/180)*D0/2];
    line_m = (P_T(2)-P_B(2))./(P_T(1)-P_B(1)); %


    for ix = 1:length(ind_x)
        for iz = 1:length(ind_z)
            lx = ind_x(ix);
            lz = ind_z(iz);
            [z_mid,dz,z_top,z_bottom] = compute_ribe_angles(A.Xpart(:,lx,lz),C,theta,Slab.L0,D0,weak_flag,Slab.tk_WZ);
            if weak_flag == 0 
                z1 =  z_bottom-0.6*D0;
                z2 =  z_top   + 0.6*D0;
            else
                z1 = z_bottom-D0; 
                z2 = z_top+ D0; 
            end
            if A.Zpart(:,lx,lz)>= z1 & A.Zpart(1,lx,lz)<=z2 & A.Zpart(1,lx,lz) >= line_m.*(A.Xpart(1,lx,lz)-P_B(1))+P_B(2)
            
                x_t = (z_mid-A.Zpart(:,lx,lz))./sqrt(1+dz.^2);
                Layout(:,lx,lz)        =x_t;
                Layout(:,lx,lz)        = -(x_t+D0/2);
               
                if (Layout(:,lx,lz)>0 | Layout(:,lx,lz)<-D0) & weak_flag==0  
                    Layout(:,lx,lz)        = nan;
                elseif (Layout(:,lx,lz)<0 |Layout(:,lx,lz)>Slab.tk_WZ) & weak_flag==1
                    Layout(:,lx,lz)        = nan;
                end
     

            end
        end
    end

end

end
function [z,dz,ztop,zbot] = compute_ribe_angles(xp,C,teta_0,L,D0,weak_flag,Tk)
%=========================================================================%
% TO DO A complete explanation 
%=========================================================================%
teta_0 = teta_0*pi/180;

x = xp-C(1);
if xp<C(1) | xp>C(1)+L+L*0.3
    teta_s = nan;
    z      = nan;
    dz     = nan;
    z_bot  = nan;
    z_top  = nan;
else
    teta_s               = teta_0.*x.^2.*((3.*L-2.*x))./(L^3);
    z                    = -D0/2-tan(teta_s).*x;
    dz                   = (6.*teta_0.*x.^2.*(x-L).*sec(teta_s).^2)./(L^3)-tan(teta_s);
    zbot                 = z-D0/2.*cos(teta_s);
    ztop                 = z+D0/2.*cos(teta_s);
    if weak_flag == 1
        zbot                 = z+D0/2.*cos(teta_s);
        ztop                 = z+(D0/2+Tk).*cos(teta_s);
    end

end
end

function [d] = find_distance_linear(P,p1,p2)

% Ok, wiki saves my day because I was too lazy: 
A  = abs((p2(1)-p1(1)).*(p1(2)-P(2,:))-(p1(1)-P(1,:)).*(p2(2)-p1(2)));
B  = (p2(1)-p1(1)).^2+(p2(2)-p1(2)).^2;
d  = A./sqrt(B);

end



function [A,surf] = displace_phase_isostasy(ph,A,Gr,Gen) 
% =========================================================================
% Function that compute the relative displacement of the phase if we
% consider them in Isostatic equilibrium. 
%==========================================================================
% Theoretical background: I choose the bottom of the model as reference
% depth. This reference depth represents the compensation depth (i.e.,
% where the lithostatic pressure should be equal and do not show any
% horizontal gradient). Then, using the leftest colum i compute the
% refernce density of this column. The density is computed using a DataBase
% of densities (in the later stage I introduce a simplified
% pressure-temperature dependence on the density) associated to the phase,
% and the average density is assumed to be proportional to the fraction of
% each phase along z direction. 
% Then the topography is computed w.r.t. density reference colum. After
% this simple computation, I use a moving average to smooth the topography
% and use the topography as displacement vector for each marker along z. 
% Then, as last step, I interpolate the marker population and temperature
% column wise using the nearest to the original grid, and fill the the nan
% value emerging with the bottom phase or the air phase. 
%==========================================================================
% Input argument:
% A = Structure of particles 
% ph = Phase data base 
% Gr = numerical grid information
% Output Argument: 
% A = update structure of particles 
%==========================================================================

[topo,~] = compute_topography_(A,ph);

Ph = squeeze(A.Phase(:,1,:));
T  = squeeze(A.Temp(:,1,:));
Z  = squeeze(A.Zpart(:,1,:));
x  = squeeze(A.Xpart(:,1,:));
ilx = length(x(:,1));
ily = length(squeeze(A.Ypart(1,:,1)));
Ph2 = 0.*Ph;
T2  = 0.*Ph; 
topo_M = movmean(topo,500); 

for i = 1:ilx
    z = Z(i,:)+topo_M(i);
    Ph2(i,:)=interp1(z,Ph(i,:),Z(i,:),'nearest');
    T2(i,:)=interp1(z,T(i,:),Z(i,:),'nearest'); 
    z = [];
end
 
Ph2(isnan(Ph2)& Z >0.0)=ph.Ph_Ar(1);
Ph2(isnan(Ph2)& Z <0.0)=ph.Ph_UM(1);
T2(isnan(T2) & Z >0.0)=Gen.T_S; 
T2(isnan(T2)& Z <0.0)=Gen.T_P; 

for iy = 1:ily 
    A.Phase(:,iy,:) = Ph2; 
    A.Temp(:,iy,:) = T2; 
end


[s]=save_topography(A,topo_M,Gr); 
surf(1,:)=x(:,1);
surf(2,:)=topo_M; 
disp(s); 

end

function [topo,rho] = compute_topography_(A,ph)

 x = (squeeze(A.Xpart(:,1,1)));
 rho_x = 0.*x; 
 topo  = 0.*x;

 for i = 1:length(x)
     Ph = A.Phase(i,:,:);
     Ph = Ph(Ph>0);
     l_ph = unique(Ph);

     rho = 0;
     for ip = 1:length(l_ph)

         Per_ph = length(Ph(Ph==l_ph(ip)))/length(Ph);

         rho    =rho+ Per_ph*phase_density(ph,l_ph(ip));

     end
     if i == 1
         rho_0 = rho;
     end
     topo(i) = 1000.*(rho_0-rho)./rho;
     rho_x(i)  = (rho);
     Ph = [];
     l_ph = [];
     rho = [];
 end
end

function [rho] = phase_density(ph,Ph)
field  = fieldnames(ph); 

for is = 1:numel(field)
    A=ph.(field{is});
    if A(1)==Ph
        rho = A(2);
        break
    end
end
end

function [Phase,Temp] = fill_layer(A,Terranes,Phase,Temp,Gen)
% Select the layer and portion of the model:
if strcmp(Terranes.y_lim,'none')
    indx = A.Xpart >= Terranes.x_lim(1) & A.Xpart < Terranes.x_lim(2);
else
    
        if iscell(Terranes.y_lim)
            c = {Terranes.y_lim{1},isa(Terranes.y_lim{1},'function_handle')};
            d = {Terranes.y_lim{2},isa(Terranes.y_lim{2},'function_handle')};
        else
            c = Terranes.y_lim(1);
            d = Terranes.y_lim(2); 
        end
        if iscell(Terranes.x_lim)
            a = {Terranes.x_lim{1},isa(Terranes.x_lim{1},'function_handle')};
            b = {Terranes.x_lim{2},isa(Terranes.x_lim{2},'function_handle')};
        else
            a = Terranes.x_lim(1);
            b = Terranes.x_lim(2);
        end
        %(var,less_higher,A,axis)
        indx  = 0.0*A.Zpart;
        x_a = check_belonging_layer(a,1,A,1);
        x_b = check_belonging_layer(b,0,A,1);
        y_a = check_belonging_layer(c,1,A,2);
        y_b = check_belonging_layer(d,0,A,2); 
        indx(x_a==1 & x_b ==1 & y_a == 1 & y_b ==1) = 1.0; 
    
    
end
% Compute the thermal profile
[Temp] = compute_temperature_profile(A,Temp,1,Gen,Terranes,indx);
% Fill the phase stratigraphy: 
[Phase] = fill_stratigraphy(A.Zpart,Phase,Terranes,indx);
end

function [Phase,Temp] = fill_subduction(A,Terranes,Phase,Temp,Gen)
ind_x=[];

[Layout,angle_c] = find_slab_(A,Terranes.Trench_properties,'Slab',Terranes);
% Set the temperature
[Temp] = compute_temperature_profile(Layout,Temp,-1,Gen,Terranes,ind_x);%(A,Temp,Type,T_tk,Gen,Terranes,indx)
% Correct Temperaure and Phase
id1 = min(A.Xpart(~isnan(Layout)));
id2 = max(A.Xpart(~isnan(Layout)));
id3 = min(A.Zpart(~isnan(Layout)));
ind_x1 = find(squeeze(A.Xpart(1,:,1))>=id1,1);
ind_x2 = find(squeeze(A.Xpart(1,:,1))>=id2,1);
ind_z1 = find(squeeze(A.Zpart(1,1,:))>=id3,1);
for i= ind_x1:ind_x2
    ind_L = find((Layout(1,i,:)==min(Layout(1,i,:))),1);
    Phase(:,i,ind_z1:ind_L) = Gen.Ph_Air;
    Temp(:,i,ind_z1:ind_L) = Gen.T_P;
end
[Phase] = fill_stratigraphy(Layout,Phase,Terranes,ind_x);
ind = angle_c ==1 & Layout >= -15.0 & Layout <=0.0; 
Phase(ind) = Terranes.Trench_properties.Continent_phases(1); 
ind2 = angle_c ==1 & Layout >= -30.0 & Layout <=-15.0; 
Phase(ind2) =Terranes.Trench_properties.Continent_phases(2); 


[Layout,~] = find_slab_(A,Terranes.Trench_properties,'Weak',Terranes);
if strcmp( Terranes.Trench_properties,'Mode_1')
    ind =(Layout<=0.0 & A.Zpart >= Terranes.Trench_properties.D_WZ & Phase ~=Gen.PrismPh);
    Phase(ind) = Gen.WZ;
else
    ind =(Layout>0.0 & A.Zpart >= Terranes.Trench_properties.D_WZ & Phase ~=Gen.PrismPh & A.Zpart<0.0);
    Phase(ind) = Gen.WZ;
end
end


function [Phase] = fill_stratigraphy(Z,Phase,Terranes,indx)
T_Tk = [0.0, Terranes.Stratigraphy(end)];
T = T_Tk(1); 
for i=1:length(Terranes.Phases)
    if (i == length(Terranes.Phases))
        B=T_Tk(2);
    else
        B = Terranes.Stratigraphy(i+1);
    end
    if ~isempty(indx)
        ind = Z < T & Z>=B & indx>0;
    else
        ind = Z < T & Z>=B;
    end
    Phase(ind) = Terranes.Phases(i);
    T=B;
    ind = [];
end

end


function  [Phase] = generate_accretion_prism(A,Terranes,Phase)

T = Terranes.Trench_properties{1};
f_c = Terranes.Trench_properties{2};
xlim = Terranes.Trench_properties{3};
C_z = T.C(2);
ind_x = find(A.Xpart(1,:,1)>=xlim(1),1);
ind_x2 = find(A.Xpart(1,:,1)<=xlim(2));
ind_x2 = ind_x2(end);
x =squeeze(A.Xpart(1,:,1));
y = squeeze(A.Ypart(:,1,:));
z = squeeze(A.Zpart(:,1,:));

for ix = ind_x:ind_x2
    y_C= f_c(x(ix));
    d_p = [y_C+Terranes.position_end_prism 0.0];
    s   = (d_p(2)-C_z)./(d_p(1)-y_C);
    Phases = Terranes.Phases;
    Ph = squeeze(Phase(:,ix,:));
    for i = 1:length(Phases)
        ind2 = z>s.*(y-y_C)+C_z  & y>=y_C & z<=1.0 & z>=Terranes.Stratigraphy(end);
        Ph(ind2==1) = Terranes.prism_phase(1);
        ind2=[];
    end
    Phase(:,ix,:)=Ph;
    
end
end

function [Phase,Temp] = generate_passive_margin(A,Phase,Temp,Terranes,direction,Gen,limits)

    bla =0; 
    


ind_x = find(A.Xpart(1,:,1)>=limits(1),1);
ind_x2 = find(A.Xpart(1,:,1)<=limits(2));
ind_x2 = ind_x2(end);


    if ~isempty(Terranes.Passive_Margin_segment)

        depo_Y = Terranes.Passive_Margin_depo_center_pos;
        depo_Z = Terranes.Passive_Margin_depo_center;
        Length = Terranes.Passive_Margin_length; 
      
        lim_depo = [Terranes.y_lim(1)-20, Terranes.y_lim(1)+Length]; 
        Depo_y   = lim_depo(2)-depo_Y.*Length; 
        Depo_z   = - depo_Z; 
        y        = [lim_depo(1),Depo_y,lim_depo(2),lim_depo(2)+0.1*Length];
        z        = [1.0, Depo_z,Depo_z,1.0];
        Lithos   = Terranes.Stratigraphy(end);
       % y_l      = [lim_depo(1),lim_depo(2),lim_depo(2)];
        %z_l      = [Lithos,Lithos,Lithos+Terranes.Passive_Margin_d_lithos];
    
    end
    for i = ind_x:ind_x2
        Ay = squeeze(A.Ypart(:,i,:));
        Az = squeeze(A.Zpart(:,i,:));
        Ph = squeeze(Phase(:,i,:));
        T  = squeeze(Temp(:,i,:));
        [in,~] = inpolygon(Ay,Az,y,z);
        Ph(in) = Terranes.Passive_Margin_phase(1); 
%        [in2,~]  = inpolygon(Ay,Az,y_l,z_l);
      %  Ph(in2) = Gen.Ph_UM;
       % T(in2) = Gen.T_P; 
        Phase(:,i,:)=Ph; 
        Temp(:,i,:)=T;
    end
    
    
end



function [string] = save_topography(A,topo_M,Gr)
    
    
    dX = (max(Gr.x_g)-min(Gr.x_g))/1500;
    x_t = min(Gr.x_g):dX:max(Gr.x_g);
    topo_t=interp1(squeeze(A.Xpart(:,1,1)),topo_M,x_t,'nearest','extrap');
    
    Topo = zeros(length(x_t),length(Gr.y_g));
    Topo(:,1)        = topo_t;
    Topo(:,2)        = topo_t;
    Easting     = x_t;
    Northing    = Gr.y_g;
  
    % compute grid spacing
    dx = (max(Easting) - min(Easting))/(length(Easting)-1);
    dy = (max(Northing) - min(Northing))/(length(Northing)-1);
    
    
    % transpose Topo to be read correctly by LaMEM
    %Topo    = Topo';
    
    % write binary to be read by LaMEM
    % (FILENAME,[# of points in x-dir, # of points in y-dir, x-coord of SW corner, y_coord of SW corner,
    % grid spacing in x-dir, grid spacing in y-dir, TopoMatrix in vector form])
    PetscBinaryWrite('topo.dat', [size(Topo,1); size(Topo,2); min(Easting);min(Northing); dx; dy; Topo(:)]);
    string = 'Isostatic balance finished, and the resulted topography has been saved in Topo.dat';
end 


function [A,Gr] = Parse_LaMEM_bin(Parallel_partition,Paraview_output,LaMEM_Parallel_output,npart,Testing_3D)
   %==========================================================================
    % OUTPUT OPTIONS
    %==========================================================================
    % See model setup in Paraview 1-YES; 0-NO

    RandomNoise             =   logical(0);
    Is64BIT                 =   logical(0);
    %==========================================================================
    % LOAD MESH GRID FROM LaMEM PARTITIONING FILE
    %==========================================================================
    npart_x =   npart(1);
    npart_y =   npart(2);
    npart_z =   npart(3);

    % Load grid from parallel partitioning file
    [X,Y,Z,x,y,z, Xpart,Ypart,Zpart] = FDSTAGMeshGeneratorMatlab(npart_x,npart_y,npart_z,Parallel_partition, RandomNoise, Is64BIT);

    Gr.x_g = [min(x),max(x)]; 
    Gr.z_g =[min(z),max(z)];
    Gr.y_g = [min(y),max(y)]; 
    % Update variables (size of grid in [x,y,z] direction
    nump_x  =   size(X,2);
    nump_y  =   size(X,1);
    nump_z  =   size(X,3);  
    W       =   max(x)-min(x);
    mW      =   abs(min(x));
    L       =   max(y)-min(y);
    mL      =   abs(min(y));
    H       =   max(z)-min(z);
    mH      =   abs(min(z)); 
    Xvec    =   squeeze(X(1,:,1));
    Yvec    =   squeeze(Y(:,1,1));
    Zvec    =   squeeze(Z(1,1,:));
    % Temporary save memory
    clear Z Xpart Ypart Zpart
    % Prepare data for visualization/output
    A           =   struct('W',[],'L',[],'H',[],'nump_x',[],'nump_y',[],'nump_z',[],'Phase',[],'Temp',[],'x',[],'y',[],'z',[],'npart_x',[],'npart_y',[],'npart_z',[]);
    % Linear vectors containing coords
    [A.Xpart,A.Ypart,A.Zpart] =meshgrid(single(Xvec),single(Yvec),single(Zvec));
end


function [y]  = linear_margin(xa,theta,y0,x)

% Function to describe the linear margin of a terranes 
theta_r = theta*pi/180;
m       = tan(theta_r);
y       = m.*(x-xa)+y0; 

end


function [theta_c]  = continental_theta_angle(x,xa,xb,theta_c0,theta_c1)

% Function to describe the linear margin of a terranes 
theta_c = theta_c0+(((theta_c1-theta_c0))./(xb-xa)).*(x-xa);

end

function [y]  = parabolic_margin(xa,xv,yv,x)

% Function to describe the linear margin of a terranes 
a = -yv./(xa-xv).^2;
y = a.*(x-xv).^2+yv; 

end

function [y] = circumference_margin(xa,xb,ya,Cx,R,x)

% Check if the user had the brilliant idea to create a circumference whose
% radius is less the the distance between the two points
if R<(xb-xa)/2
    error('Just think to the pytaghora theorem, i.e., you cannot have a triangle rectangle whose hyp is less than one of the other two')
end
% I try to put more general formulation, but the stuff for this particular
% case is easy
c_ = (xa-Cx).^2-R.^2+ya^2;
b_ = 2.*ya; 
delta = sqrt(b_^2-4.*c_);
center_y = [(b_-delta)./2;(b_+delta)./2]; 
y_ind = find(center_y<=ya,1); 
y_c = center_y(y_ind);
% Compute the semi-circumference  (only positive)
y = (R^2-(x-Cx).^2).^(0.5)+y_c; 
end

function [xa,xb,xc,yc] = find_terrane_size(R,arc_length)
theta_ = arc_length./(R); 
xa     = -R.*sin(theta_./2);
xb     = R.*sin(theta_./2);
xc     = 0.0; 
yc     = -R.*cos(theta_./2); 
end

function [s] = arc_length_x(x,y,xa,ya,R)
% compute the arc length of a given point using the Pa - Px points 
d = sqrt((x-xa).^2+(y-ya).^2);
theta = 2*asin(d./(2*R));
s = R.*theta; 
end


function [W] = length_parabolic_margin(xa,xv,yv,x)
a = -yv./(xa-xv).^2;
xb  = -xa; 
t  = (xa-xv);
t2 = (xb-xv); 
W1 = (2.*a.*t.*(4.*a.*a.*t.^2+1)+asinh(2.*a.*t))./(4.*a);
W2 = (2.*a.*t2.*(4.*a.*a.*t2.^2+1)+asinh(2.*a.*t2))./(4.*a);
W = W2-W1;
end

function [yes_no] = check_belonging_layer(var,less_higher,A,axis)
% Function to check if specific ensemble of particles belong to the layer 

if iscell(var)
    if(var{2}==0)
        buf = var{1};
        var =[]; 
        var = buf; 
    else
        fun = var{1};
        if less_higher == 0 
            if axis == 1 
                yes_no = A.Xpart <= fun(A.Ypart);
            else 
                yes_no = A.Ypart <= fun(A.Xpart);
            end
        else
           if axis == 1 
                yes_no = A.Xpart >= fun(A.Ypart);
            else 
                yes_no = A.Ypart >= fun(A.Xpart);
            end
        end
    end
end
if iscell(var) == 0 
   if less_higher == 0 
       if axis == 1 
            yes_no = A.Xpart <= var;
       else 
            yes_no = A.Ypart <= var;
       end

   else
       if axis == 1 
            yes_no = A.Xpart >= var;
       else 
            yes_no = A.Ypart >= var;
       end
   end
end

end

