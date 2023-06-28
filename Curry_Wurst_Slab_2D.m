
clear all;
close all;
addpath matlab   %Here -> Folder2LaMEM/matlab => all the function that handle matlab files are there


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
T.theta    = 30;   % curvature slab
T.tk_WZ    = 0;   % thickness of the weak zone
T.L0       = 500;  % length of the slab from the bottom of the lithosphere
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
ind = Phase == 0 & A.Zpart<0.0;
Phase(ind)  = Gen.Ph_UM;
Temp(ind)   = Gen.T_P;
ind = Phase == 0 & A.Zpart>0.0;
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

[A,surf] = displace_phase_isostasy(ph,A,Gr,Gen);

plot_initial_setup2D(A,surf);

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
cont = [];
A_fill_layer = cputime;
[Phase,Temp] = fill_layer(A,Terranes,Phase,Temp,Gen,cont);
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
        [Phase,Temp] = generate_passive_margin(A,Phase,Temp,Terranes,direction,Gen);
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
Tsl   = Gen.AvTS;


if Type == 1
    Zpart = A.Zpart; 
    ind_z = find(A.Zpart<T_tk(1) & A.Zpart>=T_tk(2) & indx == 1);
    ind = A.Zpart< T_tk(2);

else
    ind_z = A < T_tk(1) & A>=T_tk(2);
    ind = A < T_tk(2);
    Zpart = A; 
    

end
if isnan(T_Age)
    Temp(ind_z)=Tsl;
else

    erf_function = erf(Zpart(ind_z).*1000/2/(k*T_Age)^0.5);
    Temp(ind_z) = T_S - (T_P-T_S).*erf_function;
    Temp(ind) = T_P;

end




end

function [Temp]= compute_temperature_profile_McKenzie(Layout,Temp,Type,Gen,Terranes,indx,A)

T_tk = [0.0, Terranes.Stratigraphy(end)];
k = Terranes.K;
T_Age = Terranes.Age.*Terranes.secMyrsyear;
T_P   = Gen.T_P;
T_S   = Gen.T_S;
Tsl   = Gen.AvTS;
vl    = 10; 
vl    = vl./1e2./365.25./24./60./60; 
Re    = (3300*Terranes.Cp.*vl*Terranes.Trench_properties.D0.*1000)./2./k;
% Compute the lenght of the slab starting from the trench
%1) Prepare the vector
z     = A.Zpart(:);
x     = A.Xpart(:);
%2) Select index that belongs to slab 
d     = Layout(:);                    %layout filled with d
slab  = abs(d)>=0 & abs(d)<=Terranes.Trench_properties.D0; % index layout
l     = d.*0.0;  %arc_lenght of the slab 
%3) Extract the variable concerning the center of circumference 
C     = Terranes.Trench_properties.C;   %center of circumference
R     = Terranes.Trench_properties.r;%Radius
%4) Extract the properties of the bending angle 
theta = Terranes.Trench_properties.theta*pi/180;% theta bending
%5) Find the midplane radius 
r    = sum(R)./length(R);% radius of mid plane 
%6) Compute the angle 
angle = 0.0*d; % angle layout
D     = angle; % buffer vector containing the distance w.r.t. to the radius of curvature
D(slab==1) = sqrt((x(slab==1)-C(1)).^2+(z(slab==1)- C(2)).^2);
angle(slab==1)=acos((z(slab==1)-C(2))./D(slab==1)); % angle 
% Empty the garbage 
D      = []; 
%7) Compute the length along the curved area 
l(slab==1 & angle<=theta)=r.*angle(slab==1 & angle<=theta); 
%8) Compute the length along the linear part. 
% Coordinate of the point of midsurface after curvature
PA(1)=C(1)+r*sin(theta);
PA(2)=C(2)+r*cos(theta); 
%Coordinate of the midsurface at the depth of slab 
PB(1)=PA(1)+Terranes.Trench_properties.L0*cos(theta);
PB(2)=PA(2)-Terranes.Trench_properties.L0*sin(theta); 
% Compute the slope of the line: 
m=(PB(2)-PA(2))./(PB(1)-PA(1)); 
%m = 1./m;
% project all the x points into the line (from:
% https://math.stackexchange.com/questions/62633/orthogonal-projection-of-a-point-onto-a-line)
% 
xprojection = angle.*0.0; 
zprojection = angle.*0.0;
b = PA(2)-PA(1).*m;
xprojection(slab==1 & l==0.0) = (x(slab==1 & l == 0.0) + m.*z(slab==1 & l==0.0)-m.*b)./(1+m.^2); 
zprojection(slab==1 & l==0.0) = (m.*x(slab==1 & l == 0.0) + m^2.*z(slab==1 & l==0.0)+b)./(1+m.^2);
l(slab==1 & l==0.0) = r*theta+sqrt((xprojection(slab==1 & l==0.0)-PA(1)).^2+(zprojection(slab==1 & l==0.0)-PA(2)).^2);

LL = max(l(slab==1)); 
% make dimensionless d and l; 
l = l./Terranes.Trench_properties.D0; 
d = d./Terranes.Trench_properties.D0;
%======== Compute the temperature profile
n = 50; 
t_prov = Temp(:).*0.0; 
Sigma    = t_prov; 
Sigma(slab==1) = 0.0; 
% Dividi et impera
dn = 1-abs(d(slab==1));
for i=1:n
    a = (-1).^(i)./(i.*pi);
    b = (Re-(Re^2+i^2*pi^2).^(0.5)).*l(slab==1);
    c = sin(i.*pi.*(1-abs(d(slab==1))));
    e = exp(b);
    Sigma(slab==1) = Sigma(slab==1)+ a.*e.*c;  
end
    t_prov(slab==1) = (T_P)+2.*(T_P-T_S).*Sigma(slab==1); 
    t_prov(slab==1 & t_prov<0) = T_S; 
    weight = abs(z./Terranes.Trench_properties.D0);
    weight(weight>1) = 0.9; 
    Temp(slab==1)=(Temp(slab==1).*(1-weight(slab==1))+t_prov(slab==1)).*weight(slab==1); 
    Temp = reshape(Temp,size(A.Xpart)); 
end




function [Layout,arc_angleS,continent] = find_slab_(A,Slab,type)
%
%
% Really convoluted, but it works:
% 1) create a polygon for the slab
% 2) find the curvature
% 3) save the relative distannce w.t.r the top
%
%
% find the area belonging to the curvature slab:
Xvec = A.Xpart(1,:,1);
Yvec = A.Ypart(:,1,1);
Zvec = A.Zpart(1,1,:);
C     = Slab.C;
D0    = abs(Slab.D0);




theta = Slab.theta;
theta_c = Slab.theta_c;
ang_  = (theta)*pi/180;
ang_2 = (theta+90)*pi/180;
L0    = Slab.L0 - Slab.R;


r     = Slab.r;
L0    = Slab.L0 - Slab.R;
sl    = 1               ;
arc_angleS = [];
if ~strcmp(type,'Mode_1')
    ind = (A.Xpart(:)>=Slab.C(1) & A.Xpart(:)<=Slab.C(1)+2.*L0) & A.Zpart(:)>=-1.5*L0 & A.Zpart(:)<=1.0;
    A_time = cputime;
    x     = A.Xpart(:);
    z     = A.Zpart(:);
    y     = A.Ypart(:);
    %x     = x(ind==1);
    %y     = y(ind==1);
    %z     = z(ind==1);
    % Vector version
    Layout = x*nan;
    d      = Layout;

    arc_angleS = x.*NaN;
    continent = arc_angleS;
    angle_c = continent;
    %s = fhandle_arc_length(x,fhandle(x));

    u = zeros(2,length(x));
    v = u;
    u(1,ind==1) = x(ind==1)';
    u(2,ind==1) = z(ind==1)';
    d(ind==1) = sqrt((u(1,ind==1)-Slab.C(1)).^2+(u(2,ind==1)- Slab.C(2)).^2);
    arc_angleS(ind==1) = acos((z(ind==1)-Slab.C(2))./d(ind==1)).*180/pi;
    c_slab = d(ind==1)>=r(1) & d(ind==1)<=(r(2)) & arc_angleS(ind==1)<=theta & arc_angleS(ind==1)>=0;
    angle_c(ind==1) = d(ind==1)>=r(2)-Slab.Cont+((Slab.Cont-0.0)./(theta_c)).*arc_angleS(ind==1) & d(ind==1)<=(r(2)); %& arc_angleS(ind==1)<=theta_c ;
    continent(angle_c==1)=1.0;
    Layout= d;
    Layout(d<r(1) | d>r(2)) = NaN;
    Layout(arc_angleS>theta | arc_angleS<0)=NaN;
    Layout = Layout-r(2);
    if strcmp(type,'Weak')
        r(1) = r(2);
        r(2) = r(2)+Slab.tk_WZ;
        if Slab.theta == 90
            sl = 0;
        end
        %=============

    end

    C(1)     =Slab.C(1);
    C(2)     = Slab.C(2);
    p1(1)    = C(1)+r(1)*sin(ang_);
    p1(2)    = C(2)+r(1)*cos(ang_);
    p2(1) = C(1)+r(2)*sin(ang_);
    p2(2) = C(2)+r(2)*cos(ang_);
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
    [in,~] = inpolygon(x,z,Py,Pz);
    P = [x(in), z(in)];
    Layout(in) = -(find_distance_linear(P,p2,p4))';
    B_time = cputime;
    disp(['Time Loop = ', num2str(B_time-A_time)]);
    Layout = reshape(Layout,size(A.Xpart));
    continent = reshape(continent,size(A.Xpart));

end

if strcmp(Slab.Type,'Ribe_Mode')
    x_p    = A.Xpart(:); % Vector of X
    z_p    = A.Zpart(:); % Vector of Y
    x_t = x_p.*nan;      % Vector X position along ribe curve.
    Layout = x_p*nan;
    arc_angleS =[];
    continent = Layout*nan;
    D_c  = Slab.Depth_C;
    weak_flag = 0;

    if strcmp(type,'Weak')
        weak_flag = 1;
    end
    % New method: 
    % I was not the brightest candle in the chandelor, and I believed to
    % getting rid of rotation matrix by using a convoluted method. However
    % this methods miserably fail in certain condition

    % [A] : coordinate slab tip of the slab 
    SA = [C(1)+L0;0];
    SB = [C(1)+L0;-D0]; 
    % [B]: Find the coordinate of the bended slab 
    theta_r = theta*pi/180; 
    rotation_matrix_0 = [cos(theta_r) sin(theta_r); -sin(theta_r) cos(theta_r)];
    SAB = rotation_matrix_0*SA;
    SBB = rotation_matrix_0*SB; 
    %[C]: Find the limit of the points that are worth to check
    z_lim = min(SBB(2),SAB(2));
    x_lim = max(SAB(1),SBB(1));
    % [D]: create index array of the choosen point
    ip = z_p >= -L0 & z_p<=0.0 & x_p>=C(1) & x_p < C(1)+L0;
    % [E] Find the angle of these points 
    x_pb = x_p(ip==1);
    z_pb = z_p(ip==1);
    % [F] Compute the bending angle. 
    [~,~,~,~,teta_s] = compute_ribe_angles(x_pb,C,theta,L0,D0,weak_flag,Slab.tk_WZ,C(1)+L0);
 %   P                = [x_pb;z_pb]; 
    x_rot = x_pb.*cos(teta_s)-z_pb.*sin(teta_s);
    z_rot = x_pb.*sin(teta_s)+z_pb.*cos(teta_s); 
    % Select the point that are worth of the god love
    ichosen =  x_rot < C(1)+L0 & z_rot>=-80 & z_rot<=0;
    % Since i rotate the plate accordingly i save the depth of the plate 
    z_pb(ichosen)=z_rot(ichosen);
    z_pb(ichosen==0)=nan; 
    x_pb(ichosen)=x_rot(ichosen);
    x_pb(ichosen==0)=nan;
    LLL = Layout; 
    LLL(ip==1)=x_pb;
    Layout(ip==1) = z_pb; 
    continent(z_p>=D_c)  = 1.0;


    Layout = reshape(Layout,size(A.Xpart));
    continent = reshape(continent,size(A.Xpart));

end

end
function [z,dz,ztop,zbot,teta_s] = compute_ribe_angles(xp,C,teta_0,L,D0,weak_flag,Tk,Lim)
%=========================================================================%
% TO DO A complete explanation
%=========================================================================%
teta_0 = teta_0*pi/180;

x = xp-C(1);

teta_s               = teta_0.*x.^2.*((3.*L-2.*x))./(L^3);
z                    = -tan(teta_s).*x;
dz                   = (6.*teta_0.*x.^2.*(x-L).*sec(teta_s).^2)./(L^3)-tan(teta_s);
zbot                 = z-D0/2.*cos(teta_s);
ztop                 = z+D0/2.*cos(teta_s);
if weak_flag == 1
    zbot                 = z+D0/2.*cos(teta_s);
    ztop                 = z+(D0/2+Tk).*cos(teta_s);
end


end

function [d] = find_distance_linear(P,p1,p2)

% Ok, wiki saves my day because I was too lazy:
A  = abs((p2(1)-p1(1)).*(p1(2)-P(:,2))-(p1(1)-P(:,1)).*(p2(2)-p1(2)));
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

function [Phase,Temp] = fill_layer(A,Terranes,Phase,Temp,Gen,cont)
% Select the layer and portion of the model:
indx = A.Xpart >= Terranes.x_lim(1) & A.Xpart < Terranes.x_lim(2);
% Compute the thermal profile
[Temp] = compute_temperature_profile(A,Temp,1,Gen,Terranes,indx);
% Fill the phase stratigraphy:
[Phase] = fill_stratigraphy(A.Zpart,Phase,Terranes,indx,cont);
end

function [Phase,Temp] = fill_subduction(A,Terranes,Phase,Temp,Gen)
ind_x=[];
cont = [];

[Layout,~,cont] = find_slab_(A,Terranes.Trench_properties,'Slab');
% Set the temperature
if strcmp(Terranes.Trench_properties.Temperature,'McKenzie')
    [Temp] = compute_temperature_profile(Layout,Temp,-1,Gen,Terranes,ind_x);%(A,Temp,Type,T_tk,Gen,Terranes,indx
    [Temp] = compute_temperature_profile_McKenzie(Layout,Temp,-1,Gen,Terranes,ind_x,A);%(A,Temp,Type,T_tk,Gen,Terranes,indx)
else
    [Temp] = compute_temperature_profile(Layout,Temp,-1,Gen,Terranes,ind_x);%(A,Temp,Type,T_tk,Gen,Terranes,indx)
end
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
[Phase] = fill_stratigraphy(Layout,Phase,Terranes,ind_x,cont);
[Layout] = find_slab_(A,Terranes.Trench_properties,'Weak');
if strcmp( Terranes.Trench_properties,'Mode_1')
    ind =(Layout<=0.0 & A.Zpart >= Terranes.Trench_properties.D_WZ & Phase ~=Gen.PrismPh);
    Phase(ind) = Gen.WZ;
else
    ind =(Layout>0.0 & A.Zpart >= Terranes.Trench_properties.D_WZ & Phase ~=Gen.PrismPh & A.Zpart<0.0);
    Phase(ind) = Gen.WZ;
end

end


function [Phase] = fill_stratigraphy(Z,Phase,Terranes,indx,cont)
T_Tk = [0.0, Terranes.Stratigraphy(end)];
T = T_Tk(1);


for i=1:length(Terranes.Phases)
    if (i == length(Terranes.Phases))
        B=T_Tk(2);
    else
        B = Terranes.Stratigraphy(i+1);
    end
    if ~isempty(indx)
        ind = Z < T & Z>=B & indx>0 & Phase == 0;
    else
        ind = Z < T & Z>=B;
    end
    Phase(ind) = Terranes.Phases(i);
    T=B;
    ind = [];
end
if ~isempty(cont)
    TT = Terranes.Trench_properties.CCS;
    T = TT.Stratigraphy(1);
    for i=1:length(TT.phases)
        if (i == length(TT.phases))
            B=TT.Stratigraphy(end);
        else
            B = TT.Stratigraphy(i+1);
        end

        ind = Z < T & Z>=B & cont>0;

        Phase(ind) = TT.phases(i);
        T=B;
        ind = [];
    end
end

end


function  [Phase] = generate_accretion_prism(A,Terranes,Phase)

C = Terranes.Trench_properties.C;
d_p = [C(1)+Terranes.position_end_prism 0.0];
s   = (d_p(2)-C(2))./(d_p(1)-C(1));
Phases = Terranes.Phases;
for i = 1:length(Phases)
    ind2 = A.Zpart>s.*(A.Xpart-C(1))+C(2)  & Phase == Phases(i);
    Phase(ind2) = Terranes.prism_phase(1);
end

end

function [Phase,Temp] = generate_passive_margin(A,Phase,Temp,Terranes,direction,Gen)

bla =0;

if ~isempty(direction)
    depo_X = Terranes.Passive_Margin_depo_center_pos;
    depo_Z = Terranes.Passive_Margin_depo_center;
    Length = Terranes.Passive_Margin_length;
    if strcmp(direction,'left')

    else

        lim_depo = [Terranes.x_lim(2)-Length, Terranes.x_lim(2)];
        Depo_x   = lim_depo(1)+depo_X.*Length;
        Depo_z   = - depo_Z;
        x        = [lim_depo(1),Depo_x,lim_depo(2),lim_depo(2)+0.1*Length];
        z        = [0.0, Depo_z,Depo_z,0.0];
        Lithos   = Terranes.Stratigraphy(end);
        x_l      = [lim_depo(1),lim_depo(2),lim_depo(2)];
        z_l      = [Lithos,Lithos,Lithos+Terranes.Passive_Margin_d_lithos];

    end
    [in,~] = inpolygon(A.Xpart,A.Zpart,x,z);
    Phase(in) = Terranes.Passive_Margin_phase(1);
    [in2,~]  = inpolygon(A.Xpart,A.Zpart,x_l,z_l);
    Phase(in2) = Gen.Ph_UM;
    Temp(in2) = Gen.T_P;


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



function plot_initial_setup2D(A,surf)
figure(1)
ratioX = max(A.Xpart(:,1,1))-min(A.Xpart(:,1,1));
ratioZ = max(A.Zpart(1,1,:))- min(A.Zpart(1,1,:));
r = [1,ratioZ./ratioX,1];
x = squeeze(A.Xpart(:,1,:));
z = squeeze(A.Zpart(:,1,:));
T = squeeze(A.Temp(:,1,:));
ph = squeeze(A.Phase(:,1,:));
T(ph==0)=nan;
ph(ph==0)=nan;
hold on
plot(surf(1,:),surf(2,:),"Color",'r','LineStyle','--','LineWidth',1.2)
pcolor(x,z,T);
hold off
shading interp;
pbaspect(r);
xlabel('x, [km]', Interpreter='latex');
ylabel('z, [km]',Interpreter='latex');
title('Temperature Field [$^\circ$ C]',Interpreter='latex');
grid on;
colormap(crameri('bilbao',13));
box on
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex';
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';
c = colorbar;
caxis([0 1350])
c.Label.String = 'Temperature $[^\circ C]$';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter='latex';
print('Temperature','-dpng','-r0')


figure(2)
lv = 1:12;
hold on
plot(surf(1,:),surf(2,:),"Color",'r','LineStyle','--','LineWidth',1.2)
contourf(x,z,ph,13);

hold off
shading interp;
pbaspect(r);
xlabel('x, [km]', Interpreter='latex');
ylabel('z, [km]',Interpreter='latex');
title('Phase Field [$^\circ$ C]',Interpreter='latex');
grid on
caxis([1,13]);
box on
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex';
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';
xlim([-300 300]);
ylim([-100 10])
colormap(crameri('nuuk',13))
c = colorbar;%(['Upper Crust 1','Lower Crust 1', 'Lithospheric Mantle 1', ' Lithospheric Mantle 2', 'Upper Mantle','Oceanic Slab', 'Oceanic Crust', 'Weak Zone', 'Upper Crust 2', 'Lower Crust 2', 'Oceanic Sed.', 'Passive/Prism']);
c.Label.String = 'Phase ';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter='latex';
print('Phase','-dpng','-r0')

figure(3)
plot(surf(1,:),surf(2,:),"Color",'r','LineStyle','--','LineWidth',1.2)
xlim([-600 600])
bla = 0;
end
function [A,Gr] = Parse_LaMEM_bin(Parallel_partition,Paraview_output,LaMEM_Parallel_output,npart)
%==========================================================================
% OUTPUT OPTIONS
%==========================================================================
% See model setup in Paraview 1-YES; 0-NO

Parallel_partition     = 'ProcessorPartitioning_8cpu_4.1.2.bin';
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
function [D] = arc_length_ribe(x,teta_0,L)
teta_s  = teta_0.*x.^2.*((3.*L-2.*x))./(L^3);
dz  = (6.*teta_0.*x.^2.*(x-L).*sec(teta_s).^2)./(L^3)-tan(teta_s);
D    = (1+dz.^2).^0.5;
end
% min zero -> find L0, arc_length
function [res] = find_the_length_slab(theta,L,arc_length,C)
f_handle = @(x) arc_length_ribe(x,theta*pi/180,L);
res = integral(f_handle,C,L)-arc_length;
end
function [Length] = compute_the_arc_integral(theta,x_p,ip,L,C)
    f_handle = @(x) arc_length_ribe(x,theta*pi/180,L);
    x_pb = x_p(ip==1);
    X = unique(x_pb);
    for i = 1:length(X)
       
        Buf(i) = integral(f_handle,C,X(i));
        
    end


    bla = 0;
end