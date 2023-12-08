
function [A] = displace_phase_isostasy(ph,A,Gr,TI)
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
[rho_ij] = Compute_density(A,ph);
[topo,~] = compute_topography_(A,rho_ij);
%==========================================================================
% Compute the length scale of the moving average using a slab of 90 km as
% reference.
%==========================================================================
x  = squeeze(A.Xpart(:,:,1));
y  = squeeze(A.Ypart(:,:,1));
% My bad and ignorance: I apply moving average along x, and y, with a
% number of nodes that grossly is equivalent to 200 km (two times more or
% less a typical slab).

dx = diff(x(:,1));
dy = diff(y(1,:));
dx_MV = floor(100./mean(dx));
dy_MV = floor(100./mean(dy));
topo_M = movmean(topo,dx_MV,2);
topo_M = movmean(topo_M,dy_MV,1);
% Substract the median of the topography:
% Here I would like to open a philosophical discussion on the significance
% of the base level and how to conceive a continental lithosphere. Atm, the
% average continental topography is above the sea level. One cool thing is
% adding the water to LaMEM solving the main issue with this matter.
topo_M = topo_M-median(topo_M,"all");
Ph = squeeze(A.Phase(:,:,:));
T  = squeeze(A.Temp(:,:,:));
Z  = squeeze(A.Zpart(:,:,:));
%x  = squeeze(A.Xpart(:,1,:));
%y  = squeeze(A.Ypart(:,1,:));
ilx = length(x(:,1));
ily = length(y(1,:));
Ph2 = 0.*Ph;
T2  = 0.*Ph;
% loop over the 2D nodes, interpolate along z
% This is where most of the time is lost: i cannot do better, you need to
% port this stuff into Julia or C for having a decent performance: Matlab
% is shit.
% is shit.
for i = 1:ilx
    for j =1:ily
        z = squeeze(Z(i,j,:))+topo_M(i,j);
        Ph2(i,j,:)=interp1(z,squeeze(Ph(i,j,:)),squeeze(Z(i,j,:)),'nearest');
        T2(i,j,:)=interp1(z,squeeze(T(i,j,:)),squeeze(Z(i,j,:)),'nearest');
        z = [];
    end
end
Ph2(isnan(Ph2(:))& Z(:) >0.0)=ph.Ph_Ar(1);
Ph2(isnan(Ph2(:))& Z(:) <0.0)=ph.Ph_UM(1);
T2(isnan(T2(:)) & Z(:) >0.0)=TI.TS;
T2(isnan(T2(:))& Z(:) <0.0)=TI.TP;
A.Phase = reshape(Ph2,size(A.Xpart));
A.Temp = reshape(T2,size(A.Xpart));
plot_section(A,topo_M)


[s]=save_topography(A,topo_M,Gr);
disp(s);

end

function [topo,rho] = compute_topography_(A,rho_ij)

x = (squeeze(A.Xpart(:,1,1)));
y = (squeeze(A.Ypart(1,:,1)));
topo  = zeros(length(x),length(y));
rho= nanmean(rho_ij,[3]);
rho_0 = rho(1,1);
CPU_A = cputime;
topo = 1000.*(rho_0-rho)./rho;
CPU_B = cputime;
disp(['Raw Topography took, ',num2str(CPU_B-CPU_A,3),' sec']);
end



function [string] = save_topography(A,topo_M,Gr)
y=squeeze(A.Ypart(1,:,1));
x=squeeze(A.Xpart(:,1,1));
% %%g
dX = (max(Gr.x_g+0.1)-min(Gr.x_g-0.1))/(Gr.lx-1);
dY = (max(Gr.y_g+0.1)-min(Gr.y_g-0.1))/(Gr.ly-1);
x_t = min(Gr.x_g-0.1):dX:max(Gr.x_g+0.1);
y_t = min(Gr.y_g-0.1):dY:max(Gr.y_g+0.1);
[Y_t,X_t] = meshgrid(y_t,x_t);

[Y,X] = meshgrid(y,x);
% How to make fast this piece of art?
[topo_t] = griddata(double(Y(:)),double(X(:)),double(topo_M(:)),Y_t,X_t);
topo_t(1,:) = topo_t(2,:);
topo_t(end,:)=topo_t(end-1,:);
topo_t(:,1) = topo_t(:,2);
topo_t(:,end)=topo_t(:,end-1);
% It appears that Matlab is really lobotomized, and does not extrapolate in
% the corner, implying that I need to do yet an other loop, but hey, we
% have a lot of time to invest, waiting ages for doing simple task, because
% a multimilionare company is too lazy to implement this feature.
% So the best solution is to linear interpolation of the marker to a
% regular grid:
Topo = topo_t;
Easting     = x_t;
Northing    = y_t;

% compute grid spacing
dx = (max(Easting) - min(Easting))/(length(Easting)-1);
dy = (max(Northing) - min(Northing))/(length(Northing)-1);


% transpose Topo to be read correctly by LaMEM
Topo    = Topo;

% write binary to be read by LaMEM
% (FILENAME,[# of points in x-dir, # of points in y-dir, x-coord of SW corner, y_coord of SW corner,
% grid spacing in x-dir, grid spacing in y-dir, TopoMatrix in vector form])
PetscBinaryWrite('Topo.dat', [size(Topo,1); size(Topo,2); min(Easting);min(Northing); dx; dy; Topo(:)]);
string = 'Isostatic balance finished, and the resulted topography has been saved in Topo.dat';
end


function [rho_ij] = Compute_density(A,ph)
rho_ij = A.Xpart.*0.0;

% loop over the phases contained in ph
field_names = fields(ph);
ip = numel(field_names);
CPU_PRIME = cputime;
for i = 1:ip
    CPU_A = cputime;
    % Select Phase
    current_phase = ph.(field_names{i});
    % Extract unrelevant information such as the reference density
    Density_ref = current_phase(2);
    % Compute the density as a function of the temperature
    if current_phase(1) ~=0
        rho_ij(A.Phase==current_phase(1)) = density_node(Density_ref,A.Temp(A.Phase==current_phase(1)));
    else
        rho_ij(A.Phase==current_phase(1))=Density_ref;
        disp('Air does not deserve to be computed, but, as a quick reminder:'); disp('do not set 0 friction angle or cohesion to air, you will not notice my help.');disp(' But believe me KSP solver will be grateful and gives candy')
    end
    % Just random information for losing time on Monday, when I am
    % supposed to take care of myself.
    CPU_B = cputime;

    disp(['Phase nr. ', num2str(current_phase(1)),' has been processed in ',num2str(CPU_B-CPU_A,2), 'seconds']);
end
disp('Additional note: density is computed assuming a thermal expansion of:');
disp('3e-5 [1/K], if you want to change, just go to function density node, in displace_isostasy.m');
disp(['Average density  has been processed in ',num2str(CPU_B-CPU_PRIME,2), 'seconds']);


end

function [rho_n] = density_node(Density_ref,T)
%=================================================================%
% Hard coded stuff, but I doubt that anyone would really enjoy
% himself on the daunting task to create ad hoc alpha for each
% phase. There are worst crime in our community. For example melt
% extraction or other funny approximation (like Venus with a
% surface temperature of 0 degree celsius)
%=================================================================%
% Prior: selected Phases -> select Markers -> send T
% informationhere and phase information here: compute the density
%=================================================================%
rho_n = Density_ref.*(1-3e-5.*(T-25.0));
end

function plot_section(A,topo_M)
% This function plot 10 section along x,y it plot every 100 nodes,
% temperature
X = squeeze(A.Xpart(:,1,1));
Y = squeeze(A.Ypart(1,:,1));
Z = squeeze(A.Zpart(1,1,:));
lx = length(X);
ly = length(Y);
it = 0;
Ph = A.Phase;
Ph(Ph==0)=nan;
figure(10)
ax = gca;
s = pcolor(X,Y,topo_M');
shading interp;
ax.XColor = [0,0,0];
ax.YColor = [0,0,0];
ax.Box= 'on';
ax.TickLabelInterpreter = 'latex';
string_title = ('Topography, [km]');
ax.Title.String = string_title;
ax.Title.Interpreter = 'latex';
ax.XLabel.String = 'X, [km]';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.String = 'Y, [km]';
ax.YLabel.Interpreter = 'latex';
ax.ZLabel.String = 'H, [km]';
ax.ZLabel.Interpreter = 'latex';
ax.YLim = [-600, 600];
ax.XLim = [-600,600];
colorbar;

path_folder = 'Initial_Setup';
colorbar;
if ~isfolder(path_folder)
    mkdir(path_folder)
end

name_pic = (['Topography',num2str(it),'.png']);
filename = fullfile(path_folder,name_pic);
print(filename,'-dpng')






for i = 1:20:lx
    it = it+1;
    figure(it)
    clf;
    %======================================%
    % Section along x                      %
    %======================================%
    ax = gca;
    p1=pcolor(Y,Z,squeeze(Ph(i,:,:))');
    shading interp;
    ax.XColor = [0,0,0];
    ax.YColor = [0,0,0];
    ax.Box= 'on';
    ax.TickLabelInterpreter = 'latex';
    string_title = (['Phase Field, section x = ',num2str(X(i),2), ' km']);
    ax.Title.String = string_title;
    ax.Title.Interpreter = 'latex';
    ax.XLabel.String = 'Y, [km]';
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.String = 'Z, [km]';
    ax.XLabel.Interpreter = 'latex';
    path_folder = 'Initial_Setup';
    colorbar;
    if ~isfolder(path_folder)
        mkdir(path_folder)
    end
    path_save = 'Initial_Setup/Phase';
    if ~isfolder(path_save)
        mkdir(path_save)
    end
    name_pic = (['Phase',num2str(it),'.png']);
    filename = fullfile(path_save,name_pic);
    print(filename,'-dpng')



    figure(it)
    clf;
    %======================================%
    % Section along x                      %
    %======================================%
    ax = gca;
    p1=pcolor(Y,Z,squeeze(A.Temp(i,:,:))');
    shading interp;
    ax.XColor = [0,0,0];
    ax.YColor = [0,0,0];
    ax.Box= 'on';
    ax.TickLabelInterpreter = 'latex';
    string_title = (['Temp Field, section x = ',num2str(X(i),2), ' km']);
    ax.Title.String = string_title;
    ax.Title.Interpreter = 'latex';
    ax.XLabel.String = 'Y, [km]';
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.String = 'Z, [km]';
    ax.XLabel.Interpreter = 'latex';
    path_folder = 'Initial_Setup';
    caxis([20,1300]);
    colorbar;
    if ~isfolder(path_folder)
        mkdir(path_folder)
    end
    path_save = 'Initial_Setup/Temp';
    if ~isfolder(path_save)
        mkdir(path_save)
    end
    name_pic = (['Temp',num2str(it),'.png']);
    filename = fullfile(path_save,name_pic);
    print(filename,'-dpng')






end
close all;

end

function [int_val] = linear_interpolation_marker_grid(xp,yp,x,y,I)
%====================================================================%
% Input:
% I = Quantity to interpolate
%x,y regular grid
%xp,yp marker position
% Output:
%int_val => value at the nodal point
%===================================================================%
lx = length(x);
ly = length(y);
int_val = zeros(ly,ly);
dx = diff(x);
dy = diff(y);
dx = dx(1);
dy = dy(1);
% Transform coordinate:
xp = xp-min(x);
yp = yp-min(y);
x = x-min(x);
y = y-min(y);
%======================================
[Y,X]=meshgrid(y,x);
X = X(:);
Y = Y(:);
CPU_timeA = cputime;
for i = 1:length(X)
    mx = abs(xp-X(i));
    my = abs(yp-Y(i));
    m_ind = mx<=dx & my <=dy;
    weight  = (1-mx(m_ind==1)./(dx)).*(1-my(m_ind==1)./(dy));
    int_val(i) = sum(I(m_ind==1).*weight)./sum(weight);
    weight = [];
end
CPU_timeB = cputime;
disp(['Interpolation topography took', num2str(CPU_timeB-CPU_timeA,3), ' seconds: could be worse.'])
end
