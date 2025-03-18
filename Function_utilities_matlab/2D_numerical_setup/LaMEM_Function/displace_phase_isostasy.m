
function [A,surf] = displace_phase_isostasy(ph,A,Gr,TI)
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
T2(isnan(T2) & Z >0.0)=TI.TS;
T2(isnan(T2)& Z <0.0)=TI.TP;

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


