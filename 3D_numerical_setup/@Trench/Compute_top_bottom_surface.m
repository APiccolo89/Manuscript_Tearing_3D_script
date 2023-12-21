function [obj] = Compute_top_bottom_surface(obj,theta)
%=========================================================================%
% Function to compute
% Function to compute the top and bottom surface of the slab if the bending
% angle is changing as a function of the lenght.
% Give to the function L0, D0, nseg, theta max, type of law (Ribe,linear)
% L0 is the total length of the slab; D0 is the thickness of the slab; nseg
% is the number of segment to discretize, theta max is the final angle of
% dip of the slab.
% Start at l = 0, where l is the position along the midsurface of the
% slab(l->[0,L0]). Compute ln (l+dl), where dl is L0/nseg. Compute the x,z
% position of the slab at ln, -> compute the top and bottom surface.
% Continue the procedure till l < L0
%------------------------------------------------------------------------%
%=========================================================================%
if nargin == 0
    %=====================================================================%
    % Default Debugging parameters                                        %
    %=====================================================================%
    L0 = 600;

    seg_L = 30; % number of segment

    theta = 60*pi./180; %radians

    D0 = 100;

    ftheta = @(l) ribe_angle(l,L0,theta);

    l_test = 0:L0./(1000):L0;
else

    theta = theta.*pi/180; 

    L0 =obj.L0;

    if isempty(obj.Lb)
        Lb = L0;
    else
        Lb = obj.Lb; 
    end

    D0 = obj.D0; 

    WZ_tk = obj.tk_WZ; 

    seg_L = obj.nseg; % number of segment

    % Choose the function handle
    if strcmp(obj.Type_Angle, 'linear')

        ftheta = @(l) theta_computing_linear(l,Lb,theta);

    elseif strcmp(obj.Type_Angle, 'ribe')

        ftheta = @(l) ribe_angle(l,Lb,theta);

    else 

        error('Type of angle not existent, check the grammar')

    end

end
% Prepare Top,mid and bottom surface.

Top = zeros(seg_L+1,2);
% Initialize the mid surface

MidS = Top;

MidS(1,2) = -D0./2;
% Initialize the bottom surface

Bottom = Top;

Bottom(1,2) = -D0;

% Weak zone top surface
WS_surf = Top;

WS_surf(1,2) = WZ_tk; 

%Initialize iteration
l = 0;

dl = L0./seg_L;

it = 1;
%While Loop routine 

while l<L0

    ln = l+dl;
    % Compute the mean angle within the segment

    theta_mean(it) = (ftheta(l)+ftheta(ln))./2;
    % Compute the mid surface coordinate

    MidS(it+1,1) = MidS(it,1)+dl*cos(theta_mean(it));

    MidS(it+1,2) = MidS(it,2)-dl.*sin(theta_mean(it));
    % Compute the top surface coordinate

    Top(it+1,1) = MidS(it+1,1)+0.5.*D0.*abs(sin(theta_mean(it)));

    Top(it+1,2) = MidS(it+1,2)+0.5.*D0.*abs(cos(theta_mean(it)));
    % Compute the bottom surface coordinate

    Bottom(it+1,1) = MidS(it+1,1)-0.5.*D0.*abs(sin(theta_mean(it)));

    Bottom(it+1,2) = MidS(it+1,2)-0.5.*D0.*abs(cos(theta_mean(it)));
    % Compute the top surface for the weak zone 

    WS_surf(it+1,1) = MidS(it+1,1)+(0.5.*D0+WZ_tk).*abs(sin(theta_mean(it)));

    WS_surf(it+1,2) = MidS(it+1,2)+(0.5.*D0+WZ_tk).*abs(cos(theta_mean(it)));
    % update l and it

    l = ln;

    it = it+1;

end

if nargin == 0

    % Plot the figure
    figure(1)

    clf;

    hold on

    scatter(Top(:,1),Top(:,2),'blue','filled');

    scatter(MidS(:,1),MidS(:,2),'black','filled');

    scatter(Bottom(:,1),Bottom(:,2),'red','filled');
    
    scatter(WS_surf(:,1),WS_surf(:,2),'green')

    hold on

    plot(Top(:,1),ones(length(Top(:,1)),1).*(0-100),'LineWidth',2.0)

    axis equal

end

% Import the important information about the slab surface. 

Slab_surface.Top = Top;

Slab_surface.MidS = MidS;

Slab_surface.Bottom = Bottom;

Slab_surface.WZ_surf = WS_surf; 
%Update object

obj.Slab_surface = Slab_surface; 

end

% functions for computing the bending angle of the slab

% Linear function angle

function [theta_l] = theta_computing_linear(l,L0,theta)

% slope linear function

s = (theta-0)./(L0);

theta_l = l.*s;
if l>L0
    theta_l=theta;
end

end

%Ribe angle function

function [theta_l] = ribe_angle(l,L0,theta)

theta_l = theta.*l.^2.*((3.*L0-2.*l))./(L0^3);
if l>L0
    theta_l=theta;
end

end