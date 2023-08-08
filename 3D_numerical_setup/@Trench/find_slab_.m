function [obj,Phase,Temp] = find_slab_(obj,A,Weak_Slab,Phase,Temp,Boundary)

if strcmp(obj.Type_Subduction,'Mode_1')
    [obj,Phase,Temp] = obj.find_slab_mode_1(A,Weak_Slab,Phase,Temp,Boundary,theta);
elseif strcmp(obj.Type_Subduction,'Ribe')
    [obj] = obj.find_slab_ribe(A,Weak_Slab);
else
    error('Oh pityful soul, why did you not read the not existant user guide. The two slab types are {Mode_1}, and {Ribe}')
end
end




% Ribe mode must be rethinked and properly done. This is not yet good
function [obj] = find_slab_ribe(obj,A)
x_p    = A.Xpart(:); % Vector of X
z_p    = A.Zpart(:); % Vector of Y
x_t = x_p.*nan;      % Vector X position along ribe curve.
Layout = x_p*nan;
arc_angleS =[];
continent = Layout*nan;
D_c  = obj.Depth_C;
weak_flag = 0;

if strcmp(type,'Weak')
    weak_flag = 1;
end
% New method:
% I was not the brightest candle in the chandelor, and I believed to
% getting rid of rotation matrix by using a convoluted method. However
% this methods miserably fail in certain condition

% [A] : coordinate slab tip of the slab
SA = [obj.C(1)+obj.L0;0];
SB = [obj.C(1)+obj.L0;-obj.D0];
% [B]: Find the coordinate of the bended slab
theta_r = theta*pi/180;
rotation_matrix_0 = [cos(theta_r) sin(theta_r); -sin(theta_r) cos(theta_r)];
SAB = rotation_matrix_0*SA;
SBB = rotation_matrix_0*SB;
%[obj.C]: Find the limit of the points that are worth to check
z_lim = min(SBB(2),SAB(2));
x_lim = max(SAB(1),SBB(1));
% [D]: create index array of the choosen point
ip = z_p >= -obj.L0 & z_p<=0.0 & x_p>=obj.C(1) & x_p < obj.C(1)+obj.L0;
% [E] Find the angle of these points
x_pb = x_p(ip==1);
z_pb = z_p(ip==1);
% [F] Compute the bending angle.
[~,~,~,~,teta_s] = compute_ribe_angles(x_pb,obj.C,theta,obj.L0,obj.D0,weak_flag,obj.tk_WZ,obj.C(1)+obj.L0);
%   P                = [x_pb;z_pb];
x_rot = x_pb.*cos(teta_s)-z_pb.*sin(teta_s);
z_rot = x_pb.*sin(teta_s)+z_pb.*cos(teta_s);
% Select the point that are worth of the god love
ichosen =  x_rot < obj.C(1)+obj.L0 & z_rot>=-80 & z_rot<=0;
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



function [D] = arc_length_ribe(x,teta_0,L)
teta_s  = teta_0.*x.^2.*((3.*L-2.*x))./(L^3);
dz  = (6.*teta_0.*x.^2.*(x-L).*sec(teta_s).^2)./(L^3)-tan(teta_s);
D    = (1+dz.^2).^0.5;
end
% % min zero -> find L0, arc_length
% function [res] = find_the_length_slab(theta,L,arc_length,C)
% f_handle = @(x) arc_length_ribe(x,theta*pi/180,L);
% res = integral(f_handle,C,L)-arc_length;
% end
% function [Length] = compute_the_arc_integral(theta,x_p,ip,L,C)
%     f_handle = @(x) arc_length_ribe(x,theta*pi/180,L);
%     x_pb = x_p(ip==1);
%     X = unique(x_pb);
%     for i = 1:length(X)
%         Buf(i) = integral(f_handle,C,X(i));  
%     end
% end
% 
% 
