function [Layout,arc_angleS,continent] = find_slab_(obj,A)
ang_  = (obj.theta)*pi/180;
ang_2 = (obj.theta+90)*pi/180;
r      = [obj.R T.R+obj.D0];
sl    = 1               ;
arc_angleS = [];
if ~strcmp(obj.Type_Subduction,'Mode_1')
    ind = (A.Xpart(:)>=obj.C(1) & A.Xpart(:)<=obj.C(1)+2.*obj.L0) & A.Zpart(:)>=-1.5*obj.L0 & A.Zpart(:)<=1.0;
    A_time = cputime;
    x     = A.Xpart(:);
    z     = A.Zpart(:);
    y     = A.Ypart(:);
    % Vector version
    Layout = x*nan;
    d      = Layout;
    arc_angleS = x.*NaN;
    continent = arc_angleS;
    angle_c = continent;
    u = zeros(2,length(x));
    v = u;
    u(1,ind==1) = x(ind==1)';
    u(2,ind==1) = z(ind==1)';
    d(ind==1) = sqrt((u(1,ind==1)-obj.C(1)).^2+(u(2,ind==1)- obj.C(2)).^2);
    arc_angleS(ind==1) = acos((z(ind==1)-obj.C(2))./d(ind==1)).*180/pi;
    c_slab = d(ind==1)>=r(1) & d(ind==1)<=(r(2)) & arc_angleS(ind==1)<=theta & arc_angleS(ind==1)>=0;
    angle_c(ind==1) = d(ind==1)>=r(2)-obj.Cont+((obj.Cont-0.0)./(theta_c)).*arc_angleS(ind==1) & d(ind==1)<=(r(2)); %& arc_angleS(ind==1)<=theta_c ;
    continent(angle_c==1)=1.0;
    Layout= d;
    Layout(d<r(1) | d>r(2)) = NaN;
    Layout(arc_angleS>theta | arc_angleS<0)=NaN;
    Layout = Layout-r(2);
    if strcmp(type,'Weak')
        r(1) = r(2);
        r(2) = r(2)+obj.tk_WZ;
        if obj.theta == 90
            sl = 0;
        end
        %=============

    end

    obj.C(1)     =obj.C(1);
    obj.C(2)     = obj.obj.C(2);
    p1(1)    = obj.C(1)+r(1)*sin(ang_);
    p1(2)    = obj.C(2)+r(1)*cos(ang_);
    p2(1) = obj.C(1)+r(2)*sin(ang_);
    p2(2) = obj.C(2)+r(2)*cos(ang_);
    p3(1)    = p1(1)+obj.L0*sin(ang_2);
    p3(2)    = p1(2)+obj.L0*cos(ang_2);
    p4(1)    = p2(1)+obj.L0*sin(ang_2);
    p4(2)     = p2(2)+obj.L0*cos(ang_2);
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

if strcmp(obj.Type_Subduction,'Ribe_Mode')
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
