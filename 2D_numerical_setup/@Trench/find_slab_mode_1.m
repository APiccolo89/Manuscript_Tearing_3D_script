function [obj] = find_slab_mode_1(obj,A,Weak_Slab)
% Convert the angle in radians
ang_  = (obj.theta)*pi/180;
ang_2 = (obj.theta+90)*pi/180; %needed for the linear portion of the slab
r      = [obj.R-obj.D0,obj.R]; % radius upper and lower surface
r_m    = sum(obj.R)./2; 
C      = [obj.Boundary,obj.Boundary-obj.R];
sl    = 1    ;
arc_angleS = [];
% select the point that are worth to check if they belong to the slab
ind = (A.Xpart(:)>=C(1) & A.Xpart(:)<=C(1)+2.*obj.L0) & A.Zpart(:)>=-1.5*obj.L0 & A.Zpart(:)<=1.0;
A_time = cputime;
x     = A.Xpart(:);
z     = A.Zpart(:);
y     = A.Ypart(:);
% Vector version
Layout = x*nan;
Length = Layout; 
d      = Layout;
arc_angleS = x.*NaN;
continent = arc_angleS;
angle_c = continent;
u = zeros(2,length(x));
u(1,ind==1) = x(ind==1)';
u(2,ind==1) = z(ind==1)';
d(ind==1) = sqrt((u(1,ind==1)-C(1)).^2+(u(2,ind==1)- C(2)).^2);
arc_angleS(ind==1) = acos((z(ind==1)-C(2))./d(ind==1)).*180/pi;
angle_c(ind==1) = d(ind==1)>=r(2)-obj.Depth_Continet+((obj.Depth_Continet-0.0)./(theta_c)).*arc_angleS(ind==1) & d(ind==1)<=(r(2)); %& arc_angleS(ind==1)<=theta_c ;
continent(angle_c==1)=1.0;
Layout= d;
d = []; 
Layout(d<r(1) | d>r(2)) = NaN;
Layout(arc_angleS>obj.theta | arc_angleS<0)=NaN;
Layout = Layout-obj.R(2);
% Compute the length of the curved slab as a function of the mid distance
% layer. 
Length(~isnan(Layout))= r_m.*arc_angleS(slab==1 & angle<=theta)*pi/180; % compute the length associated with the 


if strcmp(Weak_Slab,'Weak')
    r(1) = r(2);
    r(2) = r(2)+obj.tk_WZ;
    if obj.theta == 90
        sl = 0;
    end
    %=============

end


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
Pz(1) = p1(2);
Pz(2) = p2(2);
Pz(3) = p4(2);
Pz(4) = p3(2);
% Inpolygon
[in,~] = inpolygon(x,z,Py,Pz);
P = [x(in), z(in)];
Layout(in) = -(find_distance_linear(P,p2,p4))';
B_time = cputime;
disp(['Time Loop = ', num2str(B_time-A_time)]);
obj.Length = compute_length_slab();
obj.Layout = reshape(Layout,size(A.Xpart));
obj.continent = reshape(continent,size(A.Xpart));
end


function [d] = find_distance_linear(P,p1,p2)

% Ok, wiki saves my day because I was too lazy:
A  = abs((p2(1)-p1(1)).*(p1(2)-P(:,2))-(p1(1)-P(:,1)).*(p2(2)-p1(2)));
B  = (p2(1)-p1(1)).^2+(p2(2)-p1(2)).^2;
d  = A./sqrt(B);

end
