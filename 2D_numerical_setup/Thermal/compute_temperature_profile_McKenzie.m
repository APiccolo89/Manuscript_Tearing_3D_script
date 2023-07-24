
function [Temp]= compute_temperature_profile_McKenzie(Layout,Temp,Type,Gen,Terranes,indx,A)

T_tk = [0.0, Terranes.Stratigraphy(end)];
k = Terranes.K;
T_Age = Terranes.Age.*Terranes.secMyrsyear;
T_P   = Gen.T_P;
T_S   = Gen.T_S;
Tsl   = Gen.AvTS;
vl    = Terranes.Trench_properties.vl; 
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
    % Decoupling depth-length
    i = find(zprojection>=Terranes.Trench_properties.Decoupling,1);  
    weight(slab==1) = 0+(0.9-0.1)./(l(i)-0.1).*(l(slab==1)-0.1);
    weight(weight>1) =1.0; 
    weight(weight<0) =0.01; 

    weight = weight'; 
    Temp(slab==1)=(weight(slab==1)).*t_prov(slab==1)+(1-weight(slab==1)).*Temp(slab==1);
    Temp = reshape(Temp,size(A.Xpart)); 
end



