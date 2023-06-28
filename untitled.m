% Ribe 
clf;
clear all; 
close all; 
L = 500; 
x_t = 0:0.1:500;
teta_0= 10:1:90; 
[X,T] = meshgrid(x_t,teta_0);
T     = T.*pi./180;
teta=T.*X.^2.*((3.*L-2.*X))./(L^3);
z =0; 
teta_0 = 90.*pi./180;
teta_s = teta_0.*x_t.^2.*((3.*L-2.*x_t))./(L^3);
x_ = sin(teta_s.*((x_t-x_t(1))./(x_t(end)-x_t(1))));
z_ = cos(teta_s.*((x_t-x_t(1))./(x_t(end)-x_t(1))));


P = [0.0,-120]; 
theta_0 = 30; 
z_t = zeros(length(x_t),1); 
d = sqrt((x_t-P(1)).^2+(z_t'- P(2)).^2);
arc_angles = acos((z_t'-P(2))./d).*180/pi;
theta_s    = arc_angles;
theta_s(theta_s>30)=(90-theta_0);
x_correction = ones(length(teta_s),1).*0;
z_correction = -ones(length(teta_s),1).*0;
A = x_t(theta_s>30);
x_correction(theta_s>30) = (2*120*sin(theta_0./2*pi/180))*cos(theta_0*pi/180);%2*120.*(sin((theta_0./2)*pi/180)*cos(theta_0./2*pi/180));
z_correction(theta_s>30) = -2*120.*(sin((theta_0./2)*pi/180)*sin(theta_0*pi/180));


teta_s = theta_s*pi./180; 



for i=1:length(x_t)
    PA = [x_t(i)-x_correction(i);z_t(i)-z_correction(i)];
    PC = [x_t(i);z_t(i)];
    rotation = [cos(teta_s(i)) sin(teta_s(i)); -sin(teta_s(i)) cos(teta_s(i))];
    PB_ = rotation*PA;
    PD_ = rotation*PC; 
    PB(i,1) = PB_(1)+x_correction(i);%*tan(30*pi/180);
    PB(i,2) = PB_(2)+z_correction(i); 
    PD(i,1) = PD_(1);
    PD(i,2) = PD_(2); 

end

y = 

figure(1)
hold on; 
scatter(PB(:,1),PB(:,2))
plot(x_t,z_t);
bla = 0; 



function [d] = find_distance_linear(P,p1,p2)

% Ok, wiki saves my day because I was too lazy:
A  = abs((p2(1)-p1(1)).*(p1(2)-P(:,2))-(p1(1)-P(:,1)).*(p2(2)-p1(2)));
B  = (p2(1)-p1(1)).^2+(p2(2)-p1(2)).^2;
d  = A./sqrt(B);

end

