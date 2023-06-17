% Ribe 
clf;

L = 500; 
x_t = 0:1:500;
teta_0= 10:1:90; 
[X,T] = meshgrid(x_t,teta_0);
T     = T.*pi./180;
teta=T.*X.^2.*((3.*L-2.*X))./(L^3);
z =0; 
teta_0 = 45.*pi./180;
teta_s = teta_0.*x_t.^2.*((3.*L-2.*x_t))./(L^3);
x_ = sin(teta_s.*((x_t-x_t(1))./(x_t(end)-x_t(1))));
z_ = cos(teta_s.*((x_t-x_t(1))./(x_t(end)-x_t(1))));

for i=1:length(x_t)
    PA = [x_t(i);0];
    PC = [x_t(i);-80];
    rotation = [cos(teta_s(i)) sin(teta_s(i)); -sin(teta_s(i)) cos(teta_s(i))];
    PB_ = rotation*PA;
    PD_ = rotation*PC; 
    PB(i,1) = PB_(1);
    PB(i,2) = PB_(2); 
    PD(i,1) = PD_(1);
    PD(i,2) = PD_(2); 

end


plot(PD(:,1),PD(:,2))