
function [y]  = linear_margin(xa,theta,y0,x)

% Function to describe the linear margin of a terranes 
theta_r = theta*pi/180;
m       = tan(theta_r);
y       = m.*(x-xa)+y0; 

end


function [theta_c]  = continental_theta_angle(x,xa,xb,theta_c0,theta_c1)

% Function to describe the linear margin of a terranes 
theta_c = theta_c0+(((theta_c1-theta_c0))./(xb-xa)).*(x-xa);

end

function [y]  = parabolic_margin(xa,xv,yv,x)

% Function to describe the linear margin of a terranes 
a = -yv./(xa-xv).^2;
y = a.*(x-xv).^2+yv; 

end

function [y] = circumference_margin(xa,xb,ya,Cx,R,x)

% Check if the user had the brilliant idea to create a circumference whose
% radius is less the the distance between the two points
if R<(xb-xa)/2
    error('Just think to the pytaghora theorem, i.e., you cannot have a triangle rectangle whose hyp is less than one of the other two')
end
% I try to put more general formulation, but the stuff for this particular
% case is easy
c_ = (xa-Cx).^2-R.^2+ya^2;
b_ = 2.*ya; 
delta = sqrt(b_^2-4.*c_);
center_y = [(b_-delta)./2;(b_+delta)./2]; 
y_ind = find(center_y<=ya,1); 
y_c = center_y(y_ind);
% Compute the semi-circumference  (only positive)
y = (R^2-(x-Cx).^2).^(0.5)+y_c; 
end

function [xa,xb,xc,yc] = find_terrane_size(R,arc_length)
theta_ = arc_length./(R); 
xa     = -R.*sin(theta_./2);
xb     = R.*sin(theta_./2);
xc     = 0.0; 
yc     = -R.*cos(theta_./2); 
end

function [s] = arc_length_x(x,y,xa,ya,R)
% compute the arc length of a given point using the Pa - Px points 
d = sqrt((x-xa).^2+(y-ya).^2);
theta = 2*asin(d./(2*R));
s = R.*theta; 
end


function [W] = length_parabolic_margin(xa,xv,yv,x)
a = -yv./(xa-xv).^2;
xb  = -xa; 
t  = (xa-xv);
t2 = (xb-xv); 
W1 = (2.*a.*t.*(4.*a.*a.*t.^2+1)+asinh(2.*a.*t))./(4.*a);
W2 = (2.*a.*t2.*(4.*a.*a.*t2.^2+1)+asinh(2.*a.*t2))./(4.*a);
W = W2-W1;
end
