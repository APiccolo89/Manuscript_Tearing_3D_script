function [y] = circumference_margin(xa,xb,xc,ya,R,x)

% Check if the user had the brilliant idea to create a circumference wh
% radius is less the the distance between the two points
if R<(xb-xa)/2
    error('Just think to the pytaghora theorem, i.e., you cannot have a triangle rectangle whose hyp is less than one of the other two')
end
% I try to put more general formulation, but the stuff for this particular
% case is easy
c_ = (xa-xc).^2-R.^2+ya^2;
b_ = 2.*ya; 
delta = sqrt(b_^2-4.*c_);
center_y = [(b_-delta)./2;(b_+delta)./2]; 
y_ind = find(center_y<=ya,1); 
y_c = center_y(y_ind);
% Compute the semi-circumference  (only positive)
y = (R^2-(x-Cx).^2).^(0.5)+y_c; 
end

