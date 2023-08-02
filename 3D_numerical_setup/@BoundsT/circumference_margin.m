function [y] = circumference_margin(obj,Boundary,x)
%
% For being clear: the curvature IS ALWAYS directed towards the external
% portion of the terrain NOT towards its interiour. 
%
%
%
% I do not have clever idea on how to do in a simplified and elegant way. 
if strcmp(Boundary,'A')
elseif strcmp(Boundary,'B')
    B=obj.B{3}; %function handle boundary
    R   = B(3);
    coord  = obj.B{1};
    xa = coord(1);
    ya = coord(4);
    xb = coord(3);
    c = B(1); 
    sign = +1; 
elseif strcmp(Boundary,'C')
elseif strcmp(Boundary,'D')

end
% Check if the user had the brilliant idea to create a circumference wh
% radius is less the the distance between the two points
if R<(xb-xa)/2
    error('Just think to the pytaghora theorem, i.e., you cannot have a triangle rectangle whose hyp is less than one of the other two')
end
% I try to put more general formulation, but the stuff for this particular
% case is easy
c_ = (xa-c).^2-R.^2+ya^2;
b_ = 2.*ya; 
delta = sqrt(b_^2-4.*c_);
center_y = [(b_-delta)./2;(b_+delta)./2]; 
y_ind = find(center_y<=ya,1); 
y_c = center_y(y_ind);
% Compute the semi-circumference  (only positive)
y = (R^2-(x-c).^2).^(0.5)+y_c; 
end

