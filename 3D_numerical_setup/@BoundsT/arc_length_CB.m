function [s] = arc_length_CB(obj,x,y,Boundary)
% compute the arc length of a given point using the Pa - Px points 
if strcmp(Boundary,'A')
elseif strcmp(Boundary,'B')
    B=obj.D{3}; %function handle boundary
    R   = B(3);
    CC = B(1:2); 
elseif strcmp(Boundary,'C')
elseif strcmp(Boundary,'D')
   
end
d = sqrt((x-CC(1)).^2+(y-CC(2)).^2);
theta = 2*asin(d./(2*R));
s = R.*theta; 
end