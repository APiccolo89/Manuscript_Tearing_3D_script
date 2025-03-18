function [s] = compute_arc_length_circle(obj,x,Boundary)
[y] = obj.circumference_margin(Boundary,x);
%COMPUTE_ARC_LENGTH_CIRCLE Summary of this function goes here
%   Detailed explanation goes here compute the arc length of a given point using the Pa - Px points 
if strcmp(Boundary,'A')
elseif strcmp(Boundary,'B')
    B=obj.B{3}; %function handle boundary
    R   = B(3);
    xa  = obj.B{1}(1);
    ya  = obj.B{1}(2);
elseif strcmp(Boundary,'C')
elseif strcmp(Boundary,'D')
   
end
d = sqrt((x-xa).^2+(y-ya).^2);
theta = 2*asin(d./(2*R));
s = R.*theta; 
end

