function [s] = arc_length(x,y,xa,ya,R)
% compute the arc length of a given point using the Pa - Px points 
d = sqrt((x-xa).^2+(y-ya).^2);
theta = 2*asin(d./(2*R));
s = R.*theta; 
end