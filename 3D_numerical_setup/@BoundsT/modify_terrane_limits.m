function [Size,CC] = modify_terrane_limits(obj,R,arc_length,Boundary,c)
%=========================================================================%
% Argument Radius of curvature
% R = radius of curvature
% arc_length = the length of the curved slab
%=========================================================================%
<<<<<<< HEAD
theta_ = arc_length/(R);
if theta_ >= 270/pi
    error('For the given terrane, this curvature is irrealistic')
=======

theta_ = arc_length./(R);
if theta_ > pi 
    error('This algoritm is working, but there are some restriction. The maximum curvature allowed is arc_length/r=pi')
>>>>>>> Python_3D_post_process
end
xa     = c-R.*sin(theta_./2);
xb     = c+R.*sin(theta_./2);
Size = xb-xa;
if strcmp(Boundary,'A')||strcmp(Boundary,'C')
    xc     = -R.*cos(theta_./2);
    yc     = obj.c(2);

else
    xc     = obj.c(1);
    if strcmp(Boundary,'B')
        yc     = obj.y2+R.*cos(theta_./2);
    else
        yc     = obj.y1-R.*cos(theta_./2);
    end
end
CC = [xc,yc]; % Center of curvature boundary
end
