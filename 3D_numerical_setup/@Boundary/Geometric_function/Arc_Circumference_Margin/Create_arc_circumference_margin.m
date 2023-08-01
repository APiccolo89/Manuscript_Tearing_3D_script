function [obj] = Create_arc_circumference_margin(obj,R,Boundary)

if strcmp(Boundary,'A') || strcmp(Boundary,'C')
    if isempty(obj.L)
      error('You have to specify L in the class before using the function for this boundary')
    end
    [obj.W,xc,yc]= obj.modify_boundary_limits(R,arc_length,Boundary);

elseif strcmp(Boundary,'B') || strcmp(Boundary,'D')
    if isempty(obj.W)
      error('You have to specify W in the class before using the function for this boundary')
    end
    [obj.L,xc,yc]= obj.modify_boundary_limits(R,arc_length,Boundary);

else
    error('The margin is wrong A C are perpendicular to axis x, C D are perpendicular to y')
end
% Create all the other data. 
obj.compute_coordinate_boundary_composite
% 
if strcmp(Boundary,'A')
elseif strcmp(Boundary,'B')
elseif strcmp(Boundary,'C')
elseif strcmp(Boundary,'D')
    obj.D{3} = @(x) circumference_margin(obj,xc,yc,R,x); %function handle boundary
    obj.D{4} = @(x,y) arc_length(x,y,xc,yc,R);           %function handle boundary 
end



end
