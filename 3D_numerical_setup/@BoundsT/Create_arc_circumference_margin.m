function [obj] = Create_arc_circumference_margin(obj,R,Boundary,arcLength)
%=========================================================================
% Function that creates the necessary data for a terrane whose boundary. 
% Input: R the radius of curvature of the circular portion of the terrane
% Boundary: A,B,C,D => the boundary where to apply. 
% arcLength, the length of the portion NB: in order to comply to that the
% main 4 points defining the terrane might be modified, so must be done
% consistently. The current implementation serves only the goal of this
% project. The enviroment must be improved.
%==========================================================================


if strcmp(Boundary,'A') || strcmp(Boundary,'C')
    if isempty(obj.L)
        error('You have to specify L in the class before using the function for this boundary')
    end
    [obj.W,CC]= obj.modify_terrane_limits(R,arc_length,Boundary);

elseif strcmp(Boundary,'B') || strcmp(Boundary,'D')
    if isempty(obj.W) || obj.W==1.0 
        error('You have to specify W in the class before using the function for this boundary')
    end
    obj=obj.compute_coordinate_boundary_composite;
    [obj.L,CC]= obj.modify_terrane_limits(R,arcLength,Boundary,obj.c(1));

else
    error('The margin is wrong A C are perpendicular to axis x, C D are perpendicular to y')
end
% Create all the other data.
obj=obj.compute_coordinate_boundary_composite;
% This part of the code should be cleansed whenever someone has time. The
% goals justify the means, and I would like to finish the 3D project by the
% end of January with already part of the manuscript written
if strcmp(Boundary,'A')
elseif strcmp(Boundary,'B')
    obj.B{2} = 'Circular'; %function handle boundary
    obj.B{3} = [CC,R]  ; % save the important information.  
elseif strcmp(Boundary,'C')
elseif strcmp(Boundary,'D')
   
end
end
