function [Temp]= compute_temperature_profile(obj,A,ind,Boundary)
%=========================================================================
% This is a public function that is accessed by several other classes. 
% obj: -> Terrane/Boundary that needs to be filled with temperature 
% D  : -> Depth {Z coordinates or distance from the top surface in case of subduction}
% A  : -> Coordinate of the particles, 
%    : This variable changes its meaning. A can be the main structure where
%    all the properties belong, or B the transformed coordinate structure,
%    a temporary variable that is destroyed and create where it is needed.
%    
% ind : where to fill the temperature field. % Optional in many function,
% in case of slab or weak zone. obj.d_slab is already containing the
% information where are the point. 
% Temp: Temperature array
%==========================================================================
% OUTPUT: -> Temp updated
%==========================================================================

Type = obj.Thermal_type.Type; % type of continental geotherm
if isa((obj),'Terrane')
    % If the object is a simple terrane, the main stratigraphication is a
    % planar stratification. In the future could be cool introduce bending
    % and other feature but for now, it is ok this state of art. 
    D = A.Zpart;  % For Layered structure Z coordinate make its job
elseif isa((obj),'Trench')
    % If the object is Trench, the actual depth is given by the
    % perpendicalr distance from the top surface towards the interior of
    % the slab. 
    ind = ~isnan(obj.d_slab); 
    D = obj.d_slab; % Layout is the member of the trench that collects all the distance from the top surface of the slab
end

if strcmp(Type,'ContinentalGeotherm') % Is continental
    % Compute continental Geoterhm
    [Temp] = Continental_Geotherm(obj,D,ind,A.Temp);

elseif strcmp(Type,'HalfSpace')
    [Temp] = HalfSpaceCooling(obj,D,ind,A.Temp);
elseif strcmp(Type,'Ridge')
    [Temp] = ridge_oceanic(obj,D,A.Xpart,A.Ypart,ind,A.Temp);
elseif strcmp(Type,'McKenzie')
    % McKenzie type require to run the temperature field of Halfspace
    % cooling model before. As the quick hack to have a smooth temperature
    % field is to introduce a weighted average as a function of the
    % decoupling depth. 
    [A.Temp] = HalfSpaceCooling(obj,D,ind,A.Temp);
    [Temp] = compute_temperature_profile_McKenzie(obj,A,A.Temp,Boundary);
elseif strcmp(Type,'Average_T')
    [Temp] = Average_T(obj,D,ind,A.Temp);
else
    error('You do not provide a temperature type:{Continental},{HalfSpace},{McKenzie},{Average_T}, please correct your ways and repent yourself.')
end
end
