function [Temp]= compute_temperature_profile(obj,A,ind,Temp)
%=========================================================================
% This is a public function that is accessed by several other classes. 
% obj: -> Terrane/Boundary that needs to be filled with temperature 
% D  : -> Depth {Z coordinates or distance from the top surface in case of subduction}
% A  : -> Coordinate of the particles, 
% ind : where to fill the temperature field. 
% Temp: Temperature array
%==========================================================================
% OUTPUT: -> Temp updated.
%==========================================================================

Type = obj.Thermal_type.Type; % type of continental geotherm
if isa((obj),'Terrane')
    D = A.Zpart;  % For Layered structure Z coordinate make its job
elseif isa((obj),'Trench')
    ind = ~isnan(obj.d_slab); 
    D = obj.d_slab; % Layout is the member of the trench that collects all the distance from the top surface of the slab
end

if strcmp(Type,'ContinentalGeotherm') % Is continental
    [Temp] = Continental_Geotherm(obj,D,ind,Temp);

elseif strcmp(Type,'HalfSpace')
    [Temp] = HalfSpaceCooling(obj,D,ind,Temp);
elseif strcmp(Type,'McKenzie')
    % McKenzie type require to run the temperature field of Halfspace
    % cooling model before. As the quick hack to have a smooth temperature
    % field is to introduce a weighted average as a function of the
    % decoupling depth. 
    [Temp] = HalfSpaceCooling(obj,D,ind,Temp);
    [Temp] = compute_temperature_profile_McKenzie(obj,Temp,A);
elseif strcmp(Type,'Average_T')
    [Temp] = Average_T(obj,D,ind,Temp);
else
    error('You do not provide a temperature type:{Continental},{HalfSpace},{McKenzie},{Average_T}, please correct your ways and repent yourself.')
end
end
