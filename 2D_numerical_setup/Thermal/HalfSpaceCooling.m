function [Temp] = HalfSpaceCooling(obj,D,ind,Temp)
%HALFSPACECOOLING Function that compute the half space cooling 
%INPUT: 
% obj -> Geological object
% D   -> Depth {Z coordinates,distance from the surface}
% Temp -> Temperature array
% ind -> What are the choosen particles
% ==== To Do in the future ==== 
% Ridge -> Position {To do in the future}
% vel_spread -> Velocity at which the ocean spread {To do in the future}
% =============================
% OUTPUT:updated temperature field. 

T_tk = [0.0, obj.Stratigraphy.Tk(end)]; % Stratigraphy of the object
kappa  = obj.kappa; %Thermal conductivity of the object 
T_Age = obj.Age;    % Age of the object in [s]
T_P   = obj.Thermal_information.TP; % Mantle potential temperature in C
T_S   = obj.Thermal_information.TS; % Surface temperature in C

ind_z = D < T_tk(1) & D>=T_tk(2) & ind==1;
ind_o = D < T_tk(2) & ind==1;
    
erf_function = erf(D(ind_z).*1000/2/(kappa*T_Age)^0.5);
Temp(ind_z) = T_S - (T_P-T_S).*erf_function;
Temp(ind_o) = T_P;
end

