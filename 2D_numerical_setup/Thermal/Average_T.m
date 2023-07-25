function [Temp] = Average_T(obj,D,ind,Temp)
%AVERAGE_T Fill up the layer with a constant temperature
%   Detailed explanation goes here
T_tk = [0.0, obj.Stratigraphy.Tk(end)]; % Stratigraphy of the object

Tsl   = obj.Thermal_type.AverageT; % Average Temperature in C 

ind_z = (D<T_tk(1) & D>=T_tk(2) & ind == 1);

Temp(ind_z==1)=Tsl;
end

