classdef Thermal_Type
    properties
        Type = 'HalfSpaceCooling' % 'ContinentalGeotherm', 'Average_Temperature', 'McKenzie'
        Age  = 100;  %Default Value Age
        vel  = {[1,1],'none'}  ; % Default Value velocity
        Moho = 500 ; % Moho Temperaure
        Moho_d = -35; % Default valu Moho. 
    end
    
end

