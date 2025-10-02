classdef Thermal_Type
    properties
        Type = 'HalfSpaceCooling' % 'ContinentalGeotherm', 'Average_Temperature', 'McKenzie','Ridge'
        Age  = 100;  %Default Value Age
        vel  = 0 ; % Default Value velocity
        Moho = 500 ; % Moho Temperaure
        Moho_d = -35; % Default valu Moho. 
        Ridge_position = []
        v_spreading = 1.0 
        
    end
    
end

