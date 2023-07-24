classdef Passive_Margin
    %PASSIVE_MARGIN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
                Boundary       = []; 
                Boundary_C     = []; %Boundary connected continent
                Stratigraphy         %Stratigraphy bounded continent   
                Phases               %sedimentary phase 
                depo_center     = 12;
                center_pos = 0.2; 
                shape         = 'triangular'; %'rectangular'
                length          = 200; 
                d_lithos        = 20; 
                Age             = 80; 
                Thermal_type   % Thermal type of the attached continent        
    end
    
    methods
        
    end
end

