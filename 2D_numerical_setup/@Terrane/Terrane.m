classdef Terrane
    properties
        Boundary                    = [0,0,0,0];  %x1 x2 y1 y2 
        Phases 
        Stratigraphy
        Thermal_information 
        Thermal_type  
    end
    methods
        function [Phase,Temp] = fill_terranes(obj,A,Phase,Temp)
            [Phase,Temp] = fill_layer(obj,A,Phase,Temp);
        end

    end

end

