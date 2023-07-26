classdef Terrane
    properties
        name = 'Nameless Terrane'
        Boundary                    = [0,0,0,0];  %x1 x2 y1 y2 
        Phases 
        Stratigraphy
        Thermal_information 
        Thermal_type  
    end
    methods
        function [Phase,Temp] = fill_terranes(obj,A,Phase,Temp)
            disp([obj.name, 'object is filling....'])
            A_time = cputime; 
            [Phase,Temp] = fill_layer(obj,A,Phase,Temp);
            B_time = cputime; 
            disp(['and took ', num2str(B_time-A_time,3), 'seconds']);
        end

    end

end

