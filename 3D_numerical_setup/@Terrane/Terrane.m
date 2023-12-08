classdef Terrane
    properties
        name = 'Nameless Terrane' % Name Terrane
        Boundary                    = BoundsT %
        Phases
        Stratigraphy
        Thermal_information = Thermal_Information; % Default value of class
        Thermal_type   = Thermal_Type;             % Default value of class
        Passive_Margin  = Passive_Margin;          % Default value of class
    end
    methods
        function [A] = fill_terranes(obj,A)
            disp([obj.name, 'object is filling....'])
            A_time = cputime;
            [A] = fill_layer(obj,A);
            if ~isempty(obj.Passive_Margin.Direction)
                disp(['      Passive_Margins are filling....'])
                [A]= obj.Passive_Margin.generate_passive_margin(obj,A);
                P_time = cputime;
                disp(['         and took ', num2str(P_time-A_time,3), 'seconds']);
            end
            B_time = cputime;
            disp(['and took ', num2str(B_time-A_time,3), 'seconds']);
        end

    end

end

