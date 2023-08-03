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
        function [Phase,Temp] = fill_terranes(obj,A,Phase,Temp)
            disp([obj.name, 'object is filling....'])
            A_time = cputime;
            [Phase,Temp] = fill_layer(obj,A,Phase,Temp);
            if ~isempty(obj.Passive_Margin.Direction)
                disp(['      Passive_Margins are filling....'])
                [Phase,Temp]= obj.Passive_Margin.generate_passive_margin(obj,A,Phase,Temp);
                P_time = cputime;
                disp(['         and took ', num2str(P_time-A_time,3), 'seconds']);
            end
            B_time = cputime;
            disp(['and took ', num2str(B_time-A_time,3), 'seconds']);
        end

    end

end

