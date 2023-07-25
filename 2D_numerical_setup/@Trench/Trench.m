classdef Trench

    properties
        Type_Subduction % Subduction type {Mode_1}{Ribe}{To_Do: Linear _World_Builder_example => CITE IT}
        Boundary = [0,0]% x-y coordinate trench 
        R       % Curvature radius
        C       % Center of radius curvature
        theta   % Bending angle
        theta_c % angle continental subduction
        tk_WZ   % thickness of the
        L0      % Lenght of the slab
        D0      % Thickness of the slab
        Decoupling_depth % Decoupling depth
        Stratigraphy_Continental % 
        Stratigraphy_Oceanic     %
        Thermal_type        % typology thermal structure
        Thermal_information % Data required thermal structure
        Layout % store the distance from the top surface
        Lenght % store the length of the slab {useful for McKenzie temperature profile}
        continent %highlight where is the continental crust 
    end
    methods
        function [Phase,Temp] = fill_terranes(obj,A,Phase,Temp)
            [Phase,Temp] = fill_subduction(obj,A,Phase,Temp);
        end
    end
end


