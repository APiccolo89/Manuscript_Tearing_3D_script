classdef Trench

    properties
        Boundary = [0,0]% x-y coordinate trench 
        R       % Curvature radius
        C       % Center of radius curvature
        theta   % Bending angle
        theta_c % angle continental subduction
        tk_WZ   % thickness of the
        L0      % Lenght of the slab
        D0      % Thickness of the slab
        Decoupling_depth % Decoupling depth
        tk_Cont % Continental thickness
        Stratigraphy_Continental % 
        Stratigraphy_Oceanic     %
        Type_Subduction
        Thermal_type        % typology thermal structure
        Thermal_information % Data required thermal structure
        Layout
        arc_angleS
        continent
    end
    methods
        function [Phase,Temp] = Create_Subduction_zone(obj,A,Phase,Temp,Gen)
            [Phase,Temp] = fill_subduction(obj,A,Phase,Temp,Gen);
        end
    end
end


