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
        Depth_continent = -50; 
        Stratigraphy_Continental % 
        Stratigraphy_Oceanic     %
        Thermal_type        % typology thermal structure
        Thermal_information % Data required thermal structure
        Layout % store the distance from the top surface
        Lenght % store the length of the slab {useful for McKenzie temperature profile}
        continent %highlight where is the continental crust 
    end
    methods (Access = public)
        % Function that interacts with the external enviroment: Take the
        % Phase and Temp array and modify accordingly

        function [Phase,Temp] = fill_terranes(obj,A,Phase,Temp)
            [Phase,Temp] = fill_subduction(obj,A,Phase,Temp);
        end
        [Phase,Temp] = fill_subduction(obj,A,Phase,Temp);
    end
    methods (Access = private)
      % Function that modify the object only when is needed and without
      % consuming memory {everything is locally used}
      obj = find_slab_(obj,A,Weak_Slab)
      obj = find_slab_mode_1(obj,A,Weak_Slab)
   end
end


