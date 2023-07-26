classdef Trench

    properties
        name = 'Subduction zone'
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
        d_slab % store the distance from the top surface
        l_slab % store the length of the slab {useful for McKenzie temperature profile}
        continent %highlight where is the continental crust 
    end
    methods (Access = public)
        % Function that interacts with the external enviroment: Take the
        % Phase and Temp array and modify accordingly

        function [Phase,Temp] = fill_terranes(obj,A,Phase,Temp)
            disp([obj.name, 'object is filling....'])
            A_time = cputime; 
            [Phase,Temp] = fill_subduction(obj,A,Phase,Temp);
            B_time = cputime; 
            disp(['and took ', num2str(B_time-A_time,3), 'seconds']);
        end
        [Phase,Temp] = fill_subduction(obj,A,Phase,Temp);
    end
    methods (Access = private)
      % Function that modify the object only when is needed and without
      % consuming memory {everything is locally used}
      [obj,Phase,Temp] = find_slab_(obj,A,Weak_Slab,Phase,Temp)
      % Slab mode 1
      [obj,Phase,Temp]= find_slab_mode_1(obj,A,Weak_Slab,Phase,Temp)
      [l,ind_decoupling,Phase,Temp] = find_length_slab(obj,x,z,C,r,d,l,Phase,Temp) % I need the data of the object for this particular function 
   end
end


