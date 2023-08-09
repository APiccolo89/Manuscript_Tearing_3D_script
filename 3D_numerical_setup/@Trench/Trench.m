classdef Trench

    properties
        name = 'Subduction zone'
        Type_Subduction % Subduction type {Mode_1}{Ribe}{To_Do: Linear _World_Builder_example => CITE IT}
        Boundary = BoundsT; 
                   %================================================================================
        Boundaries_list = {'none'} %where to apply the subduction zone.  
        % Thermal information +/- Generic information needed
        Thermal_type        % typology thermal structure
        Thermal_information % Data required thermal structure
        % Information concerning the subduction
        R       % Curvature radius
        C       % Center of radius curvature
        theta  = {[90, 90], 'none'} % Bending angle
        tk_WZ   % thickness of the weak zone
        ph_WZ   % phase of the weak zone
        L0      % Lenght of the slab
        D0      % Thickness of the slab
        Decoupling_depth % Decoupling depth
        % Stratigraphic information
        Stratigraphy_Continental %
        Stratigraphy_Oceanic     %
        % Subducted crust
        Subducted_crust_L        % Polygon crust subduction
        % Orogenic prism information
        position_end_prism= 100;
        phase_prism
        % Field that are used internally to the class to do the proper
        % computation
        d_slab % store the distance from the top surface
        l_slab % store the length of the slab {useful for McKenzie temperature profile}
        length_continent = {[100,20],'none'} 
        continent 
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
        [Phase,Temp] = fill_subduction(obj,A,Phase,Temp,Boundary,theta);
    end
    methods (Access = private)
        % Function that modify the object only when is needed and without
        % consuming memory {everything is locally used}
        [obj,Phase,Temp] = find_slab_(obj,A,Weak_Slab,Phase,Temp,Boundary,theta)
        % Slab mode 1
        [obj,Phase,Temp]= find_slab_mode_1(obj,A,Weak_Slab,Phase,Temp,Boundary,theta)
        [l,ind_decoupling,Phase,Temp] = find_length_slab(obj,x,z,C,r,d,l,Phase,Temp) % I need the data of the object for this particular function
        [Phase] = generate_accretion_prism(obj,A,Phase)
        function [Phase] = fill_weak_zone(obj,Phase)
            ind = ~isnan(obj.d_slab(:))& Phase(:)~= obj.phase_prism;
            Phase(ind) = obj.ph_WZ;

        end
    end
end


