classdef Trench
%=========================================================================%
% Definitions, and use of this class
%=========================================================================%
    properties
        name = 'Subduction zone'

        Type_Subduction % Subduction type {Mode_1}{Mode_2}{To_Do: Linear _World_Builder_example => CITE IT} => CITE!

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

        Lb      % Bending lenght [where to apply variable bending angle]

        D0      % Thickness of the slab

        Decoupling_depth % Decoupling depth
        % Stratigraphic information

        Stratigraphy_Continental %

        Stratigraphy_Oceanic     %
        % Subducted crust

        Subducted_crust_L        % Polygon crust subduction
        % Field that are used internally to the class to do the proper
        % computation

        d_slab % store the distance from the top surface

        l_slab % store the length of the slab {useful for McKenzie temperature profile}

        length_continent = {[100,20],'none'} 

        continent 
        % Mode 2 data 
        Type_Angle % type of angle {linear, ribe}

        Slab_surface % Data structure containing top, mid and bottom surface

        nseg         % Resolution of the slab
        % Orogenic prism information

        prism_depth = -70; 

        position_end_prism= 100;

        phase_prism   

        Prism_lc_depth = -40;

        C_prism % Center of prism {if empty, it takes the center of the curvature of the slab}


    end

    methods (Access = public)
        % Function that interacts with the external enviroment: Take the
        % Phase and Temp array and modify accordingly

        function [A] = fill_terranes(obj,A)
            disp([obj.name, 'object is filling....'])
            A_time = cputime;
            [A] = fill_subduction(obj,A);
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
        [obj,Phase,Temp]= find_slab_mode_2(obj,A,Weak_Slab,Phase,Temp,Boundary,theta)
        [obj] = Compute_top_bottom_surface(obj,theta)
        [l,ind_decoupling,Phase,Temp] = find_length_slab(obj,x,z,C,r,d,l,Phase,Temp) % I need the data of the object for this particular function
        [Phase] = generate_accretion_prism(obj,A,Phase,Boundary)
        function [Phase] = fill_weak_zone(obj,Phase,Z)
            ind = ~isnan(obj.d_slab(:))& Phase(:)~= obj.phase_prism{1} &  Phase(:)~= obj.phase_prism{2};%& Z(:) <obj.crust_depth;
            Phase(ind) = obj.ph_WZ;
        end
        
    end
end


