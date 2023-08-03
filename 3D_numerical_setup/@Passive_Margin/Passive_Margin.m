classdef Passive_Margin % => Composite class with continents ~ Default value
    %               
    %PASSIVE_MARGIN========================================================
    %   Special kind of geological object that needs to modify phase and
    %   temperature field of a narrow area of an attached continent. One
    %   day I'll introduce other option to have a geological like situation
    %   (i.e. introducing the passive margin in overriding plate). For now
    %   it is simply a way to introduce a rectangular and triangular shape
    %   object, modify the continental lithosphere such that is thinning
    %   towards the oceanic plate and has a smooth transition towards the
    %   thermal structure of the ocean plate/trench.
    %======================================================================
    %
    %
    %
    %
    %======================================================================
    
    properties
                Direction     = []       % in 2D 
                Stratigraphy         %Stratigraphy bounded continent   
                ph_pas_m            %sedimentary phase 
                depo_center     = 12; %Depth of the depocenter 
                center_pos = 0.2;     %Depocenter position
                shape         = 'triangular'; %'trapezoidal'
                length          = 200;% Horizontal distance from the border
                d_lithos        = 20; % thickness difference 
                Age             = 80; % Age of the next terrane
                Thermal_type_O  = Thermal_Type;
                Thermal_type_C = Thermal_Type;    % Thermal type of the attached continent     
                Thermal_information = Thermal_Information;
                Boundary_terrane_list = 'D'; 
                Tk_X           % Information local thickness
    end
    
    methods
        [Phase,Temp] = generate_passive_margin(obj,C,A,Phase,Temp)
    end
end

