classdef Terrane
    properties
        Boundary                    = [@f_handle];  
        Phases 
        Stratigraphy
%         Passive_Margin                 = 'none';
%         Passive_Margin_phase           
%         Passive_Margin_depo_center     = 12;
%         Passive_Margin_depo_center_pos = 0.2; 
%         Passive_Margin_length          = 200; 
%         Passive_Margin_d_lithos        = 20; 
%         Passive_Margin_Age             = []; 
        Trench                         ='none';
        Trench_properties              = [];
%         Accretion_prism                = 'none';
%         position_end_prism             = 200; 
%         prism_phase                     
%         secMyrsyear                    = 365.25*60*60*24*1e6;
        Cp                             = 1050; 
        K                              = 3.0 ;
        kappa                     
    end
    methods
        function [Phase,Temp] = fill_terranes(obj,Phase,Temp)
            [Phase,Temp] = fill_layer(obj,A,Phase,Temp,Gen,cont);
        end

    end

end

