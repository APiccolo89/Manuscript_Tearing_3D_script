classdef Terrane
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    properties
        Type
        order 
        x_lim 
        y_lim                          = 'none' ;
        Phases 
        Stratigraphy
        Age 
        Passive_Margin                 = 'none';
        Passive_Margin_phase           
        Passive_Margin_depo_center     = 12;
        Passive_Margin_depo_center_pos = 0.2; 
        Passive_Margin_length          = 100; 
        Passive_Margin_d_lithos        = 20; 
        Trench                         ='none';
        Trench_properties              = [];
        Accretion_prism                = 'none';
        position_end_prism             = 200; 
        prism_phase                     
        secMyrsyear                    = 365.25*60*60*24*1e6;
        Cp                             = 1050; 
        K                              = 3.0 ;
        kappa                     
    end
    
end

