classdef Thermal_Information 
    properties
        TP      % Mantle potential temperature in C
        TS      % Surface temperature in C
        k = 3  ;    % Conductivity
        Cp = 1050  ;   % Heat Capacity
        rho  =3000;   % average density
        Hr   = 0   % radiogenic heat production 
        kappa  
        Ph_Ast     
        Ph_Air
    end
    
    methods
        function obj=kappa_calculation(obj)
            %THERMAL_INFORMATION Construct an instance of this class
            %   Detailed explanation goes here
            obj.kappa = obj.k/(obj.rho*obj.Cp);
        end
        
    end
end

