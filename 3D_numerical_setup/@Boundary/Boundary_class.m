classdef Boundary_class
    %BOUNDARY_CLASS Summary of this class goes here
    %   Detailed explanation goes here
    %======================================================================
    %             D                  % A Boundary A 
    %      _______________           y1, x1,x2 +/- function boundary
    %      |             |           y1 c(xc,yc) => y1 = yc-W/2
    %      |             |           x1  c(xc,yc) => x1 = xc-L/2
    %  A   |      c      |   C       x2  c(xc,yc)  => x2 = xc+L/2
    %      |             |           y2  c(xc,yc)  => y2 = yc+W/2
    %      |_____________|
    %      
    %             B       
    %
    %
    %  x|
    %   |_ _  
    %       y 
    %======================================================================

    properties
        name 
        type ='Composite' ;%function "i.e., circular, elipses 
        A              % 4 value {x1,x2,function_handle,Lenght_function_handle} | {}
        B              % 4 value {x1,x2,function_handle,Lenght_function_handle} | {}
        C              %
        D              %
        angle          % primary axis direction w.r.t. current axis 
        %================= if the terrane is circular or elipses 
        R              % Radius 
        r              % radius 
        C              % center
                  % major axis direction 
     
    end
    
    methods
        function obj = Boundary_class(inputArg1,inputArg2)
            %BOUNDARY_CLASS Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

