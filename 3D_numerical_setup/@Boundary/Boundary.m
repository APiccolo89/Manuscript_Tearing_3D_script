classdef Boundary
    %BOUNDARY_CLASS Summary of this class goes here
    %   Detailed explanation goes here
    %======================================================================
    %             B                  % A Boundary A 
    %      _______________           y1, x1,x2 +/- function boundary
    %      |             |           y1 c(xc,yc) => y1 = yc-W/2
    %      |             |           x1  c(xc,yc) => x1 = xc-L/2
    %  A   |      c      |   C       x2  c(xc,yc)  => x2 = xc+L/2
    %      |             |           y2  c(xc,yc)  => y2 = yc+W/2
    %      |_____________|
    %      
    %             D           For now the boundary class will consider
    %             rectangular terranes, whose angle is // to the y axis. So
    %             perpendicular to the x axis. In the future, if I have
    %             time, I would like to introduce an automatic coordinate
    %             detection for any rotation. (i.e., apply a rotation
    %             matrix to the coordinate as a function of a theta angle
    %             defined by the user) 
    %
    %
    %  x|
    %   |_ _  
    %       y 
    %======================================================================

    properties
        name 
        type (1,:){mustBeTextScalar} = 'Composite'; % elipses %circle. 
        W     {mustBeGreaterThan(W,0)} = 1.0      ; % W => can be the radius 
        L      {mustBeGreaterThan(L,0)} = 1.0     ; % L=> can be the radius 
        x1     % x1 
        x2     % x2 
        y1     %
        y2
        c    = [0,0];           % center (x,y) coordinate of the center
        A    = {[],[]};      % 4 value {x1,x2,function_handle,Lenght_function_handle} | {}
        B    = {[],[]};         % 4 value {x1,x2,function_handle,Lenght_function_handle} | {}
        C    = {[],[]};      %
        D    = {[],[]};       %
        angle          % primary axis direction w.r.t. current axis 
     
    end
    
    methods
        function [obj] = compute_coordinate_boundary_composite(obj)
            % For now, this function, deal with simple rectangular terrane
            % that has not any direction 
            % compute x1,x2,y1,y2
            if length(obj.c) ~= 2 
                error('Center must be a point within x-y plane')
            end
            obj.x1 = obj.c(1)-obj.L/2; 
            obj.x2 = obj.c(1)+obj.L/2; 
            obj.y1 = obj.c(2)-obj.W/2;
            obj.y2 = obj.c(2)+obj.W/2; 

            obj.A{1} = [obj.x1,obj.y1,obj.x1,obj.y2];
            obj.B{1} = [obj.x1,obj.y2,obj.x2,obj.y2];
            obj.C{1} = [obj.x2,obj.y2,obj.x2,obj.y1];
            obj.D{1} = [obj.x2,obj.y1,obj.x1,obj.y1]; 
        end
        
         [obj] = modify_boundary_limits(obj,arcLength)

      
    end
end

