classdef BoundsT
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
        type  = 'Composite'; % elipses %circle %linear
        W      = 1.0      ; % W => can be the radius 
        L       = 1.0     ; % L=> can be the radius 
        x1     % x1 
        x2     % x2 
        y1     %
        y2
        c    = [0,0];           % center (x,y) coordinate of the center
       
        % A{1} => Vector P1-P2 of a boundary
        % A{2} => Type 'circular', 'none', 'none'
        % A{3} => Vector containing v(1) x_c, v(2), y_c, v(3), radius of
        % curvature
        % 
        A    = {[],['none'],[]};      % 4 value {x1,x2,R,c} | {}
        B    = {[],['none'],[]};      % 4 value {x1,x2,R,c} | {}
        C    = {[],['none'],[]};      %
        D    = {[],['none'],[]};       %
        angle = NaN          % primary axis direction w.r.t. current axis 
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
        function [ind] = find_ind_terrane(obj,A)
            if strcmp(obj.type,'Composite')
                % Composite Terrane are defined only by the 4 boundaries {i.e. are rectangle with modified margin}

                  [ind] = find_composite_terrane(obj,A); 
                
            end



        end

         % Arc Circumference like terrane
         [Size,CC] = modify_terrane_limits(obj,R,arc_length,Boundary,c)
         [obj] = Create_arc_circumference_margin(obj,R,Boundary,arcLength) 
         [y] = circumference_margin(obj,Boundary,x) % Find the coordinate of the boundary
         [s] = arc_length_CB(obj,x,y,Boundary)      % For a given y and x -> find the arc_length
         [B] = transform_coordinate(obj,A,Boundary,theta) % Transform the local coordinate system as such that within the area of the boundary the coordinate are transformed. 
         [s] = compute_arc_length_circle(obj,x,Boundary); % Compute per each x,y => the length 
    end
end

