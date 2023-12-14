function [A] = fill_subduction(obj,A)
% Boundary list
boundary_list = obj.Boundaries_list;
for ib = 1:numel(boundary_list)
    % Transform coordinate:
    % Transform  the coordinate of the system
    [B] = obj.Boundary.transform_coordinate(A,obj.Boundaries_list{ib},obj.theta{1}(1)); % transform the coordinate of the system
    B.Length_along = B.Xpart*nan; % I need this temporary variable
    % unavoidable
    if strcmp(boundary_list{ib},'A')||strcmp(boundary_list{ib},'C')
        ya = obj.Boundary.(boundary_list{ib}){1}(2);
        yb = obj.Boundary.(boundary_list{ib}){1}(4);
    else
        ya = obj.Boundary.(boundary_list{ib}){1}(1);
        yb = obj.Boundary.(boundary_list{ib}){1}(3);
    end
     
    B.Length_along(B.Ypart>=ya & B.Ypart<=yb) = obj.Boundary.compute_arc_length_circle(B.Ypart(B.Ypart>=ya & B.Ypart<=yb),obj.Boundaries_list{ib});

    if strcmp(obj.theta{2},'none')
        theta = abs(obj.theta{1}(1));
        [obj,A.Phase,A.Temp] = obj.find_slab_(B,'Slab',A.Phase,A.Temp,obj.Boundaries_list{ib},theta); % Since the correction for the phase and temperature is inevitably connected to the mid plane, i use this function to correct this array
    else
        % Place holder for a loop over y direction 
        % ALG: 
        % B.y(i) => send function => fill => update 
        %=================================================================%
        % next => create a loop over the y particle and compute all the
        % necessary stuff each loop imply a different theta => is sliced up
        % and then integrated: this imply that I need to introduce an
        % additional temporary voice on the object, that collect d/l. 
        % The ending goal is having a data array that collects the distance
        % and the lenght along the slab.
        %==================================================================
        % NB:The direction of the slab is not perpendicular to the boundary
        %, but depends on the x/y position. This is the best easy peasy
        %solution that I found for the project at hand. 
        %==================================================================
        error('God the seventh day rested. I simply run out of my contract, do by yourself this functionality. With godly love AP')
    end

    A_time = cputime;
    B.Temp = A.Temp;
    [A.Temp] = compute_temperature_profile(obj,B,[],obj.Boundaries_list{ib});
    B.Temp = []; 
    B_time = cputime;
    disp(['   Temperature field of the slab took ', num2str(B_time-A_time,3), ' seconds'])
    A_time = cputime;
    [A.Phase] = fill_stratigraphy(obj,B,A.Phase,[]);
    [A.Phase] = obj.generate_accretion_prism(B,A.Phase,boundary_list{ib});
    [obj,A.Phase,A.Temp] = obj.find_slab_(B,'Weak',A.Phase,A.Temp,boundary_list{ib},theta);
    [A.Phase] = obj.fill_weak_zone(A.Phase,A.Zpart);
    B_time = cputime;
    disp(['   Phase field of the slab, prism weakzone  took ', num2str(B_time-A_time,3), ' seconds'])
end
end

