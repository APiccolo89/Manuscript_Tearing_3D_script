
function [A] = generate_passive_margin(obj,C,A)
%=========================================================================%
% Input:
% obj => the passive_margin object associated with a specific terranes
% C   => The terrane object 
% A   => The grid properties
% Phase and Temp = > the phase and temperature fields, that are updated
% within this function. 
% Approach: The actual version is only accounting rectangular terranes, with
% different kind of boundary (circular or parabolic, whatsoever). In order
% to avoid the generation of additional variable, and to communicate the
% grid properties, i create B (which is A in which i switch X-Y, making all
% the function related to boundarie being an explicit function of "x").
% Simply I rotate of 90 degree the reference system. 
% If the boundary is parabolic or circle like, i compute the deflection and
% transform locally the coordinate to do the phase filling without doing
% obnoxious loop. 
%==========================================================================

if ~isempty(obj.Boundary_terrane_list)
    list_boundary = numel(obj.Boundary_terrane_list);
    depo_X = obj.center_pos;
    depo_Z = obj.depo_center;
    for ib=1:list_boundary
        %=============================================== Left _ kind
        %========
        if strcmp(obj.Direction,'left') || strcmp(obj.Boundary_terrane_list(ib),'A') ||strcmp(obj.Boundary_terrane_list(ib),'D')
            % I'm too lazy to change most of the function, so, since we are
            % speaking of linear terrane, i will just focus on making it
            % work only for the specific setup. Matlab is unsuitable for
            % having complex geometries, due to the inherent slowness of
            % the function and the necessity to vectorize.
            if strcmp(obj.Boundary_terrane_list(ib),'D')
                lim_depo = [C.Boundary.y1+obj.length, C.Boundary.y1];
                Depo_x   = lim_depo(1)-depo_X.*obj.length;
                Depo_z   = - depo_Z;
                l1 = C.Boundary.x1;
                l2 = C.Boundary.x2;
            else

                lim_depo = [C.Boundary.x1+obj.length, C.Boundary.x1];
                Depo_x   = lim_depo(1)+depo_X.*obj.length;
                Depo_z   = - depo_Z;
                l1 = C.Boundary.y1;
                l2 = C.Boundary.y2;
            end
            [x,z,x_l,z_l] = find_coordinates_object(obj,C,lim_depo,Depo_x,Depo_z,1);

            %=============================================== Right _ kind
            %========

        elseif strcmp(obj.Direction,'right') || strcmp(obj.Boundary_terrane_list(ib),'B') ||strcmp(obj.Boundary_terrane_list(ib),'C')
            if strcmp(obj.Boundary_terrane_list(ib),'C')
                % local transformation of the variable name in the
                % structure. This is a way
                lim_depo = [C.Boundary.x2-obj.length, C.Boundary.x2];
                Depo_x   = lim_depo(1)+depo_X.*obj.length;
                Depo_z   = - depo_Z;
                l1 = C.Boundary.y1;
                l2 = C.Boundary.y2;
            else
                lim_depo = [C.Boundary.y2-obj.length, C.Boundary.y2];
                Depo_x   = lim_depo(1)+depo_X.*obj.length;
                Depo_z   = - depo_Z;
                l1 = C.Boundary.x1;
                l2 = C.Boundary.x2;
            end

            [x,z,x_l,z_l] = find_coordinates_object(obj,C,lim_depo,Depo_x,Depo_z,2);

        end
        % Transform  the coordinate of the system
        [B] = C.Boundary.transform_coordinate(A,obj.Boundary_terrane_list{ib});       
        iy =  B.Ypart>=l1 & B.Ypart<=l2;
        [in,~] = inpolygon(B.Xpart,B.Zpart,x,z);
        A.Phase(in==1 & iy==1) = obj.ph_pas_m;
        % => Thermal information
        % Compute the new thermal structure within the terrane area:
        % select the point chosen point:
        if strcmp(obj.Direction,'left') || strcmp(obj.Boundary_terrane_list(ib),'A') ||strcmp(obj.Boundary_terrane_list(ib),'D')
            ind = B.Xpart(:)>lim_depo(2) & B.Xpart(:)<lim_depo(1) & (A.Phase(:) ~= C.Thermal_information.Ph_Ast | isnan(A.Phase(:))) & iy(:)==1; % chosen particles
        else
            ind = B.Xpart(:)>lim_depo(1) & B.Xpart(:)<lim_depo(2) & (A.Phase(:) ~= C.Thermal_information.Ph_Ast | isnan(A.Phase(:))) & iy(:)==1; % chosen particles
        end
        x_chosen = abs((B.Xpart(ind==1)-lim_depo(1))./obj.length); % weight of the average
        z_chosen = B.Zpart(ind==1);
        obj.Tk_X    = obj.Stratigraphy.Tk(end)+obj.d_lithos*x_chosen;
        T_prov1   = x_chosen.*nan;
        T_prov2   = x_chosen.*nan;
        [T_prov1] = HalfSpaceCooling(obj,z_chosen,ind(ind==1),T_prov1);
        if strcmp(obj.Thermal_type_C.Type,'HalfSpaceCooling')
            [T_prov1] = HalfSpaceCooling(obj,z_chosen,ind(ind==1),T_prov1);

        else
            [T_prov2] = Continental_Geotherm(obj,z_chosen,ind(ind==1),T_prov2);
        end
        T         = T_prov1.*x_chosen + T_prov2.*(1-x_chosen);
        A.Temp(ind==1) = T;
        % Lithosphere correction
        [in2,~]  = inpolygon(B.Xpart,B.Zpart,x_l,z_l);
        A.Phase(in2(:)==1 &  iy(:)==1) = C.Thermal_information.Ph_Ast;
        A.Temp(in2(:)==1 & iy(:)==1) = C.Thermal_information.TP;
    end
end
end



function [x,z,x_l,z_l] = find_coordinates_object(obj,C,lim_depo,Depo_x,Depo_z,Direction)

if strcmp(obj.shape,'triangular')
    if Direction == 1
        x        = [lim_depo(2),Depo_x,lim_depo(1)];


    else
        x        = [lim_depo(1),Depo_x,lim_depo(2)];

    end
    z        = [0.0, Depo_z,0.0];

elseif strcmp (obj.shape,'trapezoidal')
    if Direction == 1
        x        = [lim_depo(2),Depo_x,lim_depo(1),lim_depo(1)];
    else
        x        = [lim_depo(1),Depo_x,lim_depo(2),lim_depo(2)];
    end
    z        = [0.0,Depo_z,Depo_z,0.0];
elseif strcmp(obj.shape, 'rectangular')
    if Direction == 1
        x        = [lim_depo(2),lim_depo(2),lim_depo(1),lim_depo(1)];
    else
        x        = [lim_depo(1),lim_depo(1),lim_depo(2),lim_depo(2)];
    end
    z        = [0.0,Depo_z,Depo_z,0.0];
else
    error('Dear user, it seems that you did a mistake: passive margin are {trapezoidal} or {triangular} or {rectangular}, any permutation of wrong letter is not admissible.')
end
Lithos   = C.Stratigraphy.Tk(end);
if Direction==1
    x_l      = [lim_depo(1),lim_depo(2),lim_depo(2)];
else
    x_l      = [lim_depo(1),lim_depo(2),lim_depo(2)];
end
z_l      = [Lithos,Lithos,Lithos+obj.d_lithos];
end

