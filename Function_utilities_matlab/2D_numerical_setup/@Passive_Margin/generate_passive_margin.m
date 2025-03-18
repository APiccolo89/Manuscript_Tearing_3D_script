
function [Phase,Temp] = generate_passive_margin(obj,C,A,Phase,Temp)


if ~isempty(obj.Direction)
    depo_X = obj.center_pos;
    depo_Z = obj.depo_center;
    
    if strcmp(obj.Direction,'left')
        % To DO, not now, but in the future

    else
        x_lim      = [C.Boundary(1),C.Boundary(3)];
        lim_depo = [x_lim(2)-obj.length, x_lim(2)];
        Depo_x   = lim_depo(1)+depo_X.*obj.length;
        Depo_z   = - depo_Z;
        
        if strcmp(obj.shape,'triangular')
            x        = [lim_depo(1),Depo_x,lim_depo(2),lim_depo(2)+0.1*obj.length];
            z        = [0.0, Depo_z,Depo_z,0.0];

        elseif strcmp (obj.shape,'trapezoidal')
            x        = [lim_depo(1),Depo_x,x_lim(2),x_lim(2)];
            z        = [0.0,Depo_z,Depo_z,0.0];
        else
            error('Dear user, it seems that you did a mistake: passive margin are {trapezoidal} or {triangular}, any permutation of wrong letter is not admissible.')
        end
        Lithos   = C.Stratigraphy.Tk(end);
        x_l      = [lim_depo(1),lim_depo(2),lim_depo(2)];
        z_l      = [Lithos,Lithos,Lithos+obj.d_lithos];

    end
    [in,~] = inpolygon(A.Xpart,A.Zpart,x,z);
    Phase(in) = obj.ph_pas_m;
    [in2,~]  = inpolygon(A.Xpart,A.Zpart,x_l,z_l);
    Phase(in2) = C.Thermal_information.Ph_Ast;
    Temp(in2) = C.Thermal_information.TP;
    % => Thermal information 
    % Compute the new thermal structure within the terrane area: 
    % select the point chosen point: 
    if strcmp(obj.Direction,'left')
    else
        ind = A.Xpart(:)>lim_depo(1) & A.Xpart(:)<lim_depo(2) & (Phase(:) ~= C.Thermal_information.Ph_Ast | isnan(Phase(:))); % chosen particles
        x_chosen = (A.Xpart(ind==1)-lim_depo(1))./obj.length; % weight of the average
        z_chosen = A.Zpart(ind==1);
        obj.Tk_X    = obj.Stratigraphy.Tk(end)+obj.d_lithos*x_chosen;
        T_prov1   = x_chosen.*nan;  
        T_prov2   = x_chosen.*nan; 
        [T_prov1] = HalfSpaceCooling(obj,z_chosen,ind(ind==1),T_prov1);
        [T_prov2] = Continental_Geotherm(obj,z_chosen,ind(ind==1),T_prov1);
        T         = T_prov1.*x_chosen + T_prov2.*(1-x_chosen); 
        Temp(ind==1) = T;

    end


end
end

