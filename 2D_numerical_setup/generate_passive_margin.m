
function [Phase,Temp] = generate_passive_margin(A,Phase,Temp,Terranes,direction,Gen)


if ~isempty(direction)
    depo_X = Terranes.Passive_Margin_depo_center_pos;
    depo_Z = Terranes.Passive_Margin_depo_center;
    Length = Terranes.Passive_Margin_length;
    if strcmp(direction,'left')

    else

        lim_depo = [Terranes.x_lim(2)-Length, Terranes.x_lim(2)];
        Depo_x   = lim_depo(1)+depo_X.*Length;
        Depo_z   = - depo_Z;
        x        = [lim_depo(1),Depo_x,lim_depo(2),lim_depo(2)+0.1*Length];
        z        = [0.0, Depo_z,Depo_z,0.0];
        Lithos   = Terranes.Stratigraphy(end);
        x_l      = [lim_depo(1),lim_depo(2),lim_depo(2)];
        z_l      = [Lithos,Lithos,Lithos+Terranes.Passive_Margin_d_lithos];

    end
    [in,~] = inpolygon(A.Xpart,A.Zpart,x,z);
    Phase(in) = Terranes.Passive_Margin_phase(1);
    [in2,~]  = inpolygon(A.Xpart,A.Zpart,x_l,z_l);
    Phase(in2) = Gen.Ph_UM;
    Temp(in2) = Gen.T_P;
    k = Terranes.K./(Terranes.Cp.*3300);
    dT_Age = (Terranes.Passive_Margin_Age(2)-Terranes.Passive_Margin_Age(1))./(x_l(2)-x(1));
    ind = A.Xpart>=x_l(1) & A.Xpart(2)<x_l(2);
    Age = A.Xpart(ind).*0.0;
    Age = (A.Xpart(ind == 1 & in==1)-x_l(1))*dT_Age+Terranes.Passive_Margin_Age(1); 
    Age = Age.*(365.25.*60.*60.*24.*1e6); 
    T_P   = Gen.T_P;
    T_S   = Gen.T_S;
    erf_function = erf(A.Zpart(ind == 1 & in ==1).*1000./2./(k.*Age).^0.5);
    Temp(ind == 1 & in == 1) = T_S - (T_P-T_S).*erf_function;
    Temp(in2) = T_P;

end
end

