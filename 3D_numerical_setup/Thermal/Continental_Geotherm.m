function [Temp] = Continental_Geotherm(obj,D,ind,Temp)
%CONTINENTAL_GEOTHERM Summary of this function goes here
%   Detailed explanation goes here
% Spell out the ingridient for the recepy: 
% ind : where to fill the temperature field. 

if isa(obj,'Passive_Margin')
    Moho = obj.Thermal_type_C.Moho;
    depth = obj.Thermal_type_C.Moho_d;
else
    Moho = obj.Thermal_type.Moho;
    depth = obj.Thermal_type.Moho_d;
end
T_tk = [0.0, obj.Stratigraphy.Tk(end)]; % Stratigraphy of the object
TP    = obj.Thermal_information.TP;
TS    = obj.Thermal_information.TS; 
gr1 = (Moho-TS)./(depth-0.0); 
gr2 = (TP-Moho)./(T_tk(2)-depth);
crust_ = @(x) linear_gradient(TS,gr1,x);
lith_mant = @(x) linear_gradient(Moho,gr2,x);
ind_crust = D(:)>=depth & D(:)< 0.0 & ind==1;
ind_lithosphere = D(:)<depth & D(:)>= T_tk(2) & ind==1;
ind_ast         = D(:)<T_tk(2); 
if isa(obj,'Passive_Margin')
        gr2 = (TP-Moho)./(obj.Tk_X-depth);
        ind_lithosphere = D(:)<depth & D(:)>=obj.Tk_X;

        gr2 = gr2(ind_lithosphere==1);
        lith_mant = @(x) linear_gradient(Moho,gr2,x);
end

T_C = (crust_(D(ind_crust==1)-0.0));
T_L = (lith_mant(D(ind_lithosphere==1)-depth));
Temp(ind_crust==1) = T_C; 
Temp(ind_lithosphere==1)=T_L;
Temp(ind_ast==1) = TP; 
end

function [T] = linear_gradient(A,B,C)
 T = A+B.*C;
end

