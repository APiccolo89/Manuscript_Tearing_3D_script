
function [Temp]= compute_temperature_profile_McKenzie(obj,Temp,A,Boundary)

% Spell out the necessary variables
k = obj.Thermal_information.k;
T_P   = obj.Thermal_information.TP;
T_S   = obj.Thermal_information.TS;
if strcmp(obj.Thermal_type.vel{2},'none')
    vl    = obj.Thermal_type.vel{1}(1);
    Re    = (obj.Thermal_information.rho*obj.Thermal_information.Cp.*vl.*obj.D0.*1000)./2./k;

else
    data_boundary = obj.Boundary.(Boundary);
    %Select the boundary = limits of the boundary
    [loc1, ~] = find(tril(data_boundary{1} == data_boundary{1}', -1)); % find the main boundary.
    if rem(loc1,2)==0
        ya = data_boundary{1}(1);
        yb = data_boundary{1}(3);
    else
        ya = data_boundary{1}(2);
        yb = data_boundary{1}(4);
    end
    vl    = Compute_properties_along_function(ya,yb,obj.Thermal_type.vel,A.Length_along);
    Re    = (obj.Thermal_information.rho*obj.Thermal_information.Cp.*vl.*obj.D0.*1000)./2./k;
    Re    = Re(~isnan(obj.d_slab));
end



% Compute the Thermal Reynolds number
% Compute the lenght of the slab starting from the trench
% make dimensionless d and l;
l = obj.l_slab(:)./obj.D0;
d = obj.d_slab(:)./obj.D0;
%======== Compute the temperature profile
n = 26;
t_prov = Temp(:).*0.0;
Sigma    = t_prov;
Sigma(~isnan(d)) = 0.0;
% Dividi et impera
dn = 1-abs(d(~isnan(d)));
for i=1:n
    a = (-1).^(i)./(i.*pi);
    b = (Re-(Re.^2+i^2*pi^2).^(0.5)).*l(~isnan(d));
    c = sin(i.*pi.*(1-abs(d(~isnan(d)))));
    e = exp(b);
    Sigma(~isnan(d)) = Sigma(~isnan(d))+ a.*e.*c;
end
t_prov(~isnan(d)) = (T_P)+2.*(T_P-T_S).*Sigma(~isnan(d));
t_prov(~isnan(d) & t_prov<0) = T_S;
% Decoupling depth-length
iz = obj.Decoupling_depth(2);
weight(~isnan(d)) = 0+(0.8-0.1)./((l(iz)+0.04)-0.01).*(l(~isnan(d))-0.1);
weight(weight>1) =1.0;
weight(weight<0) =0.01;

weight = weight';
Temp(~isnan(d))=(weight(~isnan(d))).*t_prov(~isnan(d))+(1-weight(~isnan(d))).*Temp(~isnan(d));
Temp = reshape(Temp,size(A.Xpart));
end
