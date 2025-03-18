
function [Temp]= compute_temperature_profile_McKenzie(obj,A,Temp,Boundary)
shape = size(Temp);
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
%=======================================================%
% Update %
% introducing scaling to use the object that I have

%l = obj.l_slab(:)./obj.D0;
%d = obj.d_slab(:)./obj.D0;
% 
sc = 1./obj.D0;
%======== Compute the temperature profile
n = 26;
t_prov = Temp(:).*0.0;
Sigma    = t_prov;
Sigma(~isnan(obj.d_slab)) = 0.0;
% Dividi et impera
for i=1:n
    a = (-1).^(i)./(i.*pi);
    b = (Re-(Re.^2+i^2*pi^2).^(0.5)).*obj.l_slab(~isnan(obj.d_slab)).*sc;
    c = sin(i.*pi.*(1-abs(obj.d_slab(~isnan(obj.d_slab)).*sc)));
    e = exp(b);
    Sigma(~isnan(obj.d_slab)) = Sigma(~isnan(obj.d_slab))+ a.*e.*c;
end
t_prov(~isnan(obj.d_slab)) = (T_P)+2.*(T_P-T_S).*Sigma(~isnan(obj.d_slab));
t_prov(~isnan(obj.d_slab(:)) & t_prov<0) = T_S;
% Decoupling depth-length
iz = obj.Decoupling_depth(2);
weight(~isnan(obj.d_slab)) = 0+(0.8-0.1)./((obj.l_slab(iz).*sc+0.04)-0.01).*(obj.l_slab(~isnan(obj.d_slab)).*sc-0.1);
weight(weight>1) =1.0;
weight(weight<0) =0.01;
weight = weight';
Temp(~isnan(obj.d_slab))=(weight(~isnan(obj.d_slab))).*t_prov(~isnan(obj.d_slab))+(1-weight(~isnan(obj.d_slab))).*Temp(~isnan(obj.d_slab));
Temp = reshape(Temp,shape);
end
