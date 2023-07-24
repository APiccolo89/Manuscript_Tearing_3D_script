function [Temp]= compute_temperature_profile(obj,D,Temp,Type,Gen,ind)


T_tk = [0.0, obj.Stratigraphy(end)];
k = obj.K./(obj.Cp.*3300);
T_Age = obj.Age.*Terranes.secMyrsyear;
T_P   = Gen.T_P;
T_S   = Gen.T_S;
Tsl   = Gen.AvTS;

if Type == 1
    ind_z = find(D<T_tk(1) & D>=T_tk(2) & ind == 1);
    ind_o = A.Zpart< T_tk(2);

else
    ind_z = D < T_tk(1) & D>=T_tk(2);
    ind_o = D < T_tk(2);
    
end

if isnan(T_Age)
    Temp(ind_z)=Tsl;
else

    erf_function = erf(D(ind_z).*1000/2/(k*T_Age)^0.5);
    Temp(ind_z) = T_S - (T_P-T_S).*erf_function;
    Temp(ind_o) = T_P;

end

end
