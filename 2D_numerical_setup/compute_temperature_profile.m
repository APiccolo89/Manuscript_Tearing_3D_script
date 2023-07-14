function [Temp]= compute_temperature_profile(A,Temp,Type,Gen,Terranes,indx)

T_tk = [0.0, Terranes.Stratigraphy(end)];
k = Terranes.K./(Terranes.Cp.*3300);
T_Age = Terranes.Age.*Terranes.secMyrsyear;
T_P   = Gen.T_P;
T_S   = Gen.T_S;
Tsl   = Gen.AvTS;


if Type == 1
    Zpart = A.Zpart; 
    ind_z = find(A.Zpart<T_tk(1) & A.Zpart>=T_tk(2) & indx == 1);
    ind = A.Zpart< T_tk(2);

else
    ind_z = A < T_tk(1) & A>=T_tk(2);
    ind = A < T_tk(2);
    Zpart = A; 
    

end
if isnan(T_Age)
    Temp(ind_z)=Tsl;
else

    erf_function = erf(Zpart(ind_z).*1000/2/(k*T_Age)^0.5);
    Temp(ind_z) = T_S - (T_P-T_S).*erf_function;
    Temp(ind) = T_P;

end




end
