function [Phase,Temp] = fill_subduction(A,Terranes,Phase,Temp,Gen)
ind_x=[];
cont = [];

[Layout,~,cont] = find_slab_(A,Terranes.Trench_properties,'Slab');
% Set the temperature
if strcmp(Terranes.Trench_properties.Temperature,'McKenzie')
    [Temp] = compute_temperature_profile(Layout,Temp,-1,Gen,Terranes,ind_x);%(A,Temp,Type,T_tk,Gen,Terranes,indx
    [Temp] = compute_temperature_profile_McKenzie(Layout,Temp,-1,Gen,Terranes,ind_x,A);%(A,Temp,Type,T_tk,Gen,Terranes,indx)
else
    [Temp] = compute_temperature_profile(Layout,Temp,-1,Gen,Terranes,ind_x);%(A,Temp,Type,T_tk,Gen,Terranes,indx)
end
% Correct Temperaure and Phase
id1 = min(A.Xpart(~isnan(Layout)));
id2 = max(A.Xpart(~isnan(Layout)));
id3 = min(A.Zpart(~isnan(Layout)));
ind_x1 = find(squeeze(A.Xpart(1,:,1))>=id1,1);
ind_x2 = find(squeeze(A.Xpart(1,:,1))>=id2,1);
ind_z1 = find(squeeze(A.Zpart(1,1,:))>=id3,1);
for i= ind_x1:ind_x2
    ind_L = find((Layout(1,i,:)==min(Layout(1,i,:))),1);
    Phase(:,i,ind_z1:ind_L) = Gen.Ph_Air;
    Temp(:,i,ind_z1:ind_L) = Gen.T_P;
end
[Phase] = fill_stratigraphy(Layout,Phase,Terranes,ind_x,cont);
[Layout] = find_slab_(A,Terranes.Trench_properties,'Weak');
if strcmp( Terranes.Trench_properties,'Mode_1')
    ind =(Layout<=0.0 & A.Zpart >= Terranes.Trench_properties.D_WZ & Phase ~=Gen.PrismPh);
    Phase(ind) = Gen.WZ;
else
    ind =(Layout>0.0 & A.Zpart >= Terranes.Trench_properties.D_WZ & Phase ~=Gen.PrismPh & A.Zpart<0.0);
    Phase(ind) = Gen.WZ;
end

end