function [Phase,Temp] = fill_subduction(obj,A,Phase,Temp,Gen)
obj = obj.find_slab_(A,'Slab');
[Temp] = compute_temperature_profile(obj,A,[],Temp);

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
% Fill up the weak zone of the slab
[Layout] = find_slab_(A,'Weak');
if strcmp( Terranes.Trench_properties,'Mode_1')
    ind =(Layout<=0.0 & A.Zpart >= Terranes.Trench_properties.D_WZ & Phase ~=Gen.PrismPh);
    Phase(ind) = Gen.WZ;
else
    ind =(Layout>0.0 & A.Zpart >= Terranes.Trench_properties.D_WZ & Phase ~=Gen.PrismPh & A.Zpart<0.0);
    Phase(ind) = Gen.WZ;
end

end