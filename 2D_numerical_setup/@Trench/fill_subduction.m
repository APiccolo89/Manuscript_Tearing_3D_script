function [Phase,Temp] = fill_subduction(obj,A,Phase,Temp)
[obj,Phase,Temp] = obj.find_slab_(A,'Slab',Phase,Temp); % Since the correction for the phase and temperature is inevitably connected to the mid plane, i use this function to correct this array
A_time = cputime;
[Temp] = compute_temperature_profile(obj,A,[],Temp);
B_time = cputime; 
disp(['   Temperature field of the slab took ', num2str(B_time-A_time,3), ' seconds'])
A_time = cputime; 
[Phase] = fill_stratigraphy(obj,A,Phase,[]);
B_time = cputime; 
disp(['   Phase field of the slab took ', num2str(B_time-A_time,3), ' seconds'])
% Fill up the weak zone of the slab
% [Layout] = find_slab_(A,'Weak');
% if strcmp( Terranes.Trench_properties,'Mode_1')
%     ind =(Layout<=0.0 & A.Zpart >= Terranes.Trench_properties.D_WZ & Phase ~=Gen.PrismPh);
%     Phase(ind) = Gen.WZ;
% else
%     ind =(Layout>0.0 & A.Zpart >= Terranes.Trench_properties.D_WZ & Phase ~=Gen.PrismPh & A.Zpart<0.0);
%     Phase(ind) = Gen.WZ;
% end
bla = 0; 

end