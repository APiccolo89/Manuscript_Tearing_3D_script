function [Phase,Temp] = fill_subduction(obj,A,Phase,Temp)
[obj,Phase,Temp] = obj.find_slab_(A,'Slab',Phase,Temp); % Since the correction for the phase and temperature is inevitably connected to the mid plane, i use this function to correct this array
A_time = cputime;
[Temp] = compute_temperature_profile(obj,A,[],Temp);
B_time = cputime; 
disp(['   Temperature field of the slab took ', num2str(B_time-A_time,3), ' seconds'])
A_time = cputime; 
[Phase] = fill_stratigraphy(obj,A,Phase,[]);
[Phase] = obj.generate_accretion_prism(A,Phase);
[obj,Phase,Temp] = obj.find_slab_(A,'Weak',Phase,Temp);
[Phase] = obj.fill_weak_zone(Phase);
B_time = cputime; 

disp(['   Phase field of the slab, prism weakzone  took ', num2str(B_time-A_time,3), ' seconds'])

end

