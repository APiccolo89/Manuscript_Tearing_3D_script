function [obj,Phase,Temp] = find_slab_(obj,A,Weak_Slab,Phase,Temp,Boundary,theta)

if strcmp(obj.Type_Subduction,'Mode_1')

    [obj,Phase,Temp] = obj.find_slab_mode_1(A,Weak_Slab,Phase,Temp,Boundary,theta);

elseif strcmp(obj.Type_Subduction,'Mode_2')
    % Find top, mid and bottom surface of the slab

    [obj] = obj.compute_top_bottom_surface(obj); 
    % Find the slab 
    
    [obj,Phase,Temp] = obj.find_slab_mode_2(A,Weak_Slab,Phase,Temp,Boundary,theta);

else
    error('Oh pityful soul, why did you not read the not existant user guide. The two slab types are {Mode_1}, and {Mode_2}, repent your self, and pray for the salvation')
end
end



