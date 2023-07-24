function [Phase,Temp] = fill_layer(A,Terranes,Phase,Temp,Gen,cont)
% Select the layer and portion of the model:
ind = find_ind(obj); 
if ~isempty(indx(indx==1))
% Compute the thermal profile
    [Temp] = compute_temperature_profile(A,Temp,1,Gen,Terranes,indx);
% Fill the phase stratigraphy:
    [Phase] = fill_stratigraphy(A.Zpart,Phase,Terranes,indx,cont);
end
end
