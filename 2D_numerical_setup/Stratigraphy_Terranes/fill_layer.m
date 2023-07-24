function [Phase,Temp] = fill_layer(obj,A,Phase,Temp,Gen,cont)
% Select the layer and portion of the model:
ind = find_ind(obj,A);
if ~isempty(indx(indx==1))
% Compute the thermal profile
    [Temp] = compute_temperature_profile(A,Temp,1,Gen,Terranes,indx);
% Fill the phase stratigraphy:
    [Phase] = fill_stratigraphy(A.Zpart,Phase,Terranes,indx,cont);
end
end
function [ind] = find_ind(obj)
% Place holder for more complex function in the futre
lim = obj.Boundary; 
ind = A.Xpart(:)>lim(1) & A.Xpart(:)>lim(2) ;



end
