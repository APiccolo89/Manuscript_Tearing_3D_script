function [A] = fill_layer(obj,A)
% Select the layer and portion of the model:
time_A = cputime; 
ind = obj.Boundary.find_composite_terrane(A);
if ~isempty(ind(ind==1))
% Compute the thermal profile
% Phase is easy to fill up. Temperature has many cases and I would like to
% avoid to 
    [A.Temp] = compute_temperature_profile(obj,A,ind,[]);
% Fill the phase stratigraphy:
    [A.Phase] = fill_stratigraphy(obj,A.Zpart,A.Phase,ind);
end
time_B = cputime; 
disp(['    Terrane temperature and stratigraphy took ', num2str((time_B-time_A),3),' seconds']);
end

% function [ind] = find_ind(obj,A)
% % Place holder for more complex function in the futre
% limx = [obj.Boundary(1),obj.Boundary(3)];
% limy = [obj.Boundary(2),obj.Boundary(4)]; % Sherlok, you confuse 2 with 3 losing 1 Hour of your time, this comment shall be your memento mori
% ind = A.Xpart(:)>=limx(1) & A.Xpart(:)<=limx(2) & A.Ypart(:)>=limy(1) & A.Ypart(:)<=limy(2);
% end
