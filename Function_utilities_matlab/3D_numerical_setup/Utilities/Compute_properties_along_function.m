function [L] = Compute_properties_along_function(xa,xb,val,lr)
    
if strcmp(val{2},'linear')
    [L] = linear_function(xa,xb,val{1}(1),val{1}(2),lr); 
else
    error('Yet again to lazy and unemployed for doing it')
end

end

