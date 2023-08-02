function [ind] = find_composite_terrane(obj,A)
% => create boundary array:
%=
if isnan(obj.angle)

    ind = A.Xpart(:)>=boundary_criteria(A.Ypart(:),'A') & A.Xpart(:)<=boundary_criteria(A.Ypart(:),'C') & A.Ypart(:)>=boundary_criteria(A.Xpart(:),'D') & ...
        A.Ypart(:)<=boundary_criteria(A.Xpart(:),'B');
end

end

function [dep_var] = boundary_criteria(ind_var,Boundary,const)
% input variable
%==========================================================================
% ind_var => indipendent variable
% Boundary -> what is the type of Boundary
% Output:
% dep_var => dependent variable
%==========================================================================

if strcmp(obj.(Boundary)(2),'none')
    dep_var = const;
elseif strcmp(obj.(Boundary)(2),'Circular')
    dep_var = obj.circumference_margin(Boundary,ind_var);
else
    error('Unknown boundary type, spell it correctly')
end

end