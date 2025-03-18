function [ind] = find_composite_terrane(obj,A)
% => create boundary array:
%=
if isnan(obj.angle)
    fx1 =  A.Xpart(:)>=boundary_criteria(obj,A.Ypart(:),'A',obj.x1);
    fx2 =  A.Xpart(:)<=boundary_criteria(obj,A.Ypart(:),'C',obj.x2);
    fy2 =  A.Ypart(:)<=boundary_criteria(obj,A.Xpart(:),'B',obj.y2);
    fy1 =  A.Ypart(:)>=boundary_criteria(obj,A.Xpart(:),'D',obj.y1);
    ind = fx1==1 & fx2 == 1 & fy1 ==1 & fy2 ==1; 

end

end

function [dep_var] = boundary_criteria(obj,ind_var,Boundary,const)
% input variable
%==========================================================================
% ind_var => indipendent variable
% Boundary -> what is the type of Boundary
% Output:
% dep_var => dependent variable
%==========================================================================

if strcmp(obj.(Boundary){2},'none')
    dep_var = const;
elseif strcmp(obj.(Boundary){2},'Circular')
    dep_var = obj.circumference_margin(Boundary,ind_var);
else
    error('Unknown boundary type, spell it correctly')
end

end