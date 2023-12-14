function  [Phase] = generate_accretion_prism(obj,A,Phase,Boundary)
x = A.Xpart(:);
z = A.Zpart(:);
data_boundary = obj.Boundary.(Boundary);
%Select the boundary = limits of the boundary
[loc1, ~] = find(tril(data_boundary{1} == data_boundary{1}', -1)); % find the main boundary.
if rem(loc1,2)==0
    ya = data_boundary{1}(1);
    yb = data_boundary{1}(3);
    xc = data_boundary{1}(2);
else
    ya = data_boundary{1}(2);
    yb = data_boundary{1}(4);
    xc = data_boundary{1}(3);

end

C = [xc,-obj.R(2)];
if ~isempty(obj.C_prism)
    CP = C; 
else
    CP = obj.C_prism; 
end
d_p = [CP(1)+obj.position_end_prism 0.0];
s   = (d_p(2)-obj.prism_depth)./(d_p(1)-C(1)+obj.position_end_prism*0.8);
ind2 = z(:)>s.*(x(:)-C(1))+C(2)  & (Phase(:) == ~isnan(Phase(:)) | Phase(:) ~= obj.Thermal_information.Ph_Ast) & isnan(obj.d_slab(:)) &z <0.0 & x>C(1) & A.Ypart(:)>=ya & A.Ypart(:)<=yb ;
Phase(ind2) = obj.phase_prism{1};
ind3 = Phase(:) == obj.phase_prism{1} & z<=obj.Prism_lc_depth;
Phase(ind3) = obj.phase_prism{2};

end
