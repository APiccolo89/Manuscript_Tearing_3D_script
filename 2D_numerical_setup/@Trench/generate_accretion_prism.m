function  [Phase] = generate_accretion_prism(obj,A,Phase)
x = A.Xpart(:);
z = A.Zpart(:);

C = [obj.Boundary(1),-obj.R];
d_p = [C(1)+obj.position_end_prism 0.0];
s   = (d_p(2)-C(2))./(d_p(1)-C(1));
ind2 = z(:)>s.*(x(:)-C(1))+C(2)  & (Phase(:) == ~isnan(Phase(:)) | Phase(:) ~= obj.Thermal_information.Ph_Ast) & isnan(obj.d_slab(:)) &z <0.0 & x>C(1);
Phase(ind2) = obj.phase_prism;

end
