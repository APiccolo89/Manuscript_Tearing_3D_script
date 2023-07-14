function  [Phase] = generate_accretion_prism(A,Terranes,Phase)

C = Terranes.Trench_properties.C;
d_p = [C(1)+Terranes.position_end_prism 0.0];
s   = (d_p(2)-C(2))./(d_p(1)-C(1));
Phases = Terranes.Phases;
for i = 1:length(Phases)
    ind2 = A.Zpart>s.*(A.Xpart-C(1))+C(2)  & Phase == Phases(i);
    Phase(ind2) = Terranes.prism_phase(1);
end

end
