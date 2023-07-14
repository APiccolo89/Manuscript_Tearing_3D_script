
function [Phase] = fill_stratigraphy(Z,Phase,Terranes,indx,cont)
T_Tk = [0.0, Terranes.Stratigraphy(end)];
T = T_Tk(1);


for i=1:length(Terranes.Phases)
    if (i == length(Terranes.Phases))
        B=T_Tk(2);
    else
        B = Terranes.Stratigraphy(i+1);
    end
    if ~isempty(indx)
        ind = Z < T & Z>=B & indx>0 & Phase == 0;
    else
        ind = Z < T & Z>=B;
    end
    Phase(ind) = Terranes.Phases(i);
    T=B;
    ind = [];
end
if ~isempty(cont)
    TT = Terranes.Trench_properties.CCS;
    T = TT.Stratigraphy(1);
    for i=1:length(TT.phases)
        if (i == length(TT.phases))
            B=TT.Stratigraphy(end);
        else
            B = TT.Stratigraphy(i+1);
        end

        ind = Z < T & Z>=B & cont>0;

        Phase(ind) = TT.phases(i);
        T=B;
        ind = [];
    end
end

end

