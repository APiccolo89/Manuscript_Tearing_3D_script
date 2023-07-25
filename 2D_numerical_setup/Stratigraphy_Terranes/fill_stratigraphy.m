
function [Phase] = fill_stratigraphy(Z,Phase,Terranes,indx)
% function that allows to fill up the phases 

% Check if the type of terranes (i.e., Trench or normal Terranes?)
if ismember(Terranes,'Stratigraphy_Oceanic')
    ph = Terranes.Stratigraphy_Oceanic.phases;  % phases of the terranes
    t_Tk = Terranes.Stratigraphy_Oceanic.Tk;    % Thickness stratigraphy
else
    ph = Terranes.Stratigraphy.phases;         
    t_Tk = Terranes.Stratigraphy.Tk; 
end
% Extract relevant information 
T_Tk = [0.0, t_Tk(end)]; 
% Extract relevant information
T = T_Tk(1);

% Loop for filling the layered structured
for i=1:length(Terranes.Phases)
    if (i == length(Terranes.Phases))
        B=T_Tk(2);
    else
        B = t_Tk(i+1);
    end
    if ~isempty(indx)
        ind = Z < T & Z>=B & indx>0 & Phase == 0;
    else
        ind = Z < T & Z>=B;
    end
    Phase(ind) = ph(i);
    T=B;
    ind = [];
end
if ismember(Terranes,'cont')
    TT = Terranes.Stratigraphy.Continental_stratigraphy;
    T = TT.Tk(1);
    for i=1:length(TT.phases)
        if (i == length(TT.phases))
            B=TT.Tk(end);
        else
            B = TT.Tk(i+1);
        end

        ind = Z < T & Z>=B & cont>0;

        Phase(ind) = TT.phases(i);
        T=B;
        ind = [];
    end
end

end

