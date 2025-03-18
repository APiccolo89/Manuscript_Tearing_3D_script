
function [Phase] = fill_stratigraphy(obj,Z,Phase,ind)
%==========================================================================
% function that allows to fill up the phases 
% Check if the type of terranes (i.e., Trench or normal Terranes?)
% This function loop over the layer of a certain formation, and just fill
% up the designed geographical area with the desidered phase. More in
% general, it also allow to fill w.r.t. to the distance of a given surface.
% ========================================================================
% Input: 
% obj   -> terrane/boundary object
% D     -> distance from the surface (i.e. Z, or d of the slab) 
% Phase -> Phase structure
% ind   -> index of the chosen one (particles that belong a certain
% geological object). 
% =========================================================================
% Output: 
% Phase -> Updated array of phase. 
%==========================================================================
if isa(obj,'Trench')
    ph = obj.Stratigraphy_Oceanic.phases;  % phases of the terranes
    t_Tk = obj.Stratigraphy_Oceanic.Tk;    % Thickness stratigraphy
    D    = obj.d_slab; 
    ind  = (~isnan(D));
else
    ph = obj.Stratigraphy.phases;    % phase vector: i.e., the phases that are required to fill up     
    t_Tk = obj.Stratigraphy.Tk;      % depth of interface 
    D    = Z; 
end
% Extract relevant information 
T_Tk = [0.0, t_Tk(end)]; % Thickness of the terrain
% Extract relevant information
T = T_Tk(1);

% Loop for filling the layered structured
for i=1:length(ph)
    if (i == length(ph))
        B=T_Tk(2);
    else
        B = t_Tk(i+1);
    end
    if ~isempty(ind)
        ind_ph = D(:) <=T & D(:)>=B & ind(:)==1;
    else
        ind_ph = D(:) <= T & D(:)>=B & ind(:) == 1 ;
    end
    Phase(ind_ph) = ph(i);
    T=B;
    ind_ph = [];
end
if isa(obj,'Trench')
    TT = obj.Stratigraphy_Continental;
    T = TT.Tk(1);
    for i=1:length(TT.phases)
        if (i == length(TT.phases))
            B=TT.Tk(end);
        else
            B = TT.Tk(i+1);
        end

        ind_c = D(:) <= T & D(:)>=B & obj.continent(:)>0;

        Phase(ind_c) = TT.phases(i);
        T=B;
        ind_c = [];
    end
end

end

