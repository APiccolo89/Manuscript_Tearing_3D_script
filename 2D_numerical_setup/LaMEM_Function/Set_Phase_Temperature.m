function [Phase,Temp] =  Set_Phase_Temperature(A,Phase,Temp,Terranes,Gen)
%=========================================================================%
% Function that fills the phases.
%=========================================================================%
% Input Parameter:
%=========================================================================%
% A = Data Structure with the grid/particles data.
% Phase, Temp Layer of phase or initial temperature field.
% Terranes: object containing information of the terranes.
% Gen     : generic data structure containing the information that are
% general (i.e. temperature of air, potential temperature, phase of prism,
% or passive margin or weak zone).
%=========================================================================%
% Output parameter:
%=========================================================================%
% Phase,Temp update.
%=========================================================================%
% Structure of the function:
% Fill layer between fixed bounds. Add subduction zone, add accretionary
% prism or passive margin and transitional terranes. Each of the generation
% of terranes is handled with specific functions.
%=========================================================================%
trench = Terranes.Trench;
passive_margin = Terranes.Passive_Margin;
accretion_prism = Terranes.Accretion_prism;
cont = [];
A_fill_layer = cputime;
[Phase,Temp] = fill_layer(A,Terranes,Phase,Temp,Gen,cont);
B_fill_layer = cputime;
t = B_fill_layer-A_fill_layer;
disp(['           1. Fill the layer and initialise Temperature field...', num2str(round(t)), 's']);
if strcmp(trench,'Subduction')
    A_fill_Sub = cputime;
    [Phase,Temp] = fill_subduction(A,Terranes,Phase,Temp,Gen);
    B_fill_Sub = cputime;
    t = B_fill_Sub-A_fill_Sub;
    disp(['        2. Fill the Subduction and initialise its Temperature field...', num2str(round(t)), 's']);

end
if strcmp(accretion_prism,'Prism')
    A_fill_Prism = cputime;
    [Phase] = generate_accretion_prism(A,Terranes,Phase);
    B_fill_Prism = cputime;
    t = B_fill_Prism-A_fill_Prism;
    disp(['         2. Fill the Prism ...', num2str(round(t)), 's']);
end
if strcmp(passive_margin,'none') == 0
    % generate left/right passive margin
    A_fill_Passive = cputime;

    for i=1:length(passive_margin)
        direction = Terranes.Passive_Margin{i};
        [Phase,Temp] = generate_passive_margin(A,Phase,Temp,Terranes,direction,Gen);
    end
    B_fill_Passive = cputime;
    t = B_fill_Passive-A_fill_Passive;
    disp(['       3. Fill the Passive Margins...', num2str(round(t)), 's']);

end

end

