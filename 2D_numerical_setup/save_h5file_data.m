function save_h5file_data(Terranes,Phase,Temp,Z)

%==========================================================================
% Save the relevant information about a setup on h5 file that needs to be
% read by a Python script later on. This function serves to prepare a
% proper workflow for data post processing for the Current project: i.e.
% effects of slab-break off on surface uplift, insights from 3D numerical
% simulation. Do not expect high degree of generalisation, dear user, but
% feel free to steal the idea for having an easy life. 
%==========================================================================
structure_fields = fieldnames(Terranes);
number_terranes = numel(structure_fields);
TB = [];
for it = 1:number_terranes
    t = Terranes.(structure_fields{it}); 
    sub_fields = fieldnames(t);
    S_data = [] %empty structure for introducing the variables 
    for is =1:length(sub_fields)
       if isa(t.(sub_fields{is}),'Passive_Margin') || isa(t.(sub_fields{is}),'Thermal_Type') || isa(t.(sub_fields{is}),'Thermal_Information') || isa(t.(sub_fields{is}),'BoundsT')
           sub_sub_Field = fieldnames(t.(sub_fields{is}));
           for is2 = 1:length(sub_sub_Field)
               S_data.(sub_fields{is}).(sub_sub_Field{is2}) = t.(sub_fields{is}).(sub_sub_Field{is2});
           end
       else
           S_data.(sub_fields{is}) = t.(sub_fields{is}); 
       end

    end
    if isa(t,"Trench")
        ip = t.Stratigraphy_Oceanic.phases; 
        for ipt = 1:length(ip)
            Phase(Phase==ip(ipt))=-1000; 
        end
        mean_TSlab = mean(Temp(Phase(:)==-1000 & Z(:)<=-100 & Z(:)>=-200));
        S_data.TS = mean_TSlab; 
    end
     TB.((structure_fields{it})) = struct(S_data); 
end

save('Test_Data_Base.mat','TB','-v7.3')
end
