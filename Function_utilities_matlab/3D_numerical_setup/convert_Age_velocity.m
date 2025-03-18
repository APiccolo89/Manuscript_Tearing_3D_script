function [converted] = convert_Age_velocity(value,type)
        Myrs_sec = 365.25*60*60*24*1e6;
        cm_to_meter = 100; 
        if type == 1
            converted = value*Myrs_sec; 
        elseif type ==2 
            converted = value/(Myrs_sec/1e6)/cm_to_meter; 
        end
end
