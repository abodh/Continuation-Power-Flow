function [bus_imp, branch_imp, bus_data, branch_data] = ...
    data_extraction (bus_path, branch_path, System)

    % the bus and branch data have been separated to two files
    % both of the files are inside IEEE_(14/30)_bus_data folder
    
    %import the data for each buses
    bus_data = importdata(bus_path);
    
    if (System=="IEEE-30")
        bus_data.data = bus_data.data(:,2:end);
    end

    % bus_imp = [conductance susceptance] for each of the bus
    bus_imp = bus_data.data(:,14:15); 

    %import the data for each branches
    branch_data = importdata(branch_path);
    
    % ensures that non of the tap setting is 0
    branch_data(find(branch_data(:,15)==0),15) = 1;

    % branch_imp = [R X B tap_setting] for each of the branches
    branch_imp = [branch_data(:,7:9), branch_data(:,15)];
end