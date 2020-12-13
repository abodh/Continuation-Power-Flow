function [Q_lim_status, bus_data] = Qlim(Volt, Angle, bus_data, G, B,...
        base_MW, pv_bus_id, n_bus)
    
    % Q-LIMITS
    % if no-limits, set it to (-inf, inf) as the limit
    bus_data.data(find(bus_data.data(:,12) == 0),12) = inf;
    bus_data.data(find(bus_data.data(:,13) == 0),13) = -inf;

    % Q_max and Q_min (extracting Q limits)
    Q_max = bus_data.data(:,12);
    Q_min = bus_data.data(:,13);

    % check the Q-limit
    Q = zeros(n_bus,1);
    for i = 1 : n_bus
        for j = 1 : n_bus
            Q(i) = Q(i) + Volt(i)*Volt(j)*(G(i,j)*sin(Angle(i)-Angle(j))...
                - B(i,j)*cos(Angle(i)-Angle(j)));
        end
    end
    
    % Qg = Qsch + Qload
    Q_g = Q + (bus_data.data(:,7)/base_MW);
    Q_g = Q_g * base_MW;
    
    % check if the Q-limit violates
    check_id_max = pv_bus_id(find(Q_g(pv_bus_id) > Q_max(pv_bus_id)));
    check_id_min = pv_bus_id(find(Q_g(pv_bus_id) < Q_min(pv_bus_id))); 
    if isempty(check_id_max) && isempty(check_id_min)
        Q_lim_status = 0;
        fprintf("Power flow ran successfully without hitting Q-limits \n")        
    else
        bus_data.data(check_id_max,9) = bus_data.data(check_id_max,12);
        bus_data.data(check_id_min,9) = bus_data.data(check_id_min,13);
        bus_data.data(check_id_min,3) = 0;
        bus_data.data(check_id_max,3) = 0;
        Q_lim_status = 1;
        fprintf("Q-limit hit at bus -> %d \n",[check_id_max check_id_min])
    end
end