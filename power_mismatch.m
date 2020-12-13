function mismatch = power_mismatch(Ps, Qs, G, B, V, delta, ...
    n_bus, pq_bus_id)

    P = zeros(n_bus,1);
    Q = zeros(n_bus,1);

    % calculating the active and reactive power at each bus
    %{
        P(i) = sum(j=1->n) |Vi||Vj|(Gij * cos(delta_i - delta_j) + 
                                    Bij * sin(delta_i - delta_j)
        Q(i) = sum(j=1->n) |Vi||Vj|(Gij * sin(delta_i - delta_j) - 
                                    Bij * cos(delta_i - delta_j)
    %}
    for i = 1 : n_bus
        for j = 1 : n_bus
            P(i) = P(i) + V(i)*V(j)*(G(i,j)*cos(delta(i)-delta(j)) + ...
                B(i,j)*sin(delta(i)-delta(j)));
            Q(i) = Q(i) + V(i)*V(j)*(G(i,j)*sin(delta(i)-delta(j)) - ...
                B(i,j)*cos(delta(i)-delta(j)));
        end
    end
    
    
    delta_P = Ps - P;
    delta_Q = Qs - Q;

    % since Q is unknown for PV bus, we only calculate for Q for PQ bus
    delta_Q = delta_Q(pq_bus_id);
    mismatch = [delta_P(2:end);delta_Q];    
end