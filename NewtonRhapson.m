function [Volt, Angle, error_avg] = ...
    NewtonRhapson(tolerance, n_bus, n_pv, n_pq, pq_bus_id, V, delta, ...
    G, B, Ps, Qs)
    
    % Power flow iteration
    iter = 0;
    mismatch = ones(2 * n_bus - n_pv - 2, 1);
    
    while any(abs(mismatch(:)) > tolerance)
        iter = iter + 1;
        
        % accumulate V and delta for all iteration
        Volt(:,iter) = V;
        Angle(:,iter) = delta;
        
        % calculate Jacobian
        Jacob = Jacobian(V, delta, n_bus, n_pq, pq_bus_id, G, B);
        
        % calculate power mismatch
        mismatch = power_mismatch(Ps, Qs, G, B, V, delta, n_bus,...
            pq_bus_id);
        
        % find the error [del_delta | del_V]
        error = croutLU(Jacob, mismatch);
        
        % update V and delta
        delta(2:end) = delta(2:end) + error(1 : n_bus-1);
        V(pq_bus_id) = V(pq_bus_id) + error(n_bus : end);
        
        % average error
        error_avg(iter) = abs(mean(error(:)));
    end
end