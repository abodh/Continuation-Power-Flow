function [Volt, Angle, iter] = ...
    NewtonRhapson(tolerance, n_bus, n_pv, n_pq, pq_bus_id, V, delta, ...
    G, B,Y_bus, Ps, Qs)
    
    % Power flow iteration
    iter = 0;
    % mismatch = ones(2 * n_bus - n_pv - 2, 1);
    mismatch = power_mismatch(Ps, Qs, G, B, V, delta, n_bus,...
            pq_bus_id);
        
    counter = 1;
    % accumulate V and delta for all iteration
    Volt(:,counter) = V;
    Angle(:,counter) = delta;
    
%     while any(abs(mismatch(:)) > tolerance)
    while (abs(max(mismatch)) >= tolerance)
        if (iter >= 10)
            break
        end        
        iter = iter + 1;
        
        % calculate Jacobian
        Jacob = Jacobian(V, delta, n_bus, n_pq, pq_bus_id, G, B,Y_bus);
        
        % find the error [del_delta | del_V]
        error = croutLU(Jacob, mismatch);
        
        % update V and delta
        delta(2:end) = delta(2:end) + error(1 : n_bus-1);
        V(pq_bus_id) = V(pq_bus_id) + error(n_bus : end);
        
        % calculate power mismatch
        mismatch = power_mismatch(Ps, Qs, G, B, V, delta, n_bus,...
            pq_bus_id);
        
        % accumulate V and delta for all iteration
        counter = counter + 1;
        Volt(:,counter) = V;
        Angle(:,counter) = delta;
        
%         % average error
%         error_avg(iter) = abs(mean(error(:)));
%         abs(max(mismatch(:)))
%         if ((abs(max(mismatch(:))) >= tolerance && iter <= 10))
%             counter = counter + 1;
%             Volt(:,counter) = V;
%             Angle(:,counter) = delta;     
%         end
    end
end