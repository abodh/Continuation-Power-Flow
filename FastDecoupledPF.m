function [Volt_FD, Angle_FD, error_avg_FD] = ...
    FastDecoupledPF(tolerance, from, to, n_branch, n_bus, n_pv, n_pq, ...
    pq_bus_id, V, delta, G, B, Ps, Qs, branch_imp, bus_imp)

    mismatch = ones(2 * n_bus - n_pv - 2, 1);
    iter = 0;

    % uncomment the line below when same B matrix is used for the Jacobian
    % B11 = B(2:end,2:end);

    B_prime = zeros(n_bus,n_bus);

    % the formulation below is similar to the Y-bus formulation -> Ybus.m

    for i = 1 : n_branch
        % Yij = -((1/(rij + j.xij)) + j.Bij/2)/tap-setting 
        B_prime(from(i),to(i)) = - (1/(0 + 1i*(branch_imp(i,2))))...
            /branch_imp(i,4);

        %Yij = Yji
        B_prime(to(i),from(i)) = B_prime(from(i),to(i));

        % considering tap-setting the admittance matrix looks like
        % Y = [Y/t^2  -Y/t
        %       -Y/t     Y]
        B_prime(from(i),from(i)) = B_prime(from(i),from(i)) + ((1/(0 + ...
            1i*(branch_imp(i,2))))/(branch_imp(i,4))^2) + ...
            1i*0.5*branch_imp(i,3);

        B_prime(to(i),to(i)) = B_prime(to(i),to(i)) + (1/(0 + ...
            1i*(branch_imp(i,2)))) + 1i*0.5*branch_imp(i,3);
    end

    for i  = 1 : n_bus
        % the individual buses will have their own shunt device 
        % It should also be included in the Y_bus; Yii = Yii + (Gi + j.Bi)
        B_prime(i,i) = B_prime(i,i) + 0 + 1i*bus_imp(i,2);
    end

    B_prime_ = imag(B_prime);
    % remove the col/row corresponding to slack bus
    B11 = B_prime_(2:end,2:end);

    % forming the jacobian J22 where only PQ bus B values are considered
    for i = 1 : length(pq_bus_id)
        for j = 1 : length(pq_bus_id)
            B22(i,j) = B(pq_bus_id(i), pq_bus_id(j));
        end
    end

    % final Jacobian with J12 = J21 = 0
    Jacob_FD = -[B11 zeros(n_bus-1, length(pq_bus_id)); ...
        zeros(length(pq_bus_id),n_bus-1) B22];
    
    while any(abs(mismatch(:)) > 0.01)
        iter = iter + 1;
        Volt_FD(:,iter) = V;
        Angle_FD(:,iter) = delta;   

        % mismatch calculation
        mismatch = power_mismatch(Ps, Qs, G, B, V, delta, ...
            n_bus, pq_bus_id);

        % error [del_delta | del_V] using Crout's LU factorization]
        error = croutLU(Jacob_FD, mismatch);

        % updating delta and V
        delta(2:end) = delta(2:end) + error(1 : n_bus-1);
        V(pq_bus_id) = V(pq_bus_id) + error(n_bus : end);

        % average error
        error_avg_FD(iter) = abs(mean(error(:)));
    end
    
    function mismatch_FD = power_mismatch(Ps, Qs, G, B, V, delta, ...
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
                P(i) = P(i) + V(i)*V(j)*(G(i,j)*cos(delta(i)-delta(j)) ...
                    + B(i,j)*sin(delta(i)-delta(j)));
                Q(i) = Q(i) + V(i)*V(j)*(G(i,j)*sin(delta(i)-delta(j)) ...
                    - B(i,j)*cos(delta(i)-delta(j)));
            end
        end

        delta_P = Ps - P;
        delta_Q = Qs - Q;
        
        % since Q is unknown for PV bus, we only calculate for Q for PQ bus
        delta_Q = delta_Q(pq_bus_id);
        mismatch_FD = [delta_P(2:end);delta_Q];    
    end
end