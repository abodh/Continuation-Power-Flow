function Jacob_matrix = Jacobian(V, delta, n_bus, n_pq, pq_bus_id, G, B)
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
            P(i) = P(i) + V(i)*V(j)*(G(i,j)*cos(delta(i)-delta(j))...
                + B(i,j)*sin(delta(i)-delta(j)));
            Q(i) = Q(i) + V(i)*V(j)*(G(i,j)*sin(delta(i)-delta(j))...
                - B(i,j)*cos(delta(i)-delta(j)));
        end
    end

    %{
    For IEEE 14 bus system
    slack bus : 1
    PV buses = 4
    PQ bus = 9
    size of Jacobian matrix = 2*n_bus - n_pv (4) - 2 = 22
    %}

    % size of each of the Jacobian sub matrices
    J11 = zeros(n_bus-1, n_bus-1);
    J12 = zeros(n_bus-1, n_pq); % remove the bus (cols) where V is given
    J21 = zeros(n_pq, n_bus-1); % remove the bus (rows) where Q is unknown
    J22 = zeros(n_pq, n_pq);
    
    %{
    refer to the link below for sub-matrices equations:
    https://bit.ly/2GU1Xx0
    
    B. Sereeter, C. Vuik, and C. Witteveen
    REPORT 17-07 
    On a comparison of Newton-Raphson solvers for power flow problems
    %}
    
    % Calculating J11
    for i = 2 : n_bus
        for k = 2 : n_bus
            if (i == k)
                J11(i-1,k-1) = - Q(i) - (V(i)^2 * B(i,i));
            else
                J11(i-1,k-1) = V(i)*V(k)*(G(i,k)*sin(delta(i)-delta(k))...
                    - B(i,k)*cos(delta(i)-delta(k)));
            end
        end
    end

    % Calculating J21
    for i = 2 : n_pq + 1 
        j = pq_bus_id(i-1);
        for k = 2 : n_bus
            if (j == k)
                J21(i-1,k-1) = P(j) - (V(j)^2 * G(j,j));
            else
                J21(i-1,k-1) = -V(j)*V(k)*(G(j,k)*cos(delta(j)-delta(k))...
                    + B(j,k)*sin(delta(j)-delta(k)));
            end
        end
    end


    % Calculating J12
    for i = 2 : n_bus  
        for k = 2 : n_pq + 1
            j = pq_bus_id(k-1);
            if (i == j)
                J12(i-1,k-1) = P(j) + (V(j)^2 * G(j,j));
            else
                J12(i-1,k-1) = V(i)*(G(i,j)*cos(delta(i)-delta(j))...
                    + B(i,j)*sin(delta(i)-delta(j)));
            end
        end
    end

    % Calculating J22
    for i = 2 : n_pq + 1
        j = pq_bus_id(i-1);
        for k = 2 : n_pq + 1
            l = pq_bus_id(k-1);
            if (j == l)
                J22(i-1,k-1) = Q(j) - (V(j)^2 * B(j,j));
            else
                J22(i-1,k-1) = V(j)*(G(j,l)*sin(delta(j)-delta(l))...
                    - B(j,l)*cos(delta(j)-delta(l)));
            end
        end
    end
    
    % combining the jacobian matrix
    Jacob_matrix = [J11 J12; J21 J22];    
end