function Y_bus = Ybus(n_bus, n_branch, branch_imp, bus_imp, from, to)
    % Y-bus calculation
    % initialize the Y-bus according to the size of nodes
    Y_bus = zeros(n_bus);

    for i = 1 : n_branch
        % Yij = -(1/(rij + j.xij))/tap-setting 
        Y_bus(from(i),to(i)) = - (1/(branch_imp(i,1) + ...
            1i*(branch_imp(i,2))))/branch_imp(i,4);

        %Yij = Yji
        Y_bus(to(i),from(i)) = Y_bus(from(i),to(i));

        % considering tap-setting the admittance matrix looks like
        % Y = [Y/t^2  -Y/t
        %       -Y/t     Y]
        Y_bus(from(i),from(i)) = Y_bus(from(i),from(i)) + ...
            ((1/(branch_imp(i,1) + 1i*(branch_imp(i,2))))...
            /(branch_imp(i,4))^2) + 1i*0.5*branch_imp(i,3);
                                 
        Y_bus(to(i),to(i)) = Y_bus(to(i),to(i)) + (1/(branch_imp(i,1) ...
            + 1i*(branch_imp(i,2)))) + 1i*0.5*branch_imp(i,3);
    end

    for i  = 1 : n_bus
        % the individual buses will have their own shunt device 
        % It should also be included in the Y_bus; Yii = Yii + (Gi + j.Bi)
        Y_bus(i,i) = Y_bus(i,i) + bus_imp(i,1) + 1i*bus_imp(i,2);
    end
end