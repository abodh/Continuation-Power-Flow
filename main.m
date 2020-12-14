%{
Continuation Power Flow
Author: Abodh Poudyal
Last updated: December 13, 2020
%}
clear;
clc;
% format short % to display less significant digits in the result 

%% Reading bus and branch data in common data format/ Initializations
% external function to extract the data from IEEE common data format

% select the system on which the power flow should be performed
% currently works for 'IEEE14' and 'IEEE30'
System = "IEEE-14";

bus_path = strcat(System,'bus_data/bus_data.txt');
branch_path = strcat(System,'bus_data/branch_data.txt');
[bus_imp, branch_imp, bus_data, branch_data] = ...
    data_extraction(bus_path, branch_path, System);

%{ 
to reduce the computational complexity, we will only compute 
for existing branches
%}

% define some important variables

% from which bus
from = branch_data(:,1);
% to which bus
to = branch_data(:,2);

% extract voltage data
V_flat = bus_data.data(:,11);  

% flat start means |V| = 1.0 pu and delta = 0
% exact values for PV and slack bus whereas flat start for the rest
V_flat(find(V_flat == 0)) = 1;
delta_flat = zeros(length(V_flat),1)*pi/180;;

% V = V_flat;
% delta = delta_flat;

% number of buses in the entire system
n_bus = length(bus_data.data(:,3));

% number of branches
n_branch = length(branch_imp);


% number of pq buses
n_pq = length(find(bus_data.data(:,3) == 0));

% number of PV buses 
n_pv = length(find(bus_data.data(:,3) == 2));

% stores an array of PQ bus IDs
pq_bus_id = find(bus_data.data(:,3) == 0);

% stores an array of PV bus IDs 
pv_bus_id = find(bus_data.data(:,3) == 2);

% iterate unless power mismatch < 0.01 (tolerance)
tolerance = 0.01;

% base power
base_MW = 100;

% scheduled power
Ps = (bus_data.data(:,8) - bus_data.data(:,6))/base_MW;
Qs = (bus_data.data(:,9) - bus_data.data(:,7))/base_MW;

% ek vector analysis
% check if a bus is pq
pq_bus_logic = (bus_data.data(:,3) ~= 2) & (bus_data.data(:,3) ~= 3);

% converts logical array to double
pq_bus_logic = double(pq_bus_logic);

% name the pq buses
pq_bus_logic(pq_bus_logic == 1) = transpose(1:length(pq_bus_logic ...
    (pq_bus_logic == 1)));

% positions and logic for ek vector
ek_positions = [transpose(1:n_bus) (pq_bus_logic == 0)  pq_bus_logic];

% K vector
K = [Ps(2:end);Qs(bus_data.data(:,3) == 0)];

% define the bus for which the CPF analysis is to be done
busCPF = 10;

%% Calculating the Y-bus matrix
Y_bus = Ybus(n_bus, n_branch, branch_imp, bus_imp, from, to);
G = real(Y_bus); % conductance (G) <- real part of admittance 
B = imag(Y_bus); % susceptance (B) <- the imaginary part of admittance 

%% computes the power flow
% lambda = 1;
% [V_flat, delta_flat, ~] = powerflow(tolerance, n_bus, bus_data, ...
%     base_MW, G, B,Y_bus, V_flat, delta_flat,n_pq, n_pv, pq_bus_id, ...
%     pv_bus_id, lambda*Ps, lambda*Qs);

%% PV curve - Part 1: Changing Lambda
sigma = 0.1;             
lambda = 0;  
[V_flat,delta_flat,~] = powerflow(tolerance, n_bus, bus_data,...
        base_MW, G, B, Y_bus, V_flat, delta_flat,n_pq, n_pv, pq_bus_id, ...
        pv_bus_id, lambda*Ps, lambda*Qs);
    
delta_CPF = delta_flat(bus_data.data(:,3) ~= 3);
V_CPF = V_flat(bus_data.data(:,3) == 0);
iter = 0;

while iter < 10
    %%%%%%%%%%%%%%%%%%%%%%%%%% PREDICTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % saves the last value when power flow diverges    
    plot_lambda = lambda;
    plot_V = V_flat;
    plot_delta = delta_flat; 
        
    % delta/V/lambda solution vector
    d_V_L = [delta_CPF; V_CPF; lambda];
    
    % Jacobian 
    [J] = Jacobian(V_flat, delta_flat, n_bus, n_pq, pq_bus_id, G, B, Y_bus);
    
    % ek vector -> since we are changing lambda we keep 1 at the last
    ek = [zeros(1,length(J)) 1];
    
    % augmented Jacobian
    aug_J = [J -K; ek];
    
    %inversion using crout's method
    delta_d_V_L = croutLU(aug_J, ek);
    
    % solution 
    d_V_L = d_V_L + sigma * delta_d_V_L;
    
    % extracting and assigning the respective parameters
    % delta for CPF
    delta_CPF = d_V_L(1:length(delta_CPF));
    % update original delta
    delta_flat(bus_data.data(:,3) ~= 3) = delta_CPF;
    
    % V for CPF
    V_CPF = d_V_L(length(delta_CPF) + (1:length(V_CPF)));
    % update original V
    V_flat(bus_data.data(:,3) == 0) = V_CPF;
    
    % lambda
    lambda = d_V_L(end);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% CORRECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [V_flat,delta_flat,iter] = powerflow(tolerance, n_bus, bus_data,...
        base_MW, G, B, Y_bus, V_flat, delta_flat,n_pq, n_pv, pq_bus_id, ...
        pv_bus_id, lambda*Ps, lambda*Qs);
    
    plot(plot_lambda,plot_V(busCPF),'or'); hold on;
    title(['CPF for Bus ' num2str(busCPF)])
    grid('on')  
end

%% PV curve - Part 2: Changing V
sigma = 0.005;             
lambda = plot_lambda;
V_flat = plot_V;
delta_flat = plot_delta;
Nose = 0;
change_factor=0.75;

while lambda > change_factor * Nose
    %%%%%%%%%%%%%%%%%%%%%%%%%% PREDICTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % delta/V/lambda solution vector
%     d_V_L = [delta_CPF; V_CPF; lambda];
    
    % Jacobian 
    [J] = Jacobian(V_flat, delta_flat, n_bus, n_pq, pq_bus_id, G, B,Y_bus);
    
    % ek vector -> since we are changing lambda we keep 1 at the last
    ek = [zeros(1,length(J)) 0];
    b = [zeros(1,length(J)) 1]; % used for crout's LU
    
    % put -1 to the bus on which CPF analyis is to be done
    ek(length(delta_CPF) + ek_positions(busCPF,3)) = -1;
        
    % augmented Jacobian
    aug_J = [J -K; ek];
    
    %inversion using crout's method
    delta_d_V_L = croutLU(aug_J, b);
    
    % solution 
    d_V_L = d_V_L + sigma * delta_d_V_L;
    
    [V_flat,delta_flat,lambda] = Update_Variables(sigma,delta_d_V_L,...
            V_flat,delta_flat,lambda,bus_data);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%% CORRECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mismatch = power_mismatch(lambda*Ps, lambda*Qs, G, B, V_flat,...
        delta_flat, n_bus, pq_bus_id);
    
    while max(abs(mismatch)) <= 0.01
        % Jacobian 
        [J] = Jacobian(V_flat, delta_flat, n_bus, n_pq, pq_bus_id, G,...
            B,Y_bus);
        
        % augmented Jacobian
        aug_J = [J -lambda*K; ek];
    
        % inversion using crout's method
        delta_d_V_L = croutLU(aug_J, [mismatch;0]);        
        
        [V_flat,delta_flat,lambda] = Update_Variables(sigma,delta_d_V_L,...
            V_flat,delta_flat,lambda,bus_data);
        
        % mismatch
        mismatch = power_mismatch(lambda*Ps, lambda*Qs, G, B, V_flat,...
        delta_flat, n_bus, pq_bus_id);    
    end
    
    if lambda>Nose
        Nose=lambda;
    end
    
    plot(lambda,V_flat(busCPF),'ob'); hold on;
    title(['CPF for Bus ' num2str(busCPF)])
    grid('on')    
end


%% Solving for power flow using NRPF algorithm

function [V_final,Angle_final,iter] = powerflow(tolerance, n_bus,...
    bus_data, base_MW, G, B, Y_bus, V_flat, delta_flat, n_pq, n_pv,...
    pq_bus_id, pv_bus_id, Ps, Qs)

    % Newton Rhapson Power Flow 
    [Volt, Angle, iter] = ...
        NewtonRhapson(tolerance, n_bus, n_pv, n_pq, pq_bus_id,...
        V_flat, delta_flat, G, B, Y_bus, Ps, Qs);

    V_final = Volt(:,end);
    Angle_final = Angle(:,end);
%     plot_states(Volt, Angle)
end
%% plots of the result
function plot_states(Volt, Angle)
    % NRLF Voltage
    figure('color', [1,1,1])
    str = "bus 1";
    for i = 1: length(Volt)
        plot(Volt(i,:), 'Linewidth', 1.5)
        hold on
        if i > 1
            str = [str , strcat('bus',' ', num2str(i))];
        end
    end

    ylabel('Voltage (pu)')
    xlabel('Number of iteration')
    title('NRPF')
    grid on
    set(gca,'XTick',(1:1:10))
    set(gca,'gridlinestyle','--','fontname','Times New Roman','fontsize',12);
    lgd = legend (str, 'NumColumns', 4);
    lgd.FontSize = 9;
    hold off

    % NRLF Angle
    figure('color', [1,1,1])
    for i = 1: length(Angle)
        plot(Angle(i,:), 'Linewidth', 1.5)
        hold on
    end
    ylabel('Angle (rad)')
    xlabel('Number of iteration')
    title('NRPF')
    grid on
    set(gca,'XTick',(1:1:10))
    set(gca,'gridlinestyle','--','fontname','Times New Roman','fontsize',12);
    lgd = legend (str, 'NumColumns', 3);
    lgd.FontSize = 9;
    hold off
end

function [V,theta,lambda] = Update_Variables(Step_Size,Error,V,theta,lambda,bus_data)
n=length(V);
Error=Step_Size*Error;
dtheta = Error(1:n-1);
dV = Error(n:end-1);
dlambda = Error(end);

theta(2:n) = dtheta + theta(2:n);
k = 1;
for i = 2:n
    if bus_data.data(i,3) == 0
        V(i) = dV(k) + V(i);
        k = k+1;
    end
end
lambda = lambda + dlambda;
end