clc
clear

load('big.mat')

%% Predirection
disp('Predirection...')

cvx_begin
    variables q(num_e*2);
    minimize ( (1/3)*[L,L]*pow_p(q,3) )    
    [A,-A]*q <= d
    q >= 0
cvx_end
energy_lowerbound = cvx_optval;

% change the direction of the flow
idx = q(1:num_e) < q(num_e+1:end);
A(:,idx) = -A(:,idx);

%% Convex Relaxation
cvx_begin 
    variables q(num_e) h(num_v) s(num_e);
    minimize ((1/3)*L*pow_p(q,3) + norm(h)/1000 + ones(num_e,1)'*s/100)
    A*q <= d
    q >= 0
    h >= hc
    A'*h - s == 0
    s >= diag(L)*pow_p(q,2)
cvx_end

% Model Performance
energy = q'*A'*h
energy_lowerbound
gap = norm(A'*h - diag(L)*q.^2)
min_press_constraint = sum(h >= hc)/size(hc,1)
demand_constraint = sum(A*q <= d)/size(d,1)