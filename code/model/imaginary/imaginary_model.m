clc
clear

load('town_1.mat')

%% Predirection
disp('Predirection...')

cvx_begin
    variables q(num_e*2);
    minimize ( (1/3)*[L,L]*pow_p(q,3) )
    [A,-A]*q <= d
    q >= 0
cvx_end

% change the direction of the flow
idx = q(1:num_e) < q(num_e+1:end);
A(:,idx) = -A(:,idx);


cvx_begin
    variables q(num_e);
    minimize ( (1/3)*L*pow_p(q,3) )
    A*q <= d
    q >= 0
cvx_end
energy_lowerbound = cvx_optval;

%% Add Imaginary Nodes and Edges
disp('Adding imaginary nodes and edges...')

% Add imaginary nodes
Ai = [A; zeros(1, num_e)];

% Add imaginary edges
Ei = [eye(num_v); -ones(1, num_v)];

% Add imaginary lambda and demand
epsilon = 0;
dummy_lambda = hc./epsilon^2;
dummy_lambda(terminal) = hc(terminal)./(d(terminal).^2);

di = [d; 0];
di(terminal) = 0;

qiL = ones(num_v,1)*epsilon;
qiL(terminal) = -d(terminal);

%% Solve the flow
disp('Solving for the flow...')

cvx_begin
    variables q(num_e + num_v)
    minimize ( (1/3)*L*pow_p(q(1:num_e),3) )
    [Ai, Ei]*q <= di
    q(1:num_e) >= 0
    q(num_e+1:end) >= qiL
cvx_end
% cvx_begin
%     variables q(num_e);
%     minimize ( (1/3)*[L]*pow_p(q,3) )
%     A*q <= d
%     q >= 0
% cvx_end
%% Solve pressure fixing the flow
disp('Solving for the pressure fixing the flow...')

cvx_begin
    variables h(num_v + 1);
    minimize ( norm(A'*h(1:end-1) - diag(L)*q(1:num_e).^2))
    h(1:num_v) >= hc
    h(end) == 0
    A'*h(1:end-1) - diag(L)*q(1:num_e).^2 >= 0
cvx_end

%% Model Performance
energy = q(1:num_e)'*A'*h(1:end-1)
energy_lowerbound
gap = norm(A'*h(1:end-1) - diag(L)*q(1:num_e).^2)
min_press_constraint = sum(h(1:end-1) >= hc)/size(hc,1)
demand_constraint = sum(A*q(1:num_e) <= d)/size(d,1)

%% Use the iterative method to decrease the gap
q = q(1:num_e);
h = h(1:end-1);
iter = 0;
energy = [];
gap = [];
normH = [];

sourceList = zeros(num_v, 1);
sourceList(find(d>0)) = 1;
endList = zeros(num_v, 1);
endList(find(d<0)) = 1;
while 1
    % Step 1: pressure adjustment
    cvx_begin
        variables h(num_v)
%         minimize norm( A'*h - diag(L)*q.^2) + (1/1000)*sum((h(find(d<0)) - hc((find(d<0))).*(-d((find(d<0))))))
        minimize norm( A'*h - diag(L)*pow_p(q,2) )
        h >= hc
        A'*h >= diag(L)*q.^2
    cvx_end

    energy = [energy; q'*A'*h];
    gap = [gap; norm( A'*h - diag(L)*q.^2)];
    normH = [normH; norm(h)];


    % Step 2: Maximum flow adjustment
    qU = sqrt((A'*h)./L');

    cvx_begin
        variables q(num_e)
        minimize -(sourceList'*(A*q))
        A*q <= d
        q >= 0
        q <= qU
        A'*h >= diag(L)*pow_p(q,2)
    cvx_end

    energy = [energy; q'*A'*h];
    gap = [gap; norm( A'*h - diag(L)*q.^2)];
    normH = [normH; norm(h)];

    iter = iter + 1;
    if iter >= 70
        break;
    end
end


cvx_begin
    variables h(num_v)
    minimize norm( A'*h - diag(L)*pow_p(q,2) )
    h >= hc
    A'*h >= diag(L)*q.^2
cvx_end

energy = [energy; q'*A'*h];
gap = [gap; norm( A'*h - diag(L)*q.^2)];
normH = [normH; norm(h)];


%% Performance
% q_threshold = 1e-5;
% gap_threshold = 1e-5;
%
% difference = A'*h - diag(L)*q.^2;

size(difference(edge_type==0 & difference <= 1e-4))
norm(difference(edge_type==0 & difference <= 1e-4))

size(difference(difference <= 1e-4))
norm(difference(difference <= 1e-4))




size(difference(edge_type==0 & difference > 1e-4 & q>1e-5))
norm(difference(edge_type==0 & difference > 1e-4 & q>1e-5))



size(difference(edge_type==1 & difference > 1e-4 & q>1e-5))
norm(difference(edge_type==1 & difference > 1e-4 & q>1e-5))


size(difference(edge_type==2 & difference > 1e-4 & q>1e-5))
norm(difference(edge_type==2 & difference > 1e-4 & q>1e-5))


size(difference(difference > 1e-4 & q>1e-5))
norm(difference(difference > 1e-4 & q>1e-5))



size(difference(edge_type==2))
norm(difference(edge_type==2))

size(difference)
norm(difference)
