clc
clear

%% Read the data
f = fopen('small.txt');

line = fgetl(f);
size_info = textscan(line, '%d');
num_v = str2num(num2str(size_info{1}(1)));
num_e = str2num(num2str(size_info{1}(2)));
num_p = str2num(num2str(size_info{1}(3)));

line = fgetl(f);
source_info = textscan(line, '%d');
source = source_info{1};
num_s = size(source, 1);

line = fgetl(f);
terminal_info = textscan(line, '%d');
terminal = terminal_info{1};
num_t = size(terminal, 1);

fclose(f);

%% Construct A and Lambda
f = fopen('small_edge.txt');
line = textscan(f, '%d %d %d %f %f %f %d');

E = {};
A = [];
L = [];
g = 9.8;

A = zeros(num_v, num_e);
edge_type = []; % 0: pipe; 1: pump; 2: valve.
for row = 1:num_e,
    i = line{2}(row);
    j = line{3}(row);
    length = line{4}(row);
    diameter = line{5}(row)/1000;
    roughness = line{6}(row)/1000;
    edge_type = [edge_type; line{7}(row)];

    friction = (2 * log10(roughness/(3.71*diameter)))^2;
    lambda = (8 * length * friction)/(pi^2*g*diameter^5);
    L = [L lambda];

    A(i, row) = 1;
    A(j, row) = -1;

    e = strcat(num2str(i), ',', num2str(j));
    E = [E; e];
end
size(A)
fclose(f);

%% Construct hc and d
f = fopen('small_node.txt');
line = textscan(f, '%d %f %f %d');

d = zeros(num_v,1);
hc = zeros(num_v,1);
node_type = []; % 0:normal; 1: customer; 2: source; 3: tank.

for row = 1:num_v,
    node_id = line{1}(row);
    demand = line{2}(row);
    head = line{3}(row);
    node_type = [node_type; line{4}(row)];
    d(node_id) = -demand/1000;
    hc(node_id) = head;
end

d(source) = 1000;

%% Fit the Pump Curve
f = fopen('small_pump.txt');
line = textscan(f, '%d %d %d %f %f');
num_line = size(line{1},1);
pump_info = {};
pump_curve = {};
pump_coeff = {};

for i = 1:num_p,
    pump_info{i} = [];
    pump_curve{i} = [];
end

for i = 1:num_line,
    p_id = line{1}(i);
    if p_id <= 20
        pump_info{p_id} = [line{2}(i), line{3}(i)];
        pump_curve{p_id} = [pump_curve{p_id};[line{4}(i)/100,line{5}(i)*100]];
    else

        pump_info{p_id} = [line{2}(i), line{3}(i)];
        pump_curve{p_id} = [pump_curve{p_id};[line{4}(i)/1000,line{5}(i)*100]];        
    end
end

for i = 1:num_p,
    p = polyfit(pump_curve{i}(:,1).^2,pump_curve{i}(:,2),1);
    pump_coeff{i} = p;
end


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

%% Change the L vector to fit the pump curve model
dh_max = zeros(num_e, 1);
L1 = L;
pump_edge_list = find(edge_type==1);
pump_head_list = [];
pump_tail_list = [];
for i = 1:num_p,
    pump_source = find(A(:,pump_edge_list(i))==1);
    pump_tail = find(A(:,pump_edge_list(i))==-1);
    pump_head_list = [pump_head_list,pump_source];
    pump_tail_list = [pump_tail_list,pump_tail];
    for j = 1:num_p,
        if pump_source == pump_info{j}(1) || pump_source == pump_info{j}(1),
            dh_max(pump_edge_list(i)) = pump_coeff{j}(2);
            L1(pump_edge_list(i)) = -pump_coeff{j}(1);
        end
    end
end

%% Flip the sign for L1
L2 = L1;
L2(pump_edge_list) = -L2(pump_edge_list);

%% Solve the flow
disp('Solving for the flow...')

cvx_begin
    variables q(num_e + num_v)
%     minimize ( (L1.*(edge_type'~=1))*pow_p(q(1:num_e),3) )
    minimize ( (L1)*pow_p(q(1:num_e),3) - q(pump_edge_list)'*dh_max(pump_edge_list) )
%     minimize ( (L1.*(edge_type'~=1))*pow_p(q(1:num_e),3) )
    [Ai, Ei]*q <= di
    q(1:num_e) >= 0
    q(num_e+1:end) >= qiL
cvx_end

%% Solve pressure fixing the flow
disp('Solving for the pressure fixing the flow...')

cvx_begin
    variables h(num_v + 1);
    minimize ( norm(A'*h(1:end-1) - diag(L1)*q(1:num_e).^2 + dh_max) )
%     minimize ( norm( (A'*h(1:end-1) - diag(L1)*q(1:num_e).^2 + dh_max).*(edge_type==1) ) )
    h(1:num_v) >= hc
    h(pump_head_list) == hc(pump_head_list)
    h(end) == 0
%     A'*h(1:end-1) - diag(L1)*q(1:num_e).^2 + dh_max >= 0
    A'*h(1:end-1).*(edge_type~=1) - (diag(L1)*q(1:num_e).^2).*(edge_type~=1) >= 0
    A'*h(1:end-1).*(edge_type==1) - (diag(L1)*q(1:num_e).^2).*(edge_type==1) + dh_max.*(edge_type==1) <= 0
cvx_end

%% Model Performance
energy_pump = q(pump_edge_list)'*dh_max(pump_edge_list) - L1(pump_edge_list)*pow_p(q(pump_edge_list),3);
energy_loss = L1*pow_p(q(1:num_e),3) - L1(pump_edge_list)*pow_p(q(pump_edge_list),3)
energy_lowerbound
gap = norm(A'*h(1:end-1) - diag(L1)*q(1:num_e).^2 + dh_max)
gap_pump = norm(A'*h(1:end-1).*(edge_type==1) - (diag(L1)*q(1:num_e).^2).*(edge_type==1) + dh_max.*(edge_type==1));
gap_pipe = norm(A'*h(1:end-1).*(edge_type~=1) - (diag(L1)*q(1:num_e).^2).*(edge_type~=1));
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
%         minimize norm( A'*h - diag(L)*pow_p(q,2) )
        minimize ( norm(A'*h - diag(L1)*q(1:num_e).^2 + dh_max) )
        h >= hc
        h(pump_head_list) == hc(pump_head_list)
%         A'*h - diag(L1)*q(1:num_e).^2 + dh_max >= 0
        A'*h.*(edge_type~=1) - (diag(L1)*q(1:num_e).^2).*(edge_type~=1) >= 0
        A'*h.*(edge_type==1) - (diag(L1)*q(1:num_e).^2).*(edge_type==1) + dh_max.*(edge_type==1) <= 0
    cvx_end
    
    energy = [energy; (L1.*(edge_type'~=1))*pow_p(q(1:num_e),3)];
    gap = [gap; norm(A'*h - diag(L1)*q(1:num_e).^2 + dh_max)];
    normH = [normH; norm(h)];


    % Step 2: Maximum flow adjustment
    qU = sqrt((A'*h + dh_max)./L1');

    cvx_begin
        variables q(num_e)
        minimize -(sourceList'*(A*q))
        A*q <= d
        q >= 0
        q <= qU
%         A'*h - diag(L1)*q(1:num_e).^2 + dh_max >= 0
        A'*h.*(edge_type~=1) - (diag(L1)*q(1:num_e).^2).*(edge_type~=1) >= 0
        A'*h.*(edge_type==1) - (diag(L1)*q(1:num_e).^2).*(edge_type==1) + dh_max.*(edge_type==1) <= 0
    cvx_end

    energy = [energy; (L1.*(edge_type'~=1))*pow_p(q(1:num_e),3)];
    gap = [gap; norm(A'*h - diag(L1)*q(1:num_e).^2 + dh_max)];
    normH = [normH; norm(h)];
    
    iter = iter + 1;
    if iter >= 50
        break;
    end
end

% Step 1: pressure adjustment
cvx_begin
    variables h(num_v)
%         minimize norm( A'*h - diag(L)*q.^2) + (1/1000)*sum((h(find(d<0)) - hc((find(d<0))).*(-d((find(d<0))))))
%         minimize norm( A'*h - diag(L)*pow_p(q,2) )
    minimize ( norm(A'*h - diag(L1)*q(1:num_e).^2 + dh_max) )
    h >= hc
    h(pump_head_list) == hc(pump_head_list)
%     A'*h - diag(L1)*q(1:num_e).^2 + dh_max >= 0
    A'*h.*(edge_type~=1) - (diag(L1)*q(1:num_e).^2).*(edge_type~=1) >= 0
    A'*h.*(edge_type==1) - (diag(L1)*q(1:num_e).^2).*(edge_type==1) + dh_max.*(edge_type==1) <= 0
cvx_end

energy = [energy; (L1.*(edge_type'~=1))*pow_p(q(1:num_e),3)];
gap = [gap; norm(A'*h - diag(L1)*q(1:num_e).^2 + dh_max)];
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