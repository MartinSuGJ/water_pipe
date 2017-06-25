%% Read the data
f = fopen('by.txt');

line = fgetl(f);
size_info = textscan(line, '%d');
num_v = str2num(num2str(size_info{1}(1)));
num_e = str2num(num2str(size_info{1}(2)));

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
f = fopen('by_edge.txt');
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
    %edge_type = [edge_type; line{7}(row)];

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
f = fopen('by_node.txt');
line = textscan(f, '%d %f %f %d');

d = zeros(num_v,1);
hc = zeros(num_v,1);
node_type = []; % 0:normal; 1: customer; 2: source; 3: tank.

for row = 1:num_v,
    node_id = line{1}(row);
    demand = line{2}(row);
    head = line{3}(row);
    % node_type = [node_type; line{4}(row)];
    d(node_id) = -demand/1000;
    hc(node_id) = head;
end

d(source) = 1000;


%% Predirection
cvx_solver mosek
cvx_begin
    variables q(num_e*2);
    minimize ( (1/3)*[L,L]*pow_p(q,3) )
    % minimize ( (1/3)*sum(pow_p(q,3)) )
    [A,-A]*q <= d
    q >= 0
cvx_end

% change the direction of the flow
for i = 1 : num_e
  if q(i) < q(i+num_e)
    A(:,i) = A(:,i).*-1;
  end
end


%% Add Imaginary Nodes and Edges
% Add imaginary nodes
Ai = [A; zeros(1, num_e)]; 

% Add imaginary edges
Ei = zeros(size(Ai,1), num_v);
for row = 1:num_v,
    t = row;
    im = num_v + 1;
    
    Ei(t, row) = 1;
    Ei(im, row) = -1;
end

%% Add Imaginary Lambda and Demand
epsilon = 1e-6;
dummy_lambda = hc./epsilon^2;
dummy_lambda(terminal) = hc(terminal)./(d(terminal).^2);

di = [d; sum(d(terminal))-(num_v-num_t)*epsilon];
di(terminal) = 0;

qiL = ones(num_v,1)*epsilon;
qiL(terminal) = -d(terminal);


%% Solve the imaginary model
cvx_solver mosek
cvx_begin
    variables q(num_e + num_v)
    minimize ( 100*(1/3)*[L]*pow_p(q(1:num_e),3) + (1e-8)*dummy_lambda'*q(num_e+1:end) )
    [Ai, Ei]*q <= di
    q(1:num_e) >= 0
    q(num_e+1:end) >= qiL   
cvx_end

% Solve for the head pressure by fixing the flow
cvx_begin
    variables h(num_v + 1);
    % minimize ( norm([Ai, Ei]'*h - diag([L, dummy_lambda'])*pow_p(q,2)) )    
    minimize ( 100*norm(A'*h(1:end-1) - diag(L)*q(1:num_e).^2) + q(1:num_e)'*A'*h(1:end-1) + norm(h) )
    [Ai, Ei]'*h - diag([L, dummy_lambda'])*pow_p(q,2) >= 0
    h(end) == 0
cvx_end

%% Model Performance
energy = q(1:num_e)'*A'*h(1:end-1)
gap = norm(A'*h(1:end-1) - diag(L)*q(1:num_e).^2)

