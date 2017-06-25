clc
clear

%% Read the data
f = fopen('town.txt');

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
f = fopen('town_edge.txt');
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
f = fopen('town_node.txt');
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

%% Predirection part
cvx_solver mosek
cvx_begin
    variables q(num_e*2);
    minimize ( (1/3)*[L,L]*pow_p(q,3) )
    % minimize ( (1/3)*sum(pow_p(q,3)) )
    [A,-A]*q <= d
    q >= 0
cvx_end

% for i = 1 : num_e
%   fprintf('%f, %f\n',[q(i),q(i+num_e)])
% end

% change the direction of the flow
for i = 1 : num_e
  if q(i) < q(i+num_e)
    A(:,i) = A(:,i).*-1;
  end
end

% re-solve the problem with the fixed direction
cvx_begin
    variables q(num_e);
    minimize ( (1/3)*L*pow_p(q,3) )
    A*q <= d
    q >= 0
cvx_end





%% Convex relaxation
cvx_solver mosek
cvx_save_prefs

cvx_begin
  variables q(num_e) h(num_v) s(num_e);
  % minimize ((1/3)*L*pow_p(q,3) + (1/80)*norm(h) + (2/80)*ones(num_e,1)'*s) % small network weight
  minimize ((1/3)*L*pow_p(q,3) + (1/400)*norm(h) + (1/150)*ones(num_e,1)'*s) % big network weight 
  A*q <= d
  q >= 0
  h >= hc
  A'*h - s == 0
  diag(L)*pow_p(q,2) <= s
cvx_end


%% Iterative method to adjust the pressure and flow
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
        minimize norm( A'*h - diag(L)*q.^2) + (1/1000)*sum((h(find(d<0)) - hc((find(d<0))).*(-d((find(d<0))))))
        % minimize norm( A'*h - diag(L)*pow_p(q,2) ) + (q'*A')*h + (1/1000)*norm(h)
        h >= hc
        A'*h >= diag(L)*q.^2
    cvx_end
    
    energy = [energy; q'*A'*h];
    gap = [gap; norm( A'*h - diag(L)*q.^2)];
    normH = [normH; norm(h)];


    % Step 2: Maximum flow adjustment
    qU = sqrt(diag(L)\(A'*h));

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
    if iter >= 20
        break;
    end
end




