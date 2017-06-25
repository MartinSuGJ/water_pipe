clc
clear

f = fopen('shamir.txt');

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
f = fopen('shamir_edge.txt');
line = textscan(f, '%d %d %d %f %f %f');

E = {};
A = [];
L = [];
g = 9.8;

A = zeros(num_v, num_e);
for row = 1:num_e,
    i = line{2}(row);
    j = line{3}(row);
    length = line{4}(row);
    diameter = line{5}(row)/1000;
    roughness = line{6}(row)/1000;

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
f = fopen('shamir_node.txt');
line = textscan(f, '%d %f %f');

d = zeros(num_v,1);
hc = zeros(num_v,1);

for row = 1:num_v,
    node_id = line{1}(row);
    demand = line{2}(row);
    head = line{3}(row);
    d(node_id) = -demand/1000;
    hc(node_id) = head;
end

d(source) = 1000;

%% Solve
%cvx_solver SeDuMi
cvx_solver mosek
cvx_begin
    variables q(num_e*2);
    minimize ( (1/3)*[L,L]*pow_p(q,3) )
    % minimize ( (1/3)*sum(pow_p(q,3)) )
    [A,-A]*q <= d
    q >= 0
cvx_end

for i = 1 : num_e
  fprintf('%f, %f\n',[q(i),q(i+num_e)])
end

% change the direction of the flow
for i = 1 : num_e
  if q(i) < q(i+num_e)
    A(:,i) = A(:,i).*-1;
  end
end
%%
cvx_solver mosek
cvx_save_prefs

cvx_begin
    variables q(num_e) h(num_v) s(num_e);
    minimize ((1/3)*L*pow_p(q,3) + 1/25*norm(h) + 1/5*ones(num_e,1)'*s) % big weight
    % minimize ((1/3)*L*pow_p(q,3) + norm(h) + ones(num_e,1)'*s) % original
    %minimize ((1/3)*L_new*pow_p(q_new,3) + 1/80*norm(h_new) + 1/80*ones(num_e_new,1)'*s_new) % small weight
    A*q <= d
    q >= 0
    h >= hc
    A'*h - s == 0
    diag(L)*pow_p(q,2) <= s
cvx_end

gap = norm(A'*h - diag(L)*pow_p(q,2))
opt_val = cvx_optval
resource = (1/3)*L*pow_p(q,3)
norm_h = norm(h)
sum_s = ones(num_e,1)'*s
