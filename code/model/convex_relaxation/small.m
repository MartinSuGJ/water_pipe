clc
clear

f = fopen('small.txt');

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

cvx_begin
    variables q(num_e);
    minimize ( (1/3)*L*pow_p(q,3) )
    A*q <= d
    q >= 0
cvx_end

%%
cvx_solver mosek
cvx_save_prefs


ind = 1;
statistics = zeros(19,21);
for (i = 1:0.5:10)
  cvx_begin
      variables q(num_e) h(num_v) s(num_e);
      %minimize (10*(1/3)*L*pow_p(q,3) + 1/400*norm(h) + 1/150*ones(num_e,1)'*s) % big weight
      % minimize ((1/3)*L*pow_p(q,3) + norm(h) + ones(num_e,1)'*s) % original
      minimize ((1/3)*L*pow_p(q,3) + (1/80)*norm(h) + i*(2/80)*ones(num_e,1)'*s) % big weight
      %minimize ((1/3)*L_new*pow_p(q_new,3) + 1/80*norm(h_new) + 1/80*ones(num_e_new,1)'*s_new) % small weight
      A*q <= d
      q >= 0
      h >= hc
      A'*h - s == 0
      diag(L)*pow_p(q,2) <= s
  cvx_end
  
  difference = A'*h - diag(L)*pow_p(q,2);
  
  threshold = 0.01;
  threshold_q = 0.005;
  statistics(ind, 1:end) = [numel(find(abs(difference)<=threshold & edge_type == 0)), ... 
                     numel(find(abs(difference)>threshold & edge_type == 0 & q<= threshold_q)), ...
                     numel(find(abs(difference)>threshold & edge_type == 0 & q>threshold_q)), ...
                     norm(difference(find(abs(difference)<=threshold & edge_type == 0))), ...
                     norm(difference(find(abs(difference)>threshold & edge_type == 0 & q<= threshold_q))), ...
                     norm(difference(find(abs(difference)>threshold & edge_type == 0 & q>threshold_q))), ...
                     norm(difference(find(edge_type==0))), ...
                     numel(find(abs(difference)<=threshold & edge_type == 1)), ...
                     numel(find(abs(difference)>threshold & edge_type == 1 & q<= threshold_q)), ...
                     numel(find(abs(difference)>threshold & edge_type == 1 & q>threshold_q)), ...
                     norm(difference(find(abs(difference)<=threshold & edge_type == 1))), ...
                     norm(difference(find(abs(difference)>threshold & edge_type == 1 & q<= threshold_q))), ...
                     norm(difference(find(abs(difference)>threshold & edge_type == 1 & q>threshold_q))), ...
                     norm(difference(find(edge_type==1))), ...
                     numel(find(abs(difference)<=threshold & edge_type == 2)), ...
                     numel(find(abs(difference)>threshold & edge_type == 2 & q<= threshold_q)), ...
                     numel(find(abs(difference)>threshold & edge_type == 2 & q>threshold_q)),...
                     norm(difference(find(abs(difference)<=threshold & edge_type == 2))), ...
                     norm(difference(find(abs(difference)>threshold & edge_type == 2 & q<= threshold_q))), ...
                     norm(difference(find(abs(difference)>threshold & edge_type == 2 & q>threshold_q))), ...
                     norm(difference(find(edge_type==2)))];

  gap(ind) = norm(A'*h - diag(L)*pow_p(q,2));
  opt_val(ind) = cvx_optval;
  resource(ind) = (1/3)*L*pow_p(q,3);
  norm_h(ind) = norm(h);
  sum_s(ind) = ones(num_e,1)'*s;
  ind = ind + 1;
end

figure(1);
yyaxis left;
hold on;
title('Energy Consumption and Gap Change as Weight Changes');
hold on;
ylabel('Energy Consumption');
hold on;
xlabel('Weight on Energy Consumption');
hold on;
plot(1:0.5:10, resource);
hold on;
yyaxis right;
hold on;
ylabel('Gap');
hold on;
plot(1:0.5:10, gap);
