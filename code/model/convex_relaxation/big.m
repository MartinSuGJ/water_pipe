clc
clear

f = fopen('big.txt');

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
f = fopen('big_edge.txt');
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
f = fopen('big_node.txt');
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
head_list = [];
tail_list = [];
for i = 1 : num_e
  if q(i) < q(i+num_e)
    A(:,i) = A(:,i).*-1;
  end
  head_list = [head_list; find(A(:,i)==1)];
  tail_list = [tail_list; find(A(:,i)==-1)];
end

cvx_begin
    variables q(num_e);
    minimize ( (1/3)*L*pow_p(q,3) )
    % minimize ( (1/3)*sum(pow_p(q,3)) )
    A*q <= d
    q >= 0
cvx_end

%% Solve the convexation problem
cvx_solver mosek
cvx_save_prefs


cvx_begin
  variables q(num_e) h(num_v) s(num_e);
  minimize ((1/3)*L*pow_p(q,3) + (1/400)*norm(h) + (1/150)*ones(num_e,1)'*s) % big weight
  % minimize ((1/3)*L*pow_p(q,3) + norm(h) + ones(num_e,1)'*s) % original
  %minimize ((1/3)*L_new*pow_p(q_new,3) + 1/80*norm(h_new) + 1/80*ones(num_e_new,1)'*s_new) % small weight
  A*q <= d
  q >= 0
  h >= hc
  A'*h - s == 0
  diag(L)*pow_p(q,2) <= s
cvx_end

original_gap = norm(A'*h - diag(L)*pow_p(q,2));
original_opt_val = cvx_optval;
original_resource = (1/3)*L*pow_p(q,3);
original_norm_h = norm(h);
original_sum_s = ones(num_e,1)'*s;

%% Write the profile of each edge
threshold = 0.05;
threshold_q = 0.00001;
head_diff = A'*h;
energy_use = (1/3)*diag(L)*pow_p(q,3);
gap = A'*h - diag(L)*pow_p(q,2);
edge_with_s = find(abs(gap)<=threshold);
edge_with_ns_zero_q = find(abs(gap)>threshold & q<=threshold_q);
edge_with_ns_q = find(abs(gap)>threshold & q>threshold_q);
edge_gap_profile = zeros(num_e, 1);
edge_gap_profile(edge_with_s) = 0;
edge_gap_profile(edge_with_ns_zero_q) = 1;
edge_gap_profile(edge_with_ns_q) = 2;

edge_profile = [(1:num_e)', L', q*10000000, head_list, tail_list, h(head_list)*1000, ...
                hc(head_list)*1000, h(tail_list)*1000, hc(tail_list)*1000, gap*1000, edge_type, ...
                edge_gap_profile];

filename = 'edge_profile.xlsx';
xlswrite(filename,edge_profile);



%% Remove the pipe with NS and q==0 AND remove the valve with NS and q==0
A_backup = A;
L_backup = L;
difference = A'*h - diag(L)*pow_p(q,2);
pipe_with_ns_zero_q = find(abs(difference)>threshold & edge_type == 0 & q<= threshold_q);
valve_with_ns_zero_q = find(abs(difference)>threshold & edge_type == 2 & q<= threshold_q);
remove_ns_zero_q = [pipe_with_ns_zero_q;valve_with_ns_zero_q];
A(:,remove_ns_zero_q) = [];
L(remove_ns_zero_q) = [];

size_A = size(A);

cvx_begin
  variables q(size_A(2)) h(size_A(1)) s(size_A(2));
  minimize ((1/3)*L*pow_p(q,3) + (1/400)*norm(h) + (1/150)*ones(size_A(2),1)'*s) % big weight
  % minimize ((1/3)*L*pow_p(q,3) + norm(h) + ones(num_e,1)'*s) % original
  %minimize ((1/3)*L_new*pow_p(q_new,3) + 1/80*norm(h_new) + 1/80*ones(num_e_new,1)'*s_new) % small weight
  A*q <= d
  q >= 0
  h >= hc
  A'*h - s == 0
  diag(L)*pow_p(q,2) <= s
cvx_end

remove_gap = norm(A'*h - diag(L)*pow_p(q,2));
remove_opt_val = cvx_optval;
remove_resource = (1/3)*L*pow_p(q,3);
remove_norm_h = norm(h);
remove_sum_s = ones(size_A(2),1)'*s;

edge_type(remove_ns_zero_q) = [];
difference1 = A'*h - diag(L)*pow_p(q,2);
pipe_with_ns_zero_q1 = find(abs(difference1)>threshold & edge_type == 0 & q<= threshold_q);
valve_with_ns_zero_q1 = find(abs(difference1)>threshold & edge_type == 2 & q<= threshold_q);
remove_ns_zero_q1 = [pipe_with_ns_zero_q1;valve_with_ns_zero_q1];


%% Remove the pipe with NS and q>0
pipe_with_ns_q = find(abs(difference1)>threshold & edge_type == 0 & q > threshold_q);
valve_with_ns_q = find(abs(difference1)>threshold & edge_type == 2 & q > threshold_q);
remove_ns_q = [pipe_with_ns_q; valve_with_ns_q];
[a,b] = sort(difference1(remove_ns_q), 'descend');
remove_ns_q = remove_ns_q(b);
size_remove_ns_q = size(remove_ns_q);

potential_list = [];
gaps = [];
resources = [];
solvable_profile = [];
for i=1:size_remove_ns_q(1)
    A_R = A;
    % A_R(:, [potential_list;remove_ns_q(i)]) = [];
    A_R(:, remove_ns_q(i)) = [];
    L_R = L;
    L_R([potential_list;remove_ns_q(i)]) = [];
    
    size_A_R = size(A_R);
 
    cvx_begin
      variables q(size_A_R(2)) h(size_A_R(1)) s(size_A_R(2));
      minimize ((1/3)*L_R*pow_p(q,3) + (1/400)*norm(h) + (1/150)*ones(size_A_R(2),1)'*s) % big weight
      % minimize ((1/3)*L*pow_p(q,3) + norm(h) + ones(num_e,1)'*s) % original
      %minimize ((1/3)*L_new*pow_p(q_new,3) + 1/80*norm(h_new) + 1/80*ones(num_e_new,1)'*s_new) % small weight
      A_R*q <= d
      q >= 0
      h >= hc
      A_R'*h - s == 0
      diag(L_R)*pow_p(q,2) <= s
    cvx_end
    
    if sum(size(findstr(cvx_status, 'Solved'))) > 0
        solvable_profile(i) = 1;
        gap = norm(A_R'*h - diag(L_R)*pow_p(q,2));
        resource = (1/3)*L_R*pow_p(q,3);
        gaps(i) = (gap-remove_gap)/remove_gap;
        resources(i) = (resource-remove_resource)/remove_resource;
%         if (gap < original_gap) & (resource < original_resource)
%             potential_list = [potential_list;remove_ns_q(i)];
%             gaps = [gaps;gap];
%             resources = [resources;resource];
%         end
    else
        solvable_profile(i) = 0;
        gaps(i) = nan;
        resources(i) = nan;
    end
end

%%

% 
% i = 5;
% 
% A_R = A;
% A_R(:, remove_ns_q(i)) = [];
% L_R = L;
% L_R(remove_ns_q(i)) = [];
% 
% size_A_R = size(A_R);
% 
% cvx_begin
%   variables q(size_A_R(2)) h(size_A_R(1)) s(size_A_R(2));
%   minimize ((1/3)*L_R*pow_p(q,3) + (1/400)*norm(h) + (1/150)*ones(size_A_R(2),1)'*s) % big weight
%   % minimize ((1/3)*L*pow_p(q,3) + norm(h) + ones(num_e,1)'*s) % original
%   %minimize ((1/3)*L_new*pow_p(q_new,3) + 1/80*norm(h_new) + 1/80*ones(num_e_new,1)'*s_new) % small weight
%   A_R*q <= d
%   q >= 0
%   h >= hc
%   A_R'*h - s == 0
%   diag(L_R)*pow_p(q,2) <= s
% cvx_end
% 
% 
% %%
%   threshold = 0.05;
%   threshold_q = 0.00001;
%   statistics(1, 1:end) = [numel(find(abs(difference)<=threshold & edge_type == 0)), ... 
%                      numel(find(abs(difference)>threshold & edge_type == 0 & q<= threshold_q)), ...
%                      numel(find(abs(difference)>threshold & edge_type == 0 & q>threshold_q)), ...
%                      norm(difference(find(abs(difference)<=threshold & edge_type == 0))), ...
%                      norm(difference(find(abs(difference)>threshold & edge_type == 0 & q<= threshold_q))), ...
%                      norm(difference(find(abs(difference)>threshold & edge_type == 0 & q>threshold_q))), ...
%                      norm(difference(find(edge_type==0))), ...
%                      numel(find(abs(difference)<=threshold & edge_type == 1)), ...
%                      numel(find(abs(difference)>threshold & edge_type == 1 & q<= threshold_q)), ...
%                      numel(find(abs(difference)>threshold & edge_type == 1 & q>threshold_q)), ...
%                      norm(difference(find(abs(difference)<=threshold & edge_type == 1))), ...
%                      norm(difference(find(abs(difference)>threshold & edge_type == 1 & q<= threshold_q))), ...
%                      norm(difference(find(abs(difference)>threshold & edge_type == 1 & q>threshold_q))), ...
%                      norm(difference(find(edge_type==1))), ...
%                      numel(find(abs(difference)<=threshold & edge_type == 2)), ...
%                      numel(find(abs(difference)>threshold & edge_type == 2 & q<= threshold_q)), ...
%                      numel(find(abs(difference)>threshold & edge_type == 2 & q>threshold_q)),...
%                      norm(difference(find(abs(difference)<=threshold & edge_type == 2))), ...
%                      norm(difference(find(abs(difference)>threshold & edge_type == 2 & q<= threshold_q))), ...
%                      norm(difference(find(abs(difference)>threshold & edge_type == 2 & q>threshold_q))), ...
%                      norm(difference(find(edge_type==2)))];
% 
%   gap(ind) = norm(A'*h - diag(L)*pow_p(q,2));
%   opt_val(ind) = cvx_optval;
%   resource(ind) = (1/3)*L*pow_p(q,3);
%   norm_h(ind) = norm(h);
%   sum_s(ind) = ones(num_e,1)'*s;
%   ind = ind + 1;
% end
% 
figure(1);
yyaxis left;
hold on;
title('Energy Consumption and Gap Change as Removing Edges');
hold on;
ylabel('Energy Consumption');
hold on;
xlabel('Number of Edges Being Removed');
hold on;
plot(1:23, resources);
hold on;
yyaxis right;
hold on;
ylabel('Gap');
hold on;
plot(1:23, gaps);
