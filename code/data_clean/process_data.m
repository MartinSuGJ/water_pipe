clc
clear

%% Read the data
f = fopen('big.txt');

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

%% Fit the Pump Curve
f = fopen('big_pump.txt');
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
    pump_info{p_id} = [line{2}(i), line{3}(i)];
    pump_curve{p_id} = [pump_curve{p_id};[line{4}(i),line{5}(i)]];
end

for i = 1:num_p,
    p = polyfit(pump_curve{i}(:,1).^2,pump_curve{i}(:,2),1);
    pump_coeff{i} = p;
end

% pump curve representation
for i = 1:num_p,
    x = pump_curve{i}(:,1);
    y = pump_curve{i}(:,2);
    p = pump_coeff{i};
    p = [p(1), 0.0, p(2)];
    x1 = linspace(min(x), max(x));
    y1 = polyval(p,x1);
    figure
    plot(x,y,'o')
    hold on
    plot(x1,y1)
    hold off
    title(['Pump Curve ' num2str(i)]);
end
