clc
clear
f = textread('by1.txt');

num_v = f(1,1);
num_e = f(1,2);

s = f(2,:);
s = s(s > 0);
num_s = size(s,2)

% Construct A and Lambda
E = {};
A = [];
L = [];
g = 9.8;

for row = 3:3 + num_e - 1,
    pipe = f(row,2:6);
    i = pipe(1);
    j = pipe(2);
    length = pipe(3);
    diameter = pipe(4)/1000;
    roughness = pipe(5)/1000;
    
    friction = (2 * log10(roughness/(3.71*diameter)))^2;
    lambda = (8 * length * friction)/(pi^2*g*diameter^5);
    L = [L lambda];
    
    vec = zeros(num_v, 1);
    vec(i) = 1;
    vec(j) = -1;
    A = [A vec];
    
    e = strcat(num2str(i), ',', num2str(j));
    E = [E; e];
end
size(A)

d = zeros(num_v,1);

for row = 3:3 + num_v - 1 - num_s,
    demand = f(row,7)
    if demand > 0.
        d(f(row,1)) = -demand/1000;
    end
end

d(s) = 10000;

cvx_begin
    variables q(num_e)
    minimize ( (1/3)*L*pow_p(q,3) )
    A*q <= d
cvx_end

q
h = A'\(diag(L)*(q.^2))
A'*h
diag(L)*(q.^2)
