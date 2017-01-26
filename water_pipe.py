import numpy as np
# Change the file name
# 
f = open("by2.txt", "r")
data = []
g = 9.8

for line in f:
	data.append(line.replace('\r', '').replace('\n', '').replace(" ", "\t").split('\t'))
f.close()

# retrive the number of vertices and edges
num_v = int(data[0][0])
num_e = int(data[0][1])

# retrive the information of Source Node and End Node
source = data[1]
num_source = len(source)
end = data[2]
num_end = len(end)
end_demand = []

# Construct A and Lambda
A = [[0]*(num_e + num_end) for i in range(num_v + 1)]
L = [0]*(num_e + num_end)

for i in range(num_e):
	head = int(data[i + 3][1])
	tail = int(data[i + 3][2])
	length = float(data[i + 3][3])
	diameter = float(data[i + 3][4])/1000.0
	roughness = float(data[i + 3][5])/1000.0

	friction = (2 * np.log10(roughness/(3.71*diameter)))**2
	l = (8 * length * friction)/(np.pi**2*g*diameter**5)
	L[i] = l
	A[head - 1][i] = 1
	A[tail - 1][i] = -1

for i in range(num_end):
	head = int(end[i])
	tail = num_v + 1

	l = float(data[head + 2][8])/(float(data[head + 2][7])/1000.0)**2
	L[i + num_e] = l
	A[head - 1][i + num_e] = 1
	A[tail - 1][i + num_e] = -1

# Constrcut d
h = [0]*(num_v + 1)
for i in range(num_v):
	h[i] = float(data[i + 3][8])

d = [0.0]*(num_v + 1)
for i in source:
	d[int(i) - 1] = 10000.0

for num, i in enumerate(end):
	i = int(i) - 1
	j = num_v
	end_demand.append(float(data[i + 3][7])/1000.0)
	d[num_v] -= float(data[i + 3][7])/1000.0


#----------Stage 1: Use Maximum Flow Minimum Cost to caculate an initial solution----------#
from cvxpy import *
A = np.array(A)
d = np.array(d)
L = np.array(L)

q = Variable(num_e + num_end)
objective = Minimize(1/3.0 * sum(L * q**3))
constraints = [0 <= q, A*q <= d, q[8] >= 0.01, q[9] >= 0.01, q[10] >= 0.01]
prob = Problem(objective, constraints)
result = prob.solve(solver = "CVXOPT")

q_value = q.value


#---------Stage 2:    
for i in range(20):
	print "-----------------------Iteration %d-----------------------" %i
	# First Problem: fixed q_value, and add head constraint as variables
	q = Variable(num_e + num_end)
	hc = Variable(num_v + 1)
	penalty = sum(abs(A.transpose()*hc - np.diag(L)*np.square(q_value)))
	objective = Minimize(1/3.0 * sum(L * q**3) + penalty)
	constraints = [0 <= q, A*q <= d, q[8] >= 0.01, q[9] >= 0.01, q[10] >= 0.01, 0 <= hc, hc[9] == 0, hc[8] >= h[8], hc[4] >= h[4], hc[5] >= h[5]]
	prob = Problem(objective, constraints)
	result = prob.solve(solver = "CVXOPT")

	h_value = hc.value
	q_value = q.value
	print "Objective value is %.4f; Error values Ah-Lq^2: %.4f" % (objective.value - penalty.value, sum(abs(A.transpose()*h_value - np.diag(L)*np.square(q_value))).value)
	# Second Problem: fixed h_value, and only optimize the q
	q = Variable(num_e + num_end)
	tmp = A.transpose()*h_value
	tmp.A1[tmp.A1 < 0] = 0
	anchor = np.sqrt(tmp.A1/L)
	penalty = sum(abs(q - anchor))
	objective = Minimize(1/3.0 * sum(L * q**3) + 1000*penalty)
	constraints = [0 <= q, A*q <= d, q[8] >= 0.01, q[9] >= 0.01, q[10] >= 0.01]
	prob = Problem(objective, constraints)
	result = prob.solve(solver = "CVXOPT")

	q_value = q.value

	print "Objective value is %.4f; Error values Ah-Lq^2: %.4f" % (objective.value - penalty.value, sum(abs(A.transpose()*h_value - np.diag(L)*np.square(q_value))).value)
	print "Head value is as follow:"
	print h_value

