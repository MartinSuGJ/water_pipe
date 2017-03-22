import numpy as np
# Change the file name
f = open("shamir2.txt", "r")
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


q_constraints = [0]*(num_e + num_end)
q_constraints[-num_end:] = end_demand

"""
Stage 1:
Use Maximum Flow Minimum Cost method to caculate an initial solution.
"""
from cvxpy import *
A = np.array(A)
d = np.array(d)
L = np.array(L)

q = Variable(num_e + num_end)
objective = Minimize(1/3.0 * sum(L * q**3))
constraints = [0 <= q, A*q <= d]
prob = Problem(objective, constraints)
result = prob.solve(solver = "CVXOPT")

q_value = q.value
h_value = np.linalg.lstsq(A.transpose(), np.diag(L)*np.square(q_value))[0]

"""
Stage 2:
Within this stage, we utilize an iterative procedure to find the approximate optimal solution
Sub-stage 2.1:
	Slightly adjust head pressure
Sub-stage 2.2:
	Slightly adjust flow
"""
  
for i in range(100):
	print "-----------------------Iteration %d-----------------------" %i
	# First Problem: fixed q_value, and add head constraint as variables
	hc = Variable(num_v + 1)
	objective = Minimize(15*norm2(max_elemwise(A.transpose()*hc - np.diag(L)*np.square(q_value), 0)) + 20*norm2(max_elemwise(np.diag(L)*np.square(q_value) - A.transpose()*hc, 0)) + 1*norm2(hc - h_value))
	#objective = Minimize(20*norm2(A.transpose()*hc - np.diag(L)*np.square(q_value)) + 1*norm2(hc - h_value))
	constraints = [0 <= hc, h <= hc, hc[-1] == 0]
	prob = Problem(objective, constraints)
	result = prob.solve(solver = "CVXOPT")

	h_value = hc.value
	
	# Second Problem: fixed h_value, and only optimize the q
	q = Variable(num_e + num_end)
	tmp = A.transpose()*h_value
	tmp.A1[tmp.A1 < 0] = 0
	anchor = np.sqrt(tmp.A1/L)
	penalty = sum(abs(q - anchor))
	objective = Minimize(1/3.0 * sum(L * q**3) + 1*norm2(q - anchor) + 1*norm2(q - q_value))
	constraints = [q_constraints <= q, A*q <= d]
	prob = Problem(objective, constraints)
	result = prob.solve(solver = "CVXOPT")

	q_value = q.value

	obj_adj = sum(L * q**3).value - q.value[-num_end:].transpose()*np.diag(L[-num_end:])*np.square(q.value[-num_end:]) + (q.T*(A.transpose()*h_value - np.diag(L)*np.square(q_value))).value
	head_diff = sum(abs(max_elemwise(np.diag(L)*np.square(q_value) - A.transpose()*h_value, 0))).value
	obj_diff = q.value.transpose()*max_elemwise(A.transpose()*h_value - np.diag(L)*np.square(q_value), 0).value
	print "Objective value is %.4f; Error values Ah-Lq^2: %.4f; Energy diff %.4f" % (obj_adj, head_diff, obj_diff)
	if head_diff < 0.00001:
		break


print "Head value is as follow:"
print h_value

