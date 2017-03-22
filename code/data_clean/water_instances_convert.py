# Read file from Big.inp
f = open('Big.inp', 'r')

read_format = ''

node_id = 0
edge_id = 0

source = []
outlet = []
node_map = {}

node_str = ''
edge_str = ''

for line in f:
    if (line.isspace()):
        continue

    if ('[JUNCTIONS]' in line):
        read_format = 'node'
        next(f)
        continue

    if ('[RESERVOIRS]' in line):
        read_format = 'source'
        next(f)
        continue

    if ('[TANKS]' in line):
        read_format = 'tank'
        next(f)
        continue

    if ('[PIPES]' in line):
        read_format = 'pipe'
        next(f)
        continue

    if ('[PUMPS]' in line):
        read_format = ''
        continue

    if ('[VALVES]' in line):
        read_format = 'valve'
        next(f)
        continue

    if ('[TAGS]' in line):
        break

    if (read_format is 'node'):
        node_id += 1
        data = line.split()
        id = data[0]
        head = float(data[1])
        demand = float(data[2])
        node_map[id] = node_id
        if demand > 0:
            outlet.append(node_id)
        if 'NW' in id:
            demand = -1
            source.append(node_id)
        node_str += '%d %f %f\n' % (node_id, demand, head)

    if (read_format is 'source'):
        continue

    if (read_format is 'tank'):
        node_id += 1
        data = line.split()
        id = data[0]
        head = float(data[1])
        node_map[id] = node_id
        node_str += '%d %f %f\n' % (node_id, 0, head)

    if (read_format is 'pipe'):
        edge_id += 1
        data = line.split()
        node_1 = data[1]
        node_2 = data[2]
        length = float(data[3])
        diameter = float(data[4])
        roughness = float(data[5])
        edge_str += '%d %d %d %f %f %f\n' % (edge_id, node_map[node_1], node_map[node_2], length, diameter, roughness)

    if (read_format is 'valve'):
        edge_id += 1
        data = line.split()
        node_1 = data[1]
        node_2 = data[2]
        length = 0.1
        diameter = float(data[3])
        roughness = float(data[6])
        edge_str += '%d %d %d %f %f %f\n' % (edge_id, node_map[node_1], node_map[node_2], length, diameter, roughness)

# Write general info file
info = open('big.txt', 'w')
info.write('%d %d \n' % (node_id, edge_id))
for item in source:
    info.write('%d ' % (item))
info.write('\n')
for item in outlet:
    info.write('%d ' % (item))
info.close()

# Write node file
node_file = open('big_node.txt', 'w')
node_file.write(node_str)
node_file.close()

# Write edge file
edge_file = open('big_edge.txt', 'w')
edge_file.write(edge_str)
edge_file.close()
