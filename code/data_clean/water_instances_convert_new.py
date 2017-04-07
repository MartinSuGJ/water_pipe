# Read file from Big.inp
f = open('Small.inp', 'r')
node_file = open('small_node.txt', 'w')
edge_file = open('small_edge.txt', 'w')

read_format = ''

node_id = 0
edge_id = 0

source = []
outlet = []
node_map = {}

node_str = ''
edge_str = ''
edge_map = {}

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
        read_format = 'pump'
        next(f)
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
        node_type = 0
        if demand > 0:
            outlet.append(node_id)
            node_type = 1
        node_str = '%d %f %f %d\n' % (node_id, demand, head, node_type)
        node_file.write(node_str)

    if (read_format is 'source'):
        node_id += 1
        data = line.split()
        id = data[0]
        head = float(data[1])
        demand = -1
        source.append(node_id)
        node_map[id] = node_id
        node_type = 2
        node_str = '%d %f %f %d\n' % (node_id, demand, head, node_type)
        node_file.write(node_str)

    if (read_format is 'tank'):
        node_id += 1
        data = line.split()
        id = data[0]
        head = float(data[1])
        node_map[id] = node_id
        node_type = 3
        node_str = '%d %f %f %d\n' % (node_id, 0, head, node_type)
        node_file.write(node_str)

    if (read_format is 'pipe'):
        edge_id += 1
        data = line.split()
        node_1 = data[1]
        node_2 = data[2]
        edge_map[edge_id] = node_1 + "," + node_2
        length = float(data[3])
        diameter = float(data[4])
        roughness = float(data[5])
        edge_type = 0
        edge_str = '%d %d %d %f %f %f %d\n' % (edge_id, node_map[node_1], node_map[node_2], length, diameter, roughness, edge_type)
        edge_file.write(edge_str)

    if (read_format is 'pump'):
        edge_id += 1
        data = line.split()
        node_1 = data[1]
        node_2 = data[2]
        edge_map[edge_id] = node_1 + "," + node_2
        length = 0.1
        diameter = 250.0
        roughness = 1.50
        edge_type = 1
        edge_str = '%d %d %d %f %f %f %d\n' % (edge_id, node_map[node_1], node_map[node_2], length, diameter, roughness, edge_type)
        edge_file.write(edge_str)

    if (read_format is 'valve'):
        edge_id += 1
        data = line.split()
        node_1 = data[1]
        node_2 = data[2]
        edge_map[edge_id] = node_1 + "," + node_2
        length = 0.1
        diameter = float(data[3])
        roughness = float(data[6])
        edge_type = 2
        edge_str = '%d %d %d %f %f %f %d\n' % (edge_id, node_map[node_1], node_map[node_2], length, diameter, roughness, edge_type)
        edge_file.write(edge_str)

# Write general info file
info = open('small.txt', 'w')
info.write('%d %d \n' % (node_id, edge_id))
for item in source:
    info.write('%d ' % (item))
info.write('\n')
for item in outlet:
    info.write('%d ' % (item))
info.close()


node_file.close()
edge_file.close()
