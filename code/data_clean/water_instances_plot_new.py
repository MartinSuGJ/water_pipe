# Read file from Small.inp
f = open('Big.inp', 'r')
node_file = open('big.nodes', 'w')
edge_file = open('big.tets', 'w')

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

    if ('[VERTICES]' in line):
        break

    if ('[COORDINATES]' in line):
        next(f)
        read_format = 'coordinates'
        continue

    if (read_format is 'coordinates'):
        data = line.split()
        id = data[0]
        x = float(data[1])
        y = float(data[2])
        node_map[id] = node_id
        node_str = '%f %f %f\n' % (x, y, 0.0)
        node_file.write(node_str)
        if id == 'XXXX000816':
            outlet.append(node_id)
        if 'W0' in id:
            source.append(node_id)
        node_id += 1
f.close()

f = open('Big.inp', 'r')

for line in f:
    if (line.isspace()):
            continue

    if ('[PIPES]' in line):
        read_format = 'pipe'
        next(f)
        continue

    if ('[PUMPS]' in line):
        read_format = 'pump';
        next(f)
        continue

    if ('[VALVES]' in line):
        read_format = 'valve'
        next(f)
        continue

    if ('[TAGS]' in line):
        break

    if (read_format is 'pipe'):
        edge_id += 1
        data = line.split()
        node_1 = data[1]
        node_2 = data[2]
        edge_str = '%d %d\n' % (node_map[node_1], node_map[node_2])
        edge_file.write(edge_str)

    if (read_format is 'valve'):
        edge_id += 1
        data = line.split()
        # print data
        node_1 = data[1]
        node_2 = data[2]
        edge_str = '%d %d\n' % (node_map[node_1], node_map[node_2])
        edge_file.write(edge_str)

    if (read_format is 'pump'):
        edge_id += 1
        data = line.split()
        node_1 = data[1]
        node_2 = data[2]
        edge_str = '%d %d\n' % (node_map[node_1], node_map[node_2])
        edge_file.write(edge_str)



node_file.close()
edge_file.close()

print source
print outlet
