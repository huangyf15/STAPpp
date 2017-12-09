import os
import re

TIMES = 10000


class Node():
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def dis(self, dx, dy, dz):
        self.dx = dx
        self.dy = dy
        self.dz = dz

    def __repr__(self):
        return '<Node (%f, %f, %f), dis=(%f, %f, %f)>' % (
            self.x, self.y, self.z,
            self.dx, self.dy, self.dz
        )


class Element():
    def __init__(self, nodes):
        self.nodes = nodes


def verify(filename):
    # print(filename)
    mark1 = False
    mark2 = False
    nodes = dict()
    elements = dict()

    with open(filename, 'r') as fp:
        for line in fp:
            if mark1:
                try:
                    index, _, __, ___, x, y, z = line.split()[:7]
                    nodes[int(index)] = Node(float(x), float(y), float(z))
                except:
                    mark1 = False
                    continue
            elif mark2:
                try:
                    index, dx, dy, dz = line.split()[:4]
                    nodes[int(index)].dis(float(dx), float(dy), float(dz))
                except:
                    mark2 = False
                    continue
            else:
                if '   NUMBER  CONDITION  CODES                     COORDINATES' in line:
                    mark1 = True
                elif '  NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT' in line:
                    mark2 = True

    # 判断是否通过分片试验
    return all(testPass(node) for _, node in nodes.items())


def testPass(node):
    E = 1e6
    v = 0.2
    P = 0.25e3
    l_x = 5
    l_y = 5
    l_z = 5

    def dx():
        return -v * 4 * P * node.x / (E * l_z * l_x)

    def dy():
        return -v * 4 * P * node.y / (E * l_z * l_y)

    def dz():
        return 4 * P * node.z / (E * l_z * l_z)

    res = equal(node.dx, dx()) and equal(
        node.dy, dy()) and equal(node.dz, dz())
    if not res:
        print('failed at node %s' % node)
    return res


def equal(a, b):
    res = abs(a - b) <= max(abs(a) * 1e-6, abs(b) * 1e-6, 1e-15)
    if not res:
        print('failed equal: (a, b) = %s' % str((a, b)))
    return res


def main():
    return verify('patch.out')


if __name__ == '__main__':
    main()
