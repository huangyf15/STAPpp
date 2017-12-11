import os
import re

TIMES = 10000


class Node():
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def dis(self, dx, dy, dz, tx, ty, tz):
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.tx = tx
        self.ty = ty
        self.tz = tz

    def __repr__(self):
        return '<Node (%f, %f, %f), dis=(%f, %f, %f), rot=(%f, %f, %f)>' % (
            self.x, self.y, self.z,
            self.dx, self.dy, self.dz,
            self.tx, self.ty, self.tz
        )


class Element():
    def __init__(self, nodes):
        self.nodes = nodes


def verify(filename):
    # print(filename)
    mark1 = False
    mark2 = False
    nodes = dict()

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
                    index, dx, dy, dz, tx, ty, tz = line.split()[:7]
                    nodes[int(index)].dis(float(dx), float(dy), float(
                        dz), float(tx), float(ty), float(tz))
                except:
                    mark2 = False
                    continue
            else:
                if all(item in line for item in ('NUMBER', 'CONDITION', 'CODES', 'COORDINATES')):
                    mark1 = True
                elif all(item in line for item in ('NODE', 'X-DISPLACEMENT', 'Y-DISPLACEMENT', 'Z-DISPLACEMENT')):
                    mark2 = True

    # 判断是否通过分片试验
    return all(testPass(node) for _, node in nodes.items())


def testPass(node):

    def dx():
        return 0

    def dy():
        return 0

    def dz():
        return node.x**2 - 0.2*node.y**2

    def tx():
        return -0.4 * node.y

    def ty():
        return -2 * node.x
    
    def tz():
        return 0

    res = equal(node.dx, dx()) and equal(node.dy, dy()) and \
        equal(node.dz, dz()) and equal(node.tx, tx()) and equal(node.ty, ty()) and \
        equal(node.tz, 0)
    if not res:
        print('failed at node %s' % node)
    return res


def equal(a, b):
    res = abs(a - b) <= max(abs(a) * 1e-5, abs(b) * 1e-5, 1e-13)
    if not res:
        print('failed equal: (a, b) = %s' % str((a, b)))
    return res


def main():
    return verify('patch.out')


if __name__ == '__main__':
    main()
