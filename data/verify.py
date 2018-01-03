import os
import re
import json
import copy


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


def getRoot():
    return os.path.dirname(os.path.abspath(__file__))


def getFileName(name):
    return getRoot() + os.sep + name + os.sep + 'patch.out'


def verify(name):
    mark1 = False
    mark2 = False
    nodes = dict()

    with open(getFileName(name), 'r') as fp:
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
    return all(testPass(name, node) for _, node in nodes.items())


def testPass(name, node):
    with open(getRoot() + os.sep + 'accurate.json') as f:
        js = json.load(f)[name]

    if 'args' in js:
        globals().update(js['args'])

    for item in ('dx', 'dy', 'dz', 'tx', 'ty', 'tz'):
        if item in js:
            globals()[item + 's'] = js[item]
        else:
            globals()[item + 's'] = '0'

    globals()['node'] = node
    globals()['x'] = node.x
    globals()['y'] = node.y
    globals()['z'] = node.z

    def dx(): return eval(dxs)

    def dy(): return eval(dys)

    def dz(): return eval(dzs)

    def tx(): return eval(txs)

    def ty(): return eval(tys)

    def tz(): return eval(tzs)

    res = equal(node.dx, dx()) and equal(node.dy, dy()) and \
        equal(node.dz, dz()) and equal(node.tx, tx()) and \
        equal(node.ty, ty()) and equal(node.tz, tz())
    if not res:
        print('failed at node %s' % node)
    return res


def equal(a, b):
    delta = abs(a - b)
    threshold = max(abs(a) * 1e-4, abs(b) * 1e-4, 1e-11)
    res = delta < threshold
    if not res:
        print('failed equal: (a, b) = %s' % str((a, b)))
        print('delta = %g, threshold = %g' % (delta, threshold))
    return res


if __name__ == '__main__':
    print(verify('3T'))
