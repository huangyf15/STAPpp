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


def verify(filename, plot):
    # print(filename)
    mark1 = False
    mark2 = False
    mark3 = False
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
                    index, dx, dy, dz, tx, ty, tz = line.split()[:7]
                    nodes[int(index)].dis(float(dx), float(dy), float(
                        dz), float(tx), float(ty), float(tz))
                except:
                    mark2 = False
                    continue
            elif mark3:
                try:
                    ele, n1, n2, _ = line.split()[:4]
                    elements[int(ele)] = Element([nodes[int(item)]
                                                  for item in (n1, n2)])
                except:
                    mark3 = False
                    continue
            else:
                if '   NUMBER  CONDITION  CODES                     COORDINATES' in line:
                    mark1 = True
                elif '  NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT' in line:
                    mark2 = True
                elif ' NUMBER-N      I        J       SET NUMBER' in line:
                    mark3 = True

    if plot:
        x = []
        y = []
        rx = []
        ry = []
        for _, node in nodes.items():
            x.append(node.x + TIMES * node.dx)
            y.append(node.y + TIMES * node.dy)
            rx.append(node.x + TIMES * (b / E * node.x))
            ry.append(node.y + TIMES * (-v * b / E * node.y))
        # plt.scatter(x, y)
        import matplotlib.pyplot as plt

        def f(n1, n2):
            x = [n1.x + TIMES * n1.dx, n2.x + TIMES * n2.dx]
            y = [n1.y + TIMES * n1.dy, n2.y + TIMES * n2.dy]
            plt.plot(x, y)

        lines = []
        for _, element in elements.items():
            n = element.nodes
            f(n[0], n[1])
            f(n[1], n[2])
            f(n[2], n[3])
            f(n[3], n[0])
        plt.scatter(rx, ry, color='r')
        plt.show()
        plt.savefig('patch-test-4Q.png')

    # 判断是否通过分片试验
    return all(testPass(node) for _, node in nodes.items())


def testPass(node):
    EI = 458533
    M = 100

    def dy():
        return M * (node.x - 10) * node.x / (2 * EI)

    def tz():
        return M * (node.x - 5) / EI

    res = equal(node.dx, 0) and equal(node.dy, dy()) and \
        equal(node.dz, 0) and equal(node.tx, 0) and equal(node.ty, 0) and \
        equal(node.tz, tz())
    if not res:
        print(f'failed at node {node}')
    return res


def equal(a, b):
    res = abs(a - b) <= max(abs(a) * 1e-5, abs(b) * 1e-5, 1e-15)
    if not res:
        print('failed equal: (a, b) = %s' % str((a, b)))
    return res


def main(plotFlag=False):
    return verify('patch.out', plotFlag)


if __name__ == '__main__':
    main()
