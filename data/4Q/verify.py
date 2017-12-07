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
                    index, _, __, ___, x, y, z = line.split()
                    nodes[int(index)] = Node(float(x), float(y), float(z))
                except:
                    mark1 = False
                    continue
            elif mark2:
                try:
                    index, dx, dy, dz = line.split()
                    nodes[int(index)].dis(float(dx), float(dy), float(dz))
                except:
                    mark2 = False
                    continue
            elif mark3:
                try:
                    ele, n1, n2, n3, n4, _ = line.split()
                    elements[int(ele)] = Element([nodes[int(item)]
                                                  for item in (n1, n2, n3, n4)])
                except:
                    mark3 = False
                    continue
            else:
                if '   NUMBER  CONDITION  CODES                     COORDINATES' in line:
                    mark1 = True
                elif '  NODE           X-DISPLACEMENT    Y-DISPLACEMENT    Z-DISPLACEMENT' in line:
                    mark2 = True
                elif ' NUMBER-N      I        J        K        L      SET NUMBER' in line:
                    mark3 = True

    b = 10
    E = 1.0e6
    v = 0.3

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
        # plt.show()
        # plt.savefig('patch-test-4Q.png')

    # 判断是否通过分片试验
    return all(
        equal(node.dx, b / E * node.x) and equal(node.dy, -v * b / E * node.y)
        for _, node in nodes.items())


def equal(a, b):
    res = abs(a - b) <= max(abs(a) * 1e-6, abs(b) * 1e-6, 1e-15)
    if not res:
        print('failed (a, b) = %s' % str((a, b)))
    return res


def main(plotFlag=False):
    return verify('patch.out', plotFlag)


if __name__ == '__main__':
    main()
