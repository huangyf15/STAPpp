import os
import json

from ProgressBar import ProgressBar
import numpy as np
import math
from ABAQUSparser import Vector


class Node():
    def __init__(self, pos):
        self.pos = pos


class Element():
    def __init__(self, index, nodesIndexs, materialIndex):
        self.index = index
        self.nodesIndexs = nodesIndexs
        self.materialIndex = materialIndex


class ElementGroup():
    def __init__(self, eleType, nel, nmt):
        self.eleType = eleType
        self.nel = nel
        self.nmt = nmt
        self.elements = []


class Calculator():
    def __init__(self, datname):
        self.fout = open(datname + '.2', 'w')

        self.fin1 = open(datname + '.1', 'r')

        self.fin3 = open(datname + '.3', 'r')

        self.jsonname = os.path.dirname(
            os.path.abspath(datname)) + os.sep + 'material.json'
        with open(self.jsonname) as f:
            self.materialInfo = json.load(f)

    def loadData(self):
        line = self.getLine1(3)
        non, nel, _, _ = line.split()
        self.NON = int(non)
        self.NEG = int(nel)
        print('NEN = %d, NEG = %d' % (self.NON, self.NEG))

        self.loadNodes()
        self.loadElements()

    def loadNodes(self):
        self.nodes = []
        self.getLine1()
        for i in range(self.NON):
            line = self.getLine1()
            indexNode, _, __, ___, x, y, z = line.split()
            node = Node((float(x), float(y), float(z)))
            self.nodes.append(node)

    def loadElements(self):
        self.bar = ProgressBar(1000, printTime=True, printCount=False)
        self.bar.nel = 0
        self.bar.count = 0
        self.elementGroups = []
        for eleGrpIndex in range(self.NEG):
            line = self.getLine3()
            eleType, nel, nmt = line.split()
            eleGrp = ElementGroup(int(eleType), int(nel), int(nmt))

            for i in range(eleGrp.nmt):
                self.getLine3()

            for i in range(eleGrp.nel):
                _, *args = self.getLine3().split()
                ele = Element(
                    int(_),
                    tuple(int(item) for item in args[:-1]),
                    int(args[-1])
                )
                eleGrp.elements.append(ele)
            self.bar.nel += eleGrp.nel
            self.elementGroups.append(eleGrp)

    def getLine1(self, n=1):
        for i in range(n):
            line = next(self.fin1)
        return line

    def getLine3(self):
        return next(self.fin3)

    def run(self):
        self.loadData()
        self.calc()
        self.output()
        os.remove(self.jsonname)

    def output(self):
        print('%8d %8d' % (1, len(self.forces)),
              file=self.fout)
        for nodeIndex in range(1, len(self.forces) + 1):
            force = self.forces[nodeIndex]
            print(
                '%8d %4d %20f' % (nodeIndex, 3, force),
                file=self.fout
            )

    def calc(self):
        self.forces = dict()

        self.direct = 3
        self.mag = -10

        for eleGrp in self.elementGroups:
            eleType = eleGrp.eleType
            for element in eleGrp.elements:
                # calculate a unit body force at each node
                # returns as {1:3.2, 2:2.4}
                nodeForces = self.calculateBodyForceAtElement(
                    element, eleGrp)

                for i in range(len(element.nodesIndexs)):
                    nodeIndex = element.nodesIndexs[i]

                    force = self.mag * nodeForces[nodeIndex]
                    if nodeIndex in self.forces:
                        self.forces[nodeIndex] += force
                    else:
                        self.forces[nodeIndex] = force
                self.bar.count += 1
                if self.bar.count / self.bar.nel > self.bar.currentCount / self.bar.maxCount:
                    self.bar.grow()

    def calculateBodyForceAtElement(self, element, eleGrp):
        eleType = eleGrp.eleType
        nodes = [self.nodes[index - 1] for index in element.nodesIndexs]
        args = self.materialInfo[str(eleType)][str(element.materialIndex)]
        if eleType == 7:  # Plate
            dens = args[0]
            thickness = args[1]
            res = getGaussIntegrateFor4Q(
                element, nodes, thickness * dens)
            return res

        elif eleType == 4:  # 8H
            dens = args[0]
            return getGaussIntegrateFor8H(element, nodes, dens)

        elif eleType == 5:  # Beam
            # dens, a, b, t1, t2, t3, t4
            dens = args[0]
            length = abs(Vector(nodes[0].pos) - Vector(nodes[1].pos))

            area = args[1] * args[2] - \
                (args[1] - args[3] - args[5]) * (args[2] - args[4] - args[6])

            grav = length * area * dens
            return {index: grav / 2 for index in element.nodesIndexs}

        elif eleType == 1:  # Bar
            dens = args[0]
            length = abs(Vector(nodes[0].pos) - Vector(nodes[1].pos))
            grav = length * args[1] * dens
            return {index: grav / 2 for index in element.nodesIndexs}

        else:
            raise Exception()


def getGaussIntegrateFor4Q(element, nodes, c):
    # convert w from 3d to 2d
    if all(node.pos[2] == 0 for node in nodes):
        w = np.array([list(node.pos[:2]) for node in nodes])
    else:
        # for this particular case, there is no need to calc.
        raise
    a = 1 / math.sqrt(3)
    b = (-a, a)
    s = np.array([0.0 for i in range(4)])
    for i in b:
        for j in b:
            s += getGaussIntegrateFor4QAtPos((i, j), w)

    return {element.nodesIndexs[i]: s[i] * c for i in range(4)}


def getGaussIntegrateFor8H(element, nodes, c):
    w = np.array([list(node.pos) for node in nodes])
    a = 1 / math.sqrt(3)
    b = (-a, a)
    s = np.array([0.0 for i in range(8)])
    for i in b:
        for j in b:
            for k in b:
                s += getGaussIntegrateFor8HAtPos((i, j, k), w)

    return {element.nodesIndexs[i]: s[i] * c for i in range(8)}


NGNs = dict()
NGN4Qs = dict()


def getGaussIntegrateFor4QAtPos(pos, w):
    global NGN4Qs
    if pos in NGN4Qs:
        N, GN = NGN4Qs[pos]
    else:
        a, b = pos
        N = np.array([
            (1 - a) * (1 - b),
            (1 + a) * (1 - b),
            (1 + a) * (1 + b),
            (1 - a) * (1 + b)]) / 4
        GN = np.array([
            [-1 + b, 1 - b, 1 + b, -1 - b],
            [-1 + a, -1 - a, 1 + a, 1 - a]]) / 4
        NGN4Qs[pos] = N, GN
        assert len(NGN4Qs) <= 4

    J = GN.dot(w)
    Jdet = np.linalg.det(J)
    assert(Jdet > 0)

    res = N * Jdet

    return res


def getGaussIntegrateFor8HAtPos(pos, w):
    global NGNs
    if pos in NGNs:
        N, GN = NGNs[pos]
    else:
        a, b, c = pos
        N = np.array(
            [(1 - a) * (1 - b) * (1 - c),
             (1 + a) * (1 - b) * (1 - c),
             (1 + a) * (1 + b) * (1 - c),
             (1 - a) * (1 + b) * (1 - c),
             (1 - a) * (1 - b) * (1 + c),
             (1 + a) * (1 - b) * (1 + c),
             (1 + a) * (1 + b) * (1 + c),
             (1 - a) * (1 + b) * (1 + c)]
        ) / 8
        GN = np.array(
            [[- (1 - b) * (1 - c),   (1 - b) * (1 - c),
                (1 + b) * (1 - c), - (1 + b) * (1 - c),
              - (1 - b) * (1 + c),   (1 - b) * (1 + c),
                (1 + b) * (1 + c), - (1 + b) * (1 + c)],
             [- (1 - a) * (1 - c), - (1 + a) * (1 - c),
                (1 + a) * (1 - c),   (1 - a) * (1 - c),
              - (1 - a) * (1 + c), - (1 + a) * (1 + c),
                (1 + a) * (1 + c),   (1 - a) * (1 + c)],
             [- (1 - a) * (1 - b), - (1 + a) * (1 - b),
              - (1 + a) * (1 + b), - (1 - a) * (1 + b),
                (1 - a) * (1 - b),   (1 + a) * (1 - b),
                (1 + a) * (1 + b),   (1 - a) * (1 + b)]]) / 8
        NGNs[pos] = N, GN
        assert len(NGNs) <= 8

    J = GN.dot(w)
    Jdet = np.linalg.det(J)
    assert(Jdet > 0)

    res = N * Jdet

    return res


def main():
    # Calculator('data/Job-1').run()
    a = 1 / math.sqrt(3.)
    res = np.array([0.0 for _ in range(4)])
    for i in (-a, a):
        for j in (-a, a):
            res += getGaussIntegrateFor4QAtPos(
                (i, j),
                np.array([
                    [-1, -1],
                    [2, -1],
                    [2, 1],
                    [-1, 1]
                ])
            )
            print(res)


if __name__ == '__main__':
    main()
