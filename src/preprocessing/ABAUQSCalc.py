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
        self.materials = []


class Material():
    def __init__(self, args):
        self.args = args


class Calculator():
    def __init__(self):
        self.fout = open('data/Job-1.dat.2', 'w')

        self.fin1 = open('data/Job-1.dat.1', 'r')

        self.fin3 = open('data/Job-1.dat.3', 'r')

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
        self.elementGroups = []
        for eleGrpIndex in range(self.NEG):
            line = self.getLine3()
            eleType, nel, nmt = line.split()
            eleGrp = ElementGroup(int(eleType), int(nel), int(nmt))
            print(eleGrp)

            for i in range(eleGrp.nmt):
                _, *args = self.getLine3().split()
                material = Material(args)
                eleGrp.materials.append(material)

            for i in range(eleGrp.nel):
                _, *args = self.getLine3().split()
                ele = Element(
                    int(_),
                    tuple(int(item) for item in args[:-2]),
                    int(args[-1])
                )
                eleGrp.materials.append(ele)
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

    def calc(self):
        pass


def main():
    Calculator().run()


if __name__ == '__main__':
    main()
