import re
import math

from node import *
from load import *
from material import *
from element import *
import ABAQUS


class Parser():
    def __init__(self, fin):
        self.fin = fin
        self.heading = None
        self.nodesDict = dict()
        self.eleGrp = []
        self.loads = []

    def parse(self):
        '''
        return parsed data.
        {
            'heading': headingStr,
            'nodes': [Node],
            'elementGroups': [ElementGroups],
            'loads': [Load]
        }
        '''

        with open(self.fin, 'r') as fp:
            self.lines = fp.readlines()
        self.index = -1

        self.parseHeading()
        self.parseParts()
        self.parseAssembly()
        self.parseMaterials()
        # self.parseLoad()

        del self.lines

        self.analyse()

        return self.data()

    def parseHeading(self):
        self.gotoKeyWord('*Heading')
        self.heading = self.getNextLine()[3:]

    def parseParts(self):
        self.partsDict = dict()
        while True:
            line = self.getNextLine()
            if '*Part' in line:
                print('parsing %s' % line)
                self.goBack()
                part = self.parsePart()
                self.partsDict[part.name] = part
            elif '*Assembly' in line:
                self.goBack()
                break
        print('%d part parsed.' % len(self.partsDict))

    def parsePart(self):
        part = ABAQUS.Part()

        while True:
            line = self.getNextLine()
            if '*End Part' in line:
                break
            elif '*Part' in line:
                name = re.match(r'\*Part, name=([\w\-]+)', line).groups()[0]
                part.name = name
            elif '*Node' in line:
                nodes = self.parseNode()
                part.localNodesDict = nodes
            elif '*Element' in line:
                _type, elements = self.parseElement()
                part.type = _type
                part.localElementsDict = elements
            elif '*Solid Section' in line:
                material = re.findall(r'material=(\w+)', self.getLine())[0]
                part.materialName = material

        print('part parsed')
        # print(part)
        return part

    def parseNode(self):
        'returns a dict of node {localIndex: ABAQUS.Node}'
        nodes = dict()
        while True:
            try:
                line = self.getNextLine()
                index, x, y, z = line.split(',')
                index = int(index)
                x, y, z = float(x), float(y), float(z)
                node = ABAQUS.Node()
                node.pos = (x, y, z)
                nodes[index] = node
            except:
                self.goBack()
                break
        return nodes

    def parseElement(self):
        ' returns [type, {localIndex: ABAQUS.ELEMENT}] '
        line = self.getLine()
        _type = re.match(r'\s*\*Element, type=(\w+)', line).groups()[0]
        elements = dict()
        while True:
            try:
                line = self.getNextLine()
                index, *nodeIndex = line.split(',')
                element = ABAQUS.Element()
                element.index = int(index)
                element.nodesIndexs = tuple(int(item) for item in nodeIndex)
                elements[index] = element
            except:
                self.goBack()
                break
        return [_type, elements]

    def parseAssembly(self):
        '从 *Assembly 到 *End Assembly'
        self.instancesDict = dict()
        self.nsetsDict = dict()  # {name: [nsets]}
        self.surfacesDict = dict()
        self.ties = []
        while True:
            line = self.getNextLine()
            if '*Instance' in line:
                instance = self.parseInstance()
                assert instance.name not in self.instancesDict
                self.instancesDict[instance.name] = instance
            elif '*Nset' in line:
                nset = self.parseNset()
                if nset.name in self.nsetsDict:
                    self.nsetsDict[nset.name].append(nset)
                else:
                    self.nsetsDict[nset.name] = [nset, ]
            # elif '*Elset' in line:
            #     self.parseElset()
            elif '*Surface' in line:
                surface = self.parseSurface()
                assert surface.name not in self.surfacesDict
                self.surfacesDict[surface.name] = surface
            elif '*Tie' in line:
                self.ties.append(self.parseTie())
            elif '*End Assembly' in line:
                break
        return

    def parseInstance(self):
        instance = ABAQUS.Instance()
        line = self.getLine()
        name, partName = re.match(
            r'\*Instance, name=([^,]+), part=([\w\-]+)', line).groups()

        print('parsing instance %s' % name)

        instance.name = name
        instance.part = self.partsDict[partName]

        if '*End Instance' in self.getNextLine():
            pass
        else:
            # see http://wufengyun.com:888/v6.14/books/key/default.htm?startat=ch09abk19.html#ksuperprop-rot-instance
            instance.offset = tuple(float(item)
                                    for item in self.getLine().split(','))
            if '*End Instance' in self.getNextLine():
                pass
            else:
                instance.rotation = tuple(float(item)
                                          for item in self.getLine().split(','))
                assert '*End Instance' in self.getNextLine()

        print('Instance parsed.')
        print(instance)
        return instance

    def parseMaterials(self):
        self.materialsDict = dict()
        while True:
            line = self.getNextLine()
            if '*Material' in line:
                self.goBack()
                material = self.parseMaterial()
                self.materialsDict[material.name] = material
            if '*Step' in line:
                self.goBack()
                break
        print(self.materialsDict)

    def parseMaterial(self):
        material = ABAQUS.Material()
        material.name = re.match(
            r'\*Material, name=(\w+)', self.getLine()).groups()[0]

        self.getNextLine()
        material.density = float(self.getNextLine()[:-1])

        self.getNextLine()
        material.E, material.v = (float(item)
                                  for item in self.getNextLine().split(','))

        return material

    def parseLoad(self):
        pass

    def parseNset(self):
        line = self.getLine()
        # print(line)
        nset = ABAQUS.Nset()
        res = re.match(r'\*Nset, nset=([\_\-\w]+), instance=([\w+\-]+)', line)
        assert res
        nset.name, instanceName = res.groups()
        nset.instance = self.instancesDict[instanceName]

        indexs = []
        while True:
            line = self.getNextLine()
            try:
                indexs += [int(item) for item in line.split(',') if item]
            except:
                if '*Nset' in line or '*Elset' in line or '*Surface' in line:
                    self.goBack()
                    break
                else:
                    raise
        nset.nodeIndexs = tuple(indexs)
        print('Nset parsed. %s' % nset)
        return nset

    def parseElset(self):
        line = self.getLine()
        # print(line)
        elset = ABAQUS.Elset()
        res = re.match(
            r'\*Elset, elset=([\_\-\w]+), instance=([\w+\-]+)', line)
        assert res
        elset.name, elset.instance = res.groups()

        while True:
            line = self.getNextLine()
            try:
                indexs = [int(item) for item in line.split(',') if item]
                print(indexs)
            except:
                if '*Nset' in line or '*Elset' in line or '*Surface' in line:
                    self.goBack()
                    break
                else:
                    print(line)
                    raise
        print('Elset parsed.')
        return elset

    def parseSurface(self):
        surface = ABAQUS.Surface()
        surface.name = re.match(
            r'\*Surface, type=NODE, name=([^,]+), internal', self.getLine()
        ).groups()[0]
        surface.nsetName = re.match(r'[^,]+', self.getNextLine())[0]
        print('Surface parsed. %s' % surface)
        return surface

    def parseTie(self):
        tie = ABAQUS.Tie()
        line = self.getLine()
        # Tie name is useless actually
        tie.name = re.match(r'\*Tie, name=([\w\-]+)', line).groups()[0]
        tie.rotation = 'rotation' not in line
        tie.adjust = 'adjust=yes' in line

        line = self.getNextLine()
        surf1, surf2 = line.split(',')
        tie.surfaceName1 = surf1.replace(' ', '')
        tie.surfaceName2 = surf2.replace(' ', '')

        print('tie parsed %s' % tie)
        return tie

    def data(self):
        pass
        # return {
        #     'heading': self.heading,
        #     'nodes': self.nodes,
        #     'elementGroups': self.eleGrp,
        #     'loads': self.loads
        # }

    def gotoKeyWord(self, keyword):
        while keyword not in self.getNextLine():
            pass

    def goBack(self):
        self.index -= 1

    def getNextLine(self):
        self.index += 1
        return self.getLine()

    def getLine(self):
        res = self.lines[self.index][:-1]
        if '**' in res:
            self.index += 1
            return self.getLine()
        else:
            return res

    def analyse(self):
        '''
        vars(self):
            heading, nodes, eleGrp, loads;
            partsDict, instancesDict, nsetsDict,
            surfaceDict, ties, materialDict
        '''

        # 0. bond ties-surf-nset to instance
        for tie in self.ties:
            surf1 = self.surfacesDict[tie.surfaceName1]
            surf2 = self.surfacesDict[tie.surfaceName2]

            nsets1 = self.nsetsDict[surf1.nsetName]
            nsets2 = self.nsetsDict[surf2.nsetName]

            # print('len(nsets1) = %d' % (len(nsets1)))
            # print('len(nsets2) = %d' % (len(nsets2)))

            for nset in nsets1:  # left surface node sets
                # print('len(nset) = %d' % len(nset.nodeIndexs))
                instance = nset.instance
                instance.nsets.append(nset)

            for nset in nsets2:
                instance = nset.instance
                instance.nsets.append(nset)

        # 1. indexing
        nodeCount = 0
        elementCount = 0

        for insName in self.instancesDict:
            instance = self.instancesDict[insName]
            part = instance.part
            for nodeLocalIndex in part.localNodesDict:
                localNode = part.localNodesDict[nodeLocalIndex]
                # first calculate pos
                globalPos = calculatePos(instance, nodeLocalIndex)
                #  see if node has shown in other instances, that is,
                # show up in ties
                #  if so, skip new index
                if isLocalNodeIndexInInstanceNsets(nodeLocalIndex, instance):
                    pass
                    # TODO
                else:
                    nodeCount += 1
                    globalNode = Node()
                    globalNode.index = nodeCount
                    # TODO: calculate pos
                    instance.globalNodesDict[nodeLocalIndex] = globalNode

            for elementLocalIndex in part.localElementsDict:
                pass
                # calculate gravity force
        print(nodeCount)


class Vector():
    def __init__(self):
        self.__init__((0., 0., 0.))

    def __init__(self, tup):
        self.x, self.y, self.z = tup

    def __sub__(self, vec):
        return Vector((self.x - vec.x, self.y - vec.y, self.z - vec.z))

    def __add__(self, vec):
        return Vector((self.x + vec.x, self.y + vec.y, self.z + vec.z))

    def __mul__(self, a):
        return Vector((a * self.x, a * self.y, a * self.z))

    def __abs__(self):
        return math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def normalize(self):
        return self * (1 / abs(self))

    def __repr__(self):
        return '<Vector %s>' % str((self.x, self.y, self.z))

    def dot(self, vec):
        return self.x * vec.x + self.y * vec.y + self.z * vec.z

    def cross(self, vec):
        return Vector((
            self.y * vec.z - self.z * vec.y,
            self.z * vec.x - self.x * vec.z,
            self.x * vec.y - self.y * vec.x
        ))


def calculatePos(instance, index):
    localNode = instance.part.localNodesDict[index]
    pos = tuple(localNode.pos)
    if instance.offset:
        pos = tuple(pos[i] + instance.offset[i] for i in range(3))
    if instance.rotation:
        pa = Vector(instance.rotation[:3])
        pb = Vector(instance.rotation[3:6])
        phi = instance.rotation[6] * math.pi / 180
        p1 = Vector(pos)

        p1a = p1 - pa
        i = (pb - pa).normalize()
        pc1 = p1a - i * p1a.dot(i)
        l = abs(pc1)
        if l != 0:
            j = pc1.normalize()
            n = i.cross(j)

            p2 = pa + i * p1a.dot(i) + j * l * math.cos(phi) + \
                n * l * math.sin(phi)

            pos = (p2.x, p2.y, p2.z)
        else:
            # 1 on line ab, return pos directly
            pass
    # print(pos)
    return pos


def isLocalNodeIndexInInstanceNsets(localNodeIndex, instance):
    for nset in instance.nsets:
        if localNodeIndex in nset.nodeIndexs:
            return True
    return False


if __name__ == '__main__':
    print(Parser('data/Job-1.inp').parse())
