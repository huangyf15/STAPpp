import re
import math

from node import *
from load import *
from material import *
from element import *
import ABAQUS

from ProgressBar import ProgressBar


class Parser():
    def __init__(self, fin):
        self.fin = fin
        self.heading = None
        self.globalNodesDict = dict()
        self.eleGrpDict = dict()
        self.stapppMaterialsDictByPartName = dict()
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
        self.bar = ProgressBar(1000, printCount=False, printTime=True)
        self.index = -1

        self.parseHeading()
        self.parseParts()
        self.parseAssembly()
        self.parseMaterials()
        # self.parseLoad()

        self.bar.update(self.bar.maxCount)
        del self.bar
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
                # print('parsing %s' % line)
                self.goBack()
                part = self.parsePart()
                self.partsDict[part.name] = part
            elif '*Assembly' in line:
                self.goBack()
                break
        # print('%d part parsed.' % len(self.partsDict))

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
            elif re.match(r'\*\w+ Section', line):
                section = self.parseSection()
                part.section = section

        return part

    def parseSection(self):
        line = self.getLine()
        section = ABAQUS.Section()
        section.materialName = re.findall(r'material=(\w+)', self.getLine())[0]

        while True:
            line = self.getNextLine()
            try:
                section.args += tuple(float(item)
                                      for item in line.split(',') if item)
            except Exception as ex:
                self.goBack()
                break

        return section

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

        # print('parsing instance %s' % name)

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

        # print('Instance parsed. %s' % instance)
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
        # print(self.materialsDict)

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
        # print('Nset parsed. %s' % nset)
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
        # print('Elset parsed.')
        return elset

    def parseSurface(self):
        surface = ABAQUS.Surface()
        surface.name = re.match(
            r'\*Surface, type=NODE, name=([^,]+), internal', self.getLine()
        ).groups()[0]
        surface.nsetName = re.match(r'[^,]+', self.getNextLine())[0]
        # print('Surface parsed. %s' % surface)
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

        # print('tie parsed %s' % tie)
        return tie

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
        # self.bar.update(self.index)
        if (self.index / len(self.lines)) > (self.bar.currentCount / self.bar.maxCount):
            self.bar.grow()
        if '**' in res:
            self.index += 1
            return self.getLine()
        else:
            return res

    def bondTieToInstance(self):
        for tie in self.ties:
            surf1 = self.surfacesDict[tie.surfaceName1]
            surf2 = self.surfacesDict[tie.surfaceName2]

            nsets1 = self.nsetsDict[surf1.nsetName]
            nsets2 = self.nsetsDict[surf2.nsetName]

            for nset in nsets1:  # left surface node sets
                # print('len(nset) = %d' % len(nset.nodeIndexs))
                instance = nset.instance
                instance.linkedNsets.append(nset)

            for nset in nsets2:
                instance = nset.instance
                instance.linkedNsets.append(nset)

    def indexNodes(self):
        self.linkedNodes = []
        self.nodeCount = 0

        self.localPointsSum = sum(
            len(self.instancesDict[insindex].part.localNodesDict)
            for insindex in self.instancesDict)
        localPointsCount = 0
        bar = ProgressBar(1000, printCount=False, printTime=True)
        bar.begin()
        for insName in self.instancesDict:
            instance = self.instancesDict[insName]
            part = instance.part
            for nodeLocalIndex in part.localNodesDict:
                localNode = part.localNodesDict[nodeLocalIndex]
                globalPos = calculatePos(instance, nodeLocalIndex)
                if isLocalNodeIndexInInstanceNsets(nodeLocalIndex, instance):
                    globalNode = self.matchGlobalNode(globalPos)
                    if globalNode:  # already counted
                        instance.globalNodesDict[nodeLocalIndex] = globalNode
                    else:  # new boundary node
                        self.nodeCount += 1
                        globalNode = Node(globalPos)
                        globalNode.index = self.nodeCount
                        globalNode.bounds = (0, 0, 0)
                        instance.globalNodesDict[nodeLocalIndex] = globalNode
                        self.globalNodesDict[self.nodeCount] = globalNode
                        self.linkedNodes.append(globalNode)

                else:  # inside nodes
                    self.nodeCount += 1
                    globalNode = Node(globalPos)
                    globalNode.index = self.nodeCount
                    globalNode.bounds = (0, 0, 0)
                    instance.globalNodesDict[nodeLocalIndex] = globalNode
                    self.globalNodesDict[self.nodeCount] = globalNode

                localPointsCount += 1
                if localPointsCount / self.localPointsSum > bar.currentCount / bar.maxCount:
                    bar.grow()

    def indexElements(self):
        elementTypeDict = {
            'S4R': 6,  # Plate
            'C3D8R': 4,  # 8H
            'B31': 5,  # Beam
            'T3D2': 1  # Bar
        }

        self.elementCount = 0
        for insIndex in self.instancesDict:
            instance = self.instancesDict[insIndex]
            part = instance.part
            elementType = elementTypeDict[part.type]
            elementGroup = self.getElementGroup(elementType)

            # get material for this part
            material = self.stapppMaterialsDictByPartName.get(part.name)
            if not material:
                materialIndex = len(elementGroup.materials) + 1
                material = convertSection2Material(
                    part.section, self.materialsDict, materialIndex)
                self.stapppMaterialsDictByPartName[part.name] = material
                elementGroup.materials.append(material)

            for localElementIndex in part.localElementsDict:
                localElement = part.localElementsDict[localElementIndex]
                globalElement = Element()
                self.elementCount += 1
                globalElement.index = len(elementGroup.elements) + 1
                globalElement.material = material  # same material for a instance
                globalElement.nodes = [instance.globalNodesDict[index]
                                       for index in localElement.nodesIndexs]
                elementGroup.elements.append(globalElement)

    def getElementGroup(self, elementType):
        'elementType as in number'
        if elementType not in self.eleGrpDict:
            eleGrp = ElementGroup()
            eleGrp.type = elementType
            self.eleGrpDict[elementType] = eleGrp

        return self.eleGrpDict[elementType]

    def analyse(self):
        '''
        vars(self):
            heading, nodes, eleGrp, loads;
            partsDict, instancesDict, nsetsDict,
            surfaceDict, ties, materialDict
        '''

        # 0. bond ties-surf-nset to instance
        self.bondTieToInstance()

        # 1. indexing nodes
        self.indexNodes()

        # 2. index elements
        self.indexElements()

        # 3. assembly
        for index, eleGrp in self.eleGrpDict.items():
            print(
                'Element Group %d: %d elements, %d materials' % (
                    index, len(eleGrp.elements), len(eleGrp.materials)
                )
            )

        print('totle ABAQUS nodes: %d' % self.localPointsSum)
        print('totle STAPpp nodes: %d' % self.nodeCount)
        print('  linked nodes: %d' % len(self.linkedNodes))
        print('  unlinked nodes: %d' %
              (self.nodeCount - len(self.linkedNodes)))

        # with open('mma.dat', 'w') as f:
        #     for index in range(self.nodeCount):
        #         node = self.globalNodesDict[index + 1]
        #         f.write('%f\t%f\t%f\n' % node.pos)

    def matchGlobalNode(self, gPos):
        for node in self.linkedNodes:
            pos = node.pos
            err = sum((gPos[i] - pos[i])**2 for i in range(3))
            if err < 1e-10:
                return node
        return None

    def data(self):
        return {
            'heading': self.heading,
            'nodes': tuple(self.globalNodesDict.values()),
            'elementGroups': tuple(self.eleGrpDict.values()),
            'loads': []
        }


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
    for nset in instance.linkedNsets:
        if localNodeIndex in nset.nodeIndexs:
            return nset
    return False


def convertSection2Material(section, materialsDict, index):
    material = materialsDict[section.materialName]
    res = Material()
    res.index = index
    res.attributes += (material.E, material.v) + section.args
    return res


if __name__ == '__main__':
    Parser('data/Job-1.inp').parse()
    # Job-1: 4163 nodes
    # Job-2: 37185
    # Job-3: 232720??
    # Job-4: 1910327
