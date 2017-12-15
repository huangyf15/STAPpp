import re
import math
import numpy as np

from STAPpp import *
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
        self.parseLoad()

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
                elements[element.index] = element
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
        while True:
            try:
                line = self.getNextLine()
            except:
                break
            if '*Boundary' in line:
                self.parseBoundary()
            elif '*Dload' in line:
                self.parseDLoad()

    def parseBoundary(self):
        self.boundaries = []
        while True:
            line = self.getNextLine()
            try:
                nset, dis, rot = line.split(',')
                boundary = ABAQUS.Boundary()
                boundary.nsets = self.nsetsDict[nset]
                boundary.displacementIndex = int(dis)
                boundary.rotationIndex = int(rot)
                self.boundaries.append(boundary)
            except:
                self.goBack()
                break
        # print(self.boundaries)

    def parseDLoad(self):
        while True:
            line = self.getNextLine()
            try:
                _, __, mag, x, y, z = line.split(',')
                x, y, z = float(x), float(y), float(z)
                self.dload = ABAQUS.DLoad()
                self.dload.mag = float(mag)
                self.dload.direction = (x, y, z)
            except:
                self.goBack()
                break

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
        self.bar = ProgressBar(1000, printCount=False, printTime=True)
        self.bar.begin()
        for insName in self.instancesDict:
            instance = self.instancesDict[insName]
            part = instance.part
            for nodeLocalIndex in part.localNodesDict:
                localNode = part.localNodesDict[nodeLocalIndex]
                globalPos = calculatePos(instance, localNode)
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
                if localPointsCount / self.localPointsSum > \
                        self.bar.currentCount / self.bar.maxCount:
                    self.bar.grow()

        for bound in self.boundaries:
            for nset in bound.nsets:
                for nodeIndex in nset.nodeIndexs:
                    node = nset.instance.globalNodesDict[nodeIndex]
                    b = list(node.bounds)
                    b[bound.displacementIndex - 1] = 1
                    node.bounds = tuple(b)

        del self.bar

    def indexElements(self):
        elementTypeDict = {
            'S4R': 7,  # Shell
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
            if part.type == 'B31':  # for each beam instance, create a new material
                material = None
            else:
                material = self.stapppMaterialsDictByPartName.get(part.name)

            if not material:
                materialIndex = len(elementGroup.materials) + 1
                material = convertSection2Material(
                    part.type,
                    part.section, self.materialsDict, materialIndex)
                self.stapppMaterialsDictByPartName[part.name] = material
                elementGroup.materials.append(material)

            if part.type == 'B31':  # special material processing for beam
                fakeNode = ABAQUS.Node()
                fakeNode.pos = material.attributes[-3:]
                fakeNode.pos = calculatePos(instance, fakeNode)
                if instance.offset:
                    fakeNode.pos = tuple(fakeNode.pos[i] - instance.offset[i]
                                         for i in range(3))
                material.attributes = material.attributes[:-3] + fakeNode.pos
                # print(material.attributes)

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

    def calculateForce(self):
        load = Load()

        self.bar = ProgressBar(1000, printCount=False, printTime=True)

        fs = dict()  # fs marks every direction's force as {1:2.0, 2:3.0}
        for i in range(3):
            if self.dload.direction[i]:
                fs[i + 1] = self.dload.mag * self.dload.direction[i]

        # {index: {direction: force, }}
        self.globalForcesByGlobalNodeIndexDict = dict()

        count = 0  # for progress bar
        for insIndex in self.instancesDict:
            ins = self.instancesDict[insIndex]
            part = ins.part
            for elementIndex in part.localElementsDict:
                element = part.localElementsDict[elementIndex]

                # calculate a unit body force at each node
                # returns as {1:3.2, 2:2.4}
                nodeForces = self.calculateBodyForceAtElement(element, ins)

                for localNodeIndex in nodeForces:
                    globalNode = ins.globalNodesDict[localNodeIndex]
                    for direction in fs:
                        # get forces for global node
                        forces = self.globalForcesByGlobalNodeIndexDict.get(
                            globalNode.index)
                        if not forces:
                            forces = dict()
                            self.globalForcesByGlobalNodeIndexDict[globalNode.index] = forces

                        # get direction force from this node forces
                        force = forces.get(direction)
                        if not force:
                            force = Force()
                            force.direction = direction
                            force.node = globalNode
                            forces[direction] = force
                            load.forces.append(force)  # append new force only
                        force.mag += fs[direction] * nodeForces[localNodeIndex]

                count += 1
                if count / self.elementCount > self.bar.currentCount / self.bar.maxCount:
                    self.bar.grow()

        self.loads.append(load)
        del self.bar

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

        # 3. assembly force
        self.calculateForce()

        for index, eleGrp in sorted(self.eleGrpDict.items(), key=lambda x: x[0]):
            print(
                'Element Group %4d : %8d elements, %2d materials' % (
                    index, len(eleGrp.elements), len(eleGrp.materials)
                )
            )

        print('Total ABAQUS nodes : %d' % self.localPointsSum)
        print('Total STAPpp nodes : %d' % self.nodeCount)
        print('  boundary nodes   : %d' % len(self.linkedNodes))
        print('  inside nodes     : %d' %
              (self.nodeCount - len(self.linkedNodes)))
        print('Total forces       : %d' % len(self.loads[0].forces))

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
            'loads': self.loads
        }

    def calculateBodyForceAtElement(self, element, ins):
        part = ins.part
        _type = part.type
        section = part.section
        material = self.materialsDict[section.materialName]
        nodes = [part.localNodesDict[index]
                 for index in element.nodesIndexs]

        if _type == 'S4R':  # Plate
            thickness = section.args[0]
            # print(thickness)
            res = getGaussIntegrateFor4Q(
                element, nodes, thickness * material.density)
            return res

        elif _type == 'C3D8R':  # 8H
            return getGaussIntegrateFor8H(element, nodes, material.density)

        elif _type == 'B31':  # Beam
            length = abs(Vector(nodes[0].pos) - Vector(nodes[1].pos))

            args = section.args
            area = args[0] * args[1] - \
                (args[0] - args[2] - args[4]) * (args[1] - args[3] - args[5])

            grav = length * area * material.density
            return {index: grav / 2 for index in element.nodesIndexs}

        elif _type == 'T3D2':  # Bar
            length = abs(Vector(nodes[0].pos) - Vector(nodes[1].pos))
            grav = length * section.args[0] * material.density
            return {index: grav / 2 for index in element.nodesIndexs}

        else:
            raise Exception


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


def getGaussIntegrateFor4Q(element, nodes, c):
    # convert w from 3d to 2d
    if all(node.pos[2] == 0 for node in nodes):
        w = np.array([list(node.pos[:2]) for node in nodes])
    else:
        # for this particular case, there is no need to calc.
        raise
    # print(w)
    a = 1 / math.sqrt(3)
    b = (-a, a)
    s = np.array([0.0 for i in range(4)])
    for i in b:
        for j in b:
            s += (getGaussIntegrateFor4QAtPos((i, j), w) * 4)

    return {element.nodesIndexs[i]: s[i] * c for i in range(4)}


def getGaussIntegrateFor8H(element, nodes, c):
    w = np.array([list(node.pos) for node in nodes])
    # print(w)
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
    # print(GN)

    J = GN.dot(w)
    # print(J)
    Jdet = np.linalg.det(J)
    assert(Jdet > 0)

    res = N * Jdet
    # print(res)

    return res


def calculatePos(instance, localNode):
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


def convertSection2Material(elementType, section, materialsDict, index):
    material = materialsDict[section.materialName]
    res = Material()
    res.index = index
    # 'S4R': 7,  # Shell
    # 'C3D8R': 4,  # 8H
    # 'B31': 5,  # Beam
    # 'T3D2':
    if elementType == 'B31': # Beam
        res.attributes = (material.E, material.v, *section.args)
    elif elementType == 'C3D8R':
        assert len(section.args) == 0
        res.attributes = (material.E, material.v)
    elif elementType == 'S4R':
        assert len(section.args) == 2
        res.attributes = (material.E, material.v, section.args[0])
    elif elementType == 'T3D2': # Bar
        assert len(section.args) == 1
        res.attributes = (material.E, *section.args)
    else:
        raise Exception('unknown element type')
    return res


if __name__ == '__main__':
    Parser('data/Job-1.inp').parse()
    # ABAQUS nodes:
    # Job-1: 4163
    # Job-2: 37185
    # Job-3: 259567
    # Job-4: 1910327
