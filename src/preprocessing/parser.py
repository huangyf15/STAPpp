import re

from node import *
from load import *
from material import *
from element import *
import ABAQUS


class Parser():
    def __init__(self, fin):
        self.fin = fin
        self.heading = None
        self.nodes = []
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
        while True:
            line = self.getNextLine()
            if '*Instance' in line:
                instance = self.parseInstance()
                self.instancesDict[instance.name] = instance
            elif '*Nset' in line:
                nset = self.parseNset()
                if nset.name in self.nsetsDict:
                    self.nsetsDict[nset.name].append(nset)
                else:
                    self.nsetsDict[nset.name] = [nset, ]
            elif '*Elset' in line:
                self.parseElset()
            elif '*Surface' in line:
                self.parseSurface()
            elif '*Tie' in line:
                self.parseTie()
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
        self.materials = []
        while True:
            line = self.getNextLine()
            if '*Material' in line:
                self.goBack()
                self.materials.append(self.parseMaterial())
            if '*Step' in line:
                self.goBack()
                break
        print(self.materials)

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
        surfaceName = re.match(
            r'\*Surface, type=NODE, name=([^,]+), internal', self.getLine()
        ).groups()[0]
        nset = re.match(r'[^,]+', self.getNextLine())[0]
        print(surfaceName, nset)
        print('Surface parsed.')

    def parseTie(self):
        line = self.getLine()
        ratatioin = 'rotation' not in line
        adjust = 'adjust=yes' in line
        # Tie name is useless

        line = self.getNextLine()
        surf1, surf2 = line.split(',')
        surf1 = surf1.replace(' ', '')
        surf2 = surf2.replace(' ', '')
        print('tie parsed, sf1 = %s, sf2 = %s' % (surf1, surf2))

    def data(self):
        return {
            'heading': self.heading,
            'nodes': self.nodes,
            'elementGroups': self.eleGrp,
            'loads': self.loads
        }

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


if __name__ == '__main__':
    print(Parser('data/Job-1.inp').parse())
