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
        self.parseMaterial()
        self.parseLoad()

        return self.data()

    def parseHeading(self):
        self.gotoKeyWord('*Heading')
        self.heading = self.getNextLine()[3:]

    def parseParts(self):
        self.gotoKeyWord('*Part')
        self.goBack()

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
                part.nodes = nodes
            elif '*Element' in line:
                _type, elements = self.parseElement(part.nodes)
                part.type = _type
                part.elements = elements
            elif 'Nset' in line:
                self.parseNset()
                pass
        quit()

    def parseNode(self):
        'returns a dict of node'
        nodes = dict()
        while True:
            try:
                line = self.getNextLine()
                index, x, y, z = line.split(',')
                index = int(index)
                x, y, z = float(x), float(y), float(z)
                node = Node()
                node.index = int(index)
                node.pos = x, y, z
                nodes[index] = node
            except:
                self.goBack()
                break
        return nodes

    def parseElement(self, nodes):
        ' returns [type, [elements]] '
        line = self.getLine()
        _type = re.match(r'\s*\*Element, type=(\w+)', line).groups()[0]
        elements = []
        while True:
            try:
                line = self.getNextLine()
                index, *nodeIndex = line.split(',')
                element = Element()
                element.index = int(index)
                element.nodes = [nodes[float(item)] for item in nodeIndex]
                elements.append(element)
            except:
                self.goBack()
                break
        return [_type, elements]

    def parseAssembly(self):
        '从 *Assembly 到 *End Assembly'
        self.gotoKeyWord('*Assembly')

        while True:
            line = self.getNextLine()
            if '*Instance' in line:
                self.parseInstance()
            elif '*Nset' in line:
                self.parseNset()
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
        line = self.getLine()
        res = re.match(r'\*Instance, name=([^,]+), part=([\w\-]+)', line)
        name, part = res.groups()

        try:  # 读取偏移
            line = self.getNextLine()
            dx, dy, dz = line.split(',')
            dx, dy, dz = float(dx), float(dy), float(dz)
        except:
            self.goBack()
            dx, dy, dz = 0, 0, 0

        while True:
            line = self.getNextLine()
            if '*Node' in line:
                self.parseNode()
            elif '*Element' in line:
                self.parseElement()
            elif '*Nset' in line:
                self.parseNset()
            elif '*Elset' in line:
                self.parseElset()
            elif '*Solid' in line:
                self.parseSolid()
            elif '*End Instance' in line:
                break

    def parseMaterial(self):
        pass

    def parseLoad(self):
        pass

    def parseNset(self):
        pass

    def parseElset(self):
        pass

    def parseSurface(self):
        pass

    def parseTie(self):
        pass

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
        if res[:2] == '**':
            self.index += 1
            return self.getLine()
        else:
            return res


if __name__ == '__main__':
    print(Parser('data/Job-1.inp').parse())
