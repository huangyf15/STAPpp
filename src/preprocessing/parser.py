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
        # self.parseMaterial()
        # self.parseLoad()

        return self.data()

    def parseHeading(self):
        self.gotoKeyWord('*Heading')
        self.heading = self.getNextLine()[3:]

    def parseParts(self):
        self.parts = []
        while True:
            line = self.getNextLine()
            if '*Part' in line:
                print('parsing %s' % line)
                self.goBack()
                part = self.parsePart()
                self.parts.append(part)
            elif '*Assembly' in line:
                self.goBack()
                break
        print('%d part parsed.' % len(self.parts))

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
                part.nodes = nodes
            elif '*Element' in line:
                _type, elements = self.parseElement(part.nodes)
                part.type = _type
                part.elements = elements
            # elif '*Nset' in line:
            #     nset = self.parseNset()
            #     part.nset = nset
            # elif '*Elset' in line:
            #     elset = self.parseElset()
            #     part.elset = elset
            elif '*Solid Section' in line:
                material = re.findall(r'material=(\w+)', self.getLine())[0]
                part.material = material

        print('part parsed')
        return part
        # quit()

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
        name, partName = re.match(
            r'\*Instance, name=([^,]+), part=([\w\-]+)', line).groups()

        print('parsing instance %s' % name)

        part = [part for part in self.parts if part.name == partName][0]

        if '*End Instance' in self.getNextLine():
            dx, dy, dz = 0, 0, 0
        else:
            # see http://wufengyun.com:888/v6.14/books/key/default.htm?startat=ch09abk19.html#ksuperprop-rot-instance
            dx, dy, dz = (float(item) for item in self.getLine().split(','))
            if '*End Instance' in self.getNextLine():
                pass
            else:
                rt_a_x, rt_a_y, rt_a_z, rt_b_x, rt_b_y, rt_b_z, rt_phi = (
                    float(item) for item in self.getLine().split(','))
                assert '*End Instance' in self.getNextLine()

        print('Instance parsed.')

    def parseMaterial(self):
        pass

    def parseLoad(self):
        pass

    def parseNset(self):
        line = self.getLine()
        nset = ABAQUS.Nset()
        res = re.match(r'\s*\*Nset, nset=([\_\w]+), instance=([\w+\-]+)', line)
        if res:
            nset.name, nset.instance = res.groups()
            # print(nset.name, nset.instance)
        else:
            nset.name, = re.match(r'\s*\*Nset, nset=([\_\w]+)', line).groups()
            # print(nset.name)

        line = self.getNextLine()
        args = line.split(',')
        if len(args) == 2:
            nset.indexs = (int(args[0]), int(args[1]))
        elif len(args) == 3 and int(args[2]) == 1:
            nset.indexs = range(int(args[0], int(args[1]) + 1, 1))
        else:
            raise RuntimeError('undefined nset len args')

        return nset

    def parseElset(self):
        line = self.getLine()
        elset = ABAQUS.Elset()
        res = re.match(
            r'\s*\*Elset, elset=([\_\w]+), instance=([\w+\-]+)', line)
        if res:
            elset.name, elset.instance = res.groups()
            # print(nset.name, nset.instance)
        else:
            elset.name, = re.match(
                r'\s*\*Elset, elset=([\_\w]+)', line).groups()
            # print(nset.name)

        line = self.getNextLine().replace(',', '')
        args = line.split()
        if len(args) == 1:
            elset.indexs = (int(args[0]),)
        elif len(args) == 2:
            elset.indexs = (int(args[0]), int(args[1]))
        elif len(args) == 3 and int(args[2]) == 1:
            elset.indexs = range(int(args[0], int(args[1]) + 1, 1))
        else:
            raise RuntimeError('undefined elset len args')

        return elset

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
        if '**' in res:
            self.index += 1
            return self.getLine()
        else:
            return res


if __name__ == '__main__':
    print(Parser('data/Job-1.inp').parse())
