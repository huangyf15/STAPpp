import re

from node import *
from load import *
from material import *
from element import *

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
        self.parseAssembly()
        self.parseMaterial()
        self.parseLoad()

        return self.data()

    def parseHeading(self):
        self.gotoKeyWord('*Heading')
        self.heading = self.getNextLine()[3:]

    def parseAssembly(self):
        self.gotoKeyWord('*Assembly')

        self.parseInstance()

        self.gotoKeyWord('*End Assembly')

    def parseInstance(self):
        self.gotoKeyWord('*Instance')

        # 读点信息
        self.gotoKeyWord('*Node')
        while True:
            line = self.getNextLine()
            try:
                line = line.replace(', ', '')
                index, x, y, z = line.split()
                node = Node()
                node.index = int(index)
                node.pos = (float(x), float(y), float(z))
                self.nodes.append(node)
            except:
                self.goBack()
                break

        # 读单元信息
        self.gotoKeyWord('*Element')
        line = self.getLine()
        elementType = re.findall(r'type=(\w+)', line)[0]
        element = 
        

    def parseMaterial(self):
        pass

    def parseLoad(self):
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
        return self.lines[self.index][:-1]

if __name__ == '__main__':
    print(Parser('Job-1.inp').parse())
