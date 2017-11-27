class ElementGroup():
    def __init__(self):
        self.type = None
        self.elements = []
        self.materials = []


class Element():
    def __init__(self):
        self.index = None
        self.nodes = []
        self.material = None

    def format(self):
        res = '%d' % self.index
        for node in self.nodes:
            res += '\t%d' % node.index
        res += '\t%d' % self.material.index
        return res
