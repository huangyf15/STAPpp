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

    def __repr__(self):
        return f'<Element {self.index}, material={self.material}, nodes={self.nodes}>'

    def format(self):
        res = '%5d' % self.index
        for node in self.nodes:
            res += '%8d' % node.index
        res += '%5d' % self.material.index
        return res


class Node():
    def __init__(self):
        self.index = None
        self.bounds = None
        self.pos = None

    def __init__(self, pos):
        self.index = None
        self.bounds = None
        self.pos = pos

    def format(self):
        res = '%8d' % self.index
        for item in self.bounds:
            res += ' %3d' % item
        for item in self.pos:
            res += ' %12lf' % item
        return res

    def __repr__(self):
        return '<Node %d %s %s>' % (self.index, str(self.pos), str(self.bounds))


class Load():
    def __init__(self):
        self.forces = []

    def __repr__(self):
        return f'<Load forces={str(self.forces)}>'


class Force():
    def __init__(self):
        self.node = None
        self.direction = None
        self.mag = 0.0

    def format(self):
        return '%8d %4d %20lf' % (self.node.index, self.direction, self.mag)

    def __repr__(self):
        return f'<Force node={self.node.index}, direct={self.direction}, mag={self.mag}>'

class Material():
    def __init__(self):
        self.index = None
        self.attributes = ()

    def __repr__(self):
        return f'<Material %d, attr: %s>' % (self.index, str(self.attributes))

    def format(self):
        res = '%5d' % self.index
        for attr in self.attributes:
            res += ' %12lf' % attr
        return res
