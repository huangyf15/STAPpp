class Part():
    def __init__(self):
        self.name = None
        self.type = None
        self.localNodesDict = []
        self.localElementsDict = None
        self.materialName = None

    def __repr__(self):
        return f'<Part name={self.name}, ' + \
            f'type={self.type}, nodes={self.localNodesDict}, ' + \
            f'elements={self.localElementsDict}, material={self.materialName}'


class Node():
    def __init__(self):
        self.pos = None

    def __repr__(self):
        return f'<Node {self.pos}>'


class Instance():
    def __init__(self):
        self.name = None
        self.part = None

        self.offset = None
        self.rotation = None

    def __repr__(self):
        return f'<Instance name={self.name}, part={self.part.name}, ' + \
               f'offset={self.offset}, rotation={self.rotation}'


class Element():
    def __init__(self):
        '注意这里没有 type 信息'
        # self.index = None
        self.nodesIndexs = None

    def __repr__(self):
        return f'<Element nodes={self.nodesIndexs}>'


class Nset():
    def __init__(self):
        self.name = None
        self.index = None


class Elset(Nset):
    pass


class Material():
    def __init__(self):
        self.name = None
        self.density = None
        self.E = None
        self.v = None

    def __repr__(self):
        return '<Material %s, rho=%f, E=%f, v=%f>' % (
            self.name, self.density, self.E, self.v
        )
