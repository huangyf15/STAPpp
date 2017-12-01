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

        self.nsets = []
        self.globalNodesDict = dict()

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
        self.instance = None
        self.nodeIndexs = None

    def __repr__(self):
        return f'<Nset name={self.name}, instance={self.instance.name}, indexs={self.nodeIndexs}>'


class Elset(Nset):
    pass


class Surface():
    def __init__(self):
        self.name = None
        self.nsetName = None

    def __repr__(self):
        return f'<Surface name={self.name}, nset={self.nsetName}>'


class Tie():
    def __init__(self):
        self.name = None
        self.adjust = None
        self.rotation = None
        self.surfaceName1 = None
        self.surfaceName2 = None

    def __repr__(self):
        return f'<Tie name={self.name}>, adjust={self.adjust}, ' + \
               f'rot={self.rotation}, surf1={self.surfaceName1}, surf2={self.surfaceName2}>'


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
