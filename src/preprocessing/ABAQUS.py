class Part():
    def __init__(self):
        self.name = None
        self.nodes = []
        self.Element = []
        self.Nset = None


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
