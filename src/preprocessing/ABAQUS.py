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