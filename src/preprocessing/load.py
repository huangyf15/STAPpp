class Load():
    def __init__(self):
        self.forces = []


class Force():
    def __init__(self):
        self.node = None
        self.direction = None
        self.mag = None

    def format(self):
        return '%8d%4d%20lf' % (self.node.index, self.direction, self.mag)
