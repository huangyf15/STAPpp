class Load():
    def __init__(self):
        self.forces = []
    
    def __repr__(self):
        return f'<Load forces={str(self.forces)}>'


class Force():
    def __init__(self):
        self.node = None
        self.direction = None
        self.mag = None

    def format(self):
        return '%8d %4d %20lf' % (self.node.index, self.direction, self.mag)

    def __repr__(self):
        return f'<Force node={self.node.index}, direct={self.direction}, mag={self.mag}>'
