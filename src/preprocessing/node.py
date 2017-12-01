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
