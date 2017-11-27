class Node():
    def __init__(self):
        self.index = None
        self.bounds = []
        self.pos = []

    def format(self):
        res = '%d' % self.index
        for item in self.bounds:
            res += '\t%d' % item
        for item in self.pos:
            res += '\t%lf' % item
        return res
    