class Material():
    def __init__(self):
        self.index = None
        self.attributes = dict()

    def format(self):
        res = '%d' % self.index
        for key in self.attributes:
            res += '\t%lf' % self.attributes[key]
        return res
