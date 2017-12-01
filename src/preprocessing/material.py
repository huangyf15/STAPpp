class Material():
    def __init__(self):
        self.index = None
        self.attributes = ()

    def __repr__(self):
        return f'<Material %d, attr: %s>' % (self.index, str(self.attributes))

    def format(self):
        res = '%d' % self.index
        for attr in self.attributes:
            res += '\t%lf' % attr
        return res
