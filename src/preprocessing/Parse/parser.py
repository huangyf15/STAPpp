class Parser():
    def __init__(self, fin):
        self.fin = fin

    def parse(self):
        '''
        return parsed data.
        {
            'heading': headingStr,
            'nodes': [Node],
            'elementGroups': [ElementGroups],
            'loads': [Load]
        }
        '''

        with open(fin, 'r') as fp:
            self.lines = fp.readlines()
        

