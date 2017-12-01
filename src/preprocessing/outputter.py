class Outputter():
    def __init__(self, data, fout):
        self.data = data
        self.fout = fout

    def print(self):
        self.f = open(self.fout, 'w')
        self.printHeadingLine()
        self.printControlLines()
        self.printNodalDataLines()
        self.printLoadDataLines()
        self.printElementDataLines()
        self.f.close()

    def printHeadingLine(self):
        print(self.data['heading'], file=self.f)
        print(file=self.f)

    def printControlLines(self):
        print('%d\t%d\t%d\t1' % (
            len(self.data['nodes']),
            len(self.data['elementGroups']),
            len(self.data['loads'])
        ), file=self.f)
        print(file=self.f)

    def printNodalDataLines(self):
        for node in self.data['nodes']:
            print(node.format(), file=self.f)
        print(file=self.f)

    def printLoadDataLines(self):
        count = 1
        for load in self.data['loads']:
            print('%d\t%d' % (count, len(load.forces)), file=self.f)
            count += 1
            for force in load.forces:
                print(force.format(), file=self.f)
        print(file=self.f)

    def printElementDataLines(self):
        for elementGroup in self.data['elementGroups']:
            print('%d\t%d\t%d' % (
                elementGroup.type,
                len(elementGroup.elements),
                len(elementGroup.materials)
            ), file=self.f)
            for material in elementGroup.materials:
                print(material.format(), file=self.f)
            for element in elementGroup.elements:
                print(element.format(), file=self.f)
        print(file=self.f)

