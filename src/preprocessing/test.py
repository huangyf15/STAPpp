from element import *
from load import *
from node import *
from material import *
from Output import Outputter

nodes = []
for i in range(4):
    node = Node()
    node.index = i + 1
    node.bounds = [1, 1, 1]
    node.pos = [i, i + 1, i + 2]
    nodes.append(node)

material = Material()
material.index = 1
material.attributes = {'E': 1e6, 'v': 0.3}

element = Element()
element.index = 1
element.nodes = nodes
element.material = material

eleGrp = ElementGroup()
eleGrp.type = 1
eleGrp.elements = [element]
eleGrp.materials = [material]

load = Load()
force = Force()
force.direction = 1
force.mag = 1e5
force.node = nodes[0]
load.forces = [force]

data = {
    'heading':'Hello, world!',
    'elementGroups':[eleGrp],
    'nodes':nodes,
    'loads':[load]
}

Outputter(data, 'out.dat').print()
