class Node():
    self.index = None
    self.x = None
    self.y = None
    self.z = None

    def __init__(self, index, x, y, z):
        self.index = index
        self.x = x
        self.y = y
        self.z = z

__all__ = ['Node']
