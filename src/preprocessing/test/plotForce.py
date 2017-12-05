import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math


class Ploter():
    def __init__(self, filename):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        f = open(filename)
        self.lines = f.readlines()

        self.index = -1
        ps = dict()

        line = self.getLine(3)
        pointCount = int(line.split()[0])
        self.getLine()

        for i in range(pointCount):
            ind, _, __, ___, x, y, z = self.getLine().split()
            ps[int(ind)] = float(x), float(y), float(z)

        self.getLine(2)
        x = []
        y = []
        z = []
        s = []
        for i in range(pointCount):
            ind, _, mag = self.getLine().split()
            ind = int(ind)
            mag = abs(float(mag))
            p = ps[ind]
            x.append(p[0])
            y.append(p[1])
            z.append(p[2])
            s.append(mag)

        m = min(s)
        s = [item / m for item in s]
        # s = [math.sqrt(s) for s in s]
        m = max(s)
        s = [50 * item / m for item in s]
        # s = [math.log(item) for item in s]
        # s = [1000 * s[i] / max(s) for i in range(len(s))]
        ax.scatter(x, y, z, s=s)

        plt.show()

    def getLine(self, x=1):
        self.index += x
        return self.lines[self.index]


if __name__ == '__main__':
    Ploter('../data/Job-2.dat')
