import sys
import re
import math
import os

import matplotlib.pyplot as plt
import numpy as np

b = 10
E = 1.0e6
v = 0.3

PATH = 'convergence rate data'


def calc(n):
    x = []
    y = []
    rx = []
    ry = []

    with open(PATH + os.sep + '%d.out' % n, 'r') as f:
        err = 0.0
        for line in f.readlines():
            if len(line) > 80 and 'e-' in line:
                items = line.split()
                assert(len(items) == 12)
                GPX, GPY = (float(item) for item in items[2:4])
                GPUX, GPUY = (float(item) for item in items[5:7])
                weight = float(items[11])

                realX = b * (1 - v * v) * GPX * (12 - GPX * GPX) / (6 * E)
                realY = 0

                x.append(GPX + GPUX * 10000)
                rx.append(GPX + realX * 10000)
                y.append(GPY + GPUY * 10000)
                ry.append(GPY + realY * 10000)

                err += weight * ((realX - GPUX) ** 2 + (realY - GPUY) ** 2)

    plt.scatter(x, y)
    plt.scatter(rx, ry, color='r')
    # plt.show()

    err = math.sqrt(err)
    print('{%d, %g},' % (n, err))
    return err


def main():
    ns = []
    errs = []
    for filename in os.listdir(PATH):
        if re.match(r'\d+\.dat', filename):
            ns.append(int(filename[:-4]))
    for n in sorted(ns):
        err = calc(n)
        errs.append(err)

    hs = [1 / n for n in sorted(ns)]
    loghs = [math.log(_h) for _h in hs]
    logerrs = [math.log(err) for err in errs]
    a, b = np.polyfit(loghs, logerrs, 1)

    func = 'log(error) = %.2f log(h) %s %.2f' % (
        a, '+' if b > 0 else '-', abs(b))

    plt.close()
    plt.plot(loghs, logerrs)
    plt.scatter(loghs, logerrs)
    plt.xlabel('log(h)')
    plt.ylabel('log(error)')
    plt.title(func)

    # plt.show()
    plt.savefig('3T convergence rate.png')


if __name__ == '__main__':
    main()
