import sys
import os
import random

f = None
b = 10
E = 1.0e6
v = 0.3

PATH = 'convergence rate data'


def printHeadingLine(n):
    print('4Q convergence rate test file, n = %d' % n, file=f)


def printControlLine(n):
    print('%d  1  1  1' % (
        (n + 1) * (2 * n + 1),  # 节点总数
    ), file=f)


def printNodalPointDataLines(n):
    for i in range(2 * n + 1):
        for j in range(n + 1):
            print('%-5d %3d  %d  1   %f  %f  0' % (
                (n + 1) * i + j + 1,  # 编号
                i == 0,    # x 方向是否自由，即是否在最左边一行
                # j == 0,    # y 方向是否自由，即 y=0 为 1
                1,  # 总是固定
                1 / n * i,  # x 坐标
                1 / n * j  # y 坐标
            ), file=f)


def printLoadDataLines(n):
    print('1  %d' % (
        (2 * n) * (n + 1)  # 除去左边的节点数
    ), file=f)

    for i in range(1, 2 * n + 1):
        # loadLeft = (b * i / (2 * n * n) - b / (8 * n * n)) / n
        # loadRight = (b * i / (2 * n * n) + b / (8 * n * n)) / n
        loadLeft = b * (-2 + 3 * i) / (6 * (n**3))
        loadRight = b * (2 + 3 * i) / (6 * (n**3))
        # loadLeft = b / (2 * n * n)  # 常体力
        # loadRight = b / (2 * n * n)

        if i == 2 * n:
            load = loadLeft
        else:
            load = loadLeft + loadRight

        for j in range(n + 1):
            if j == 0 or j == n:
                rload = load / 2
            else:
                rload = load
            nodeIndex = (n + 1) * i + j + 1
            print('%-5d  1   %f' % (nodeIndex, rload), file=f)


def printElementDataLines(n):
    print('2   %d   1' % (
        n * n * 2,  # 单元总数
    ), file=f)

    print('1   %f   %f' % (E, v), file=f)  # 杨氏模量和泊松比

    for i in range(2 * n):
        for j in range(n):
            print('%-5d %5d  %5d  %5d  %5d    1' % (
                n * i + j + 1,  # 单元编号
                (n + 1) * i + j + 1,  # 左下
                (n + 1) * i + j + 2,  # 左上
                (n + 1) * (i + 1) + j + 2,  # 右上
                (n + 1) * (i + 1) + j + 1  # 右下
            ), file=f)


def run(n):
    if not os.path.exists(PATH):
        os.mkdir(PATH)
    global f
    f = open(PATH + os.sep + '%d.dat' % n, 'w')
    printHeadingLine(n)
    print(file=f)
    printControlLine(n)
    print(file=f)
    printNodalPointDataLines(n)
    print(file=f)
    printLoadDataLines(n)
    print(file=f)
    printElementDataLines(n)
    print(file=f)
    f.close()

def main():
    ns = []
    ns += range(1, 33)
    ns += [64, 128]
    ns = sorted(set(ns))
    for i in ns:
        run(i)
    # print(ns)

if __name__ == '__main__':
    main()
