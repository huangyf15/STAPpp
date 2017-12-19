import os
import verify
import sys


def run(name):
    root = os.path.dirname(os.path.abspath(__file__)) + os.sep + name + os.sep
    res = os.system(sys.argv[1] + ' %spatch.dat > nul' % root)
    if res:
        print('running for patch.dat for element %s failed!' % name)
        sys.exit(1)

    res = verify.verify(name)
    if res:
        print('patch test for element %s passed.' % name)
    else:
        print('patch test failed for element %s!' % name)
        sys.exit(1)


def main():
    run('3T')
    print('patch test passed for all.')


if __name__ == '__main__':
    main()
