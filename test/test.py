import os
import sys

ProjectDir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DataDir = ProjectDir + os.sep + 'data'
STAP = None


def winmake(mkl=False):
    os.chdir(ProjectDir)
    if not os.path.exists('vsbuild'):
        os.mkdir('vsbuild')
    os.chdir(ProjectDir + os.sep + 'vsbuild')
    if os.system('cmake ../src -G "Visual Studio 15 2017 Win64" -DUSE_MKL=%s' %
                     ('ON' if mkl else 'OFF')):
        quit(1)
    if os.system('msbuild stap++.sln'):
        print('build failed.')
        quit(1)
    global STAP
    STAP = ProjectDir + os.sep + 'vsbuild' + \
        os.sep + 'Debug' + os.sep + 'stap++.exe'
    os.chdir(ProjectDir)


def unixmake(mkl=False):
    os.chdir(ProjectDir)
    if not os.path.exists('build'):
        os.mkdir('build')
    os.chdir(ProjectDir + os.sep + 'build')
    if os.system('cmake ../src -DUSE_MKL=%s' % ('ON'if mkl else 'OFF')):
        quit(1)
    if os.system('make'):
        print('build failed.')
        quit(1)
    global STAP
    STAP = ProjectDir + os.sep + 'build' + os.sep + 'stap++'
    os.chdir(ProjectDir)


def test():
    def run(name):
        if os.system(STAP + ' ' + DataDir + os.sep + name + ' > nul'):
            print('test failed for file ' + name)
            quit(2)
    # run('bar-6')
    # run('test_truss_22')
    run('truss')

    def run(name):
        os.chdir(DataDir + os.sep + name)
        if os.system(PY + ' run-patch.py ' + STAP):
            print('patch test failed for element type ' + name)
            quit(2)
    run('3T')
    run('4Q')
    run('8H')
    run('Beam')
    run('plate')
    run('shell')
    run('TimoEBMOD')
    run('TimoSRINT')


def main():
    global PY
    if len(sys.argv) == 1:
        print('Usage: py test.py <win/unix>')
        quit(1)
    if sys.argv[1] == 'win':
        PY = 'py'
        winmake(False)
        test()
        winmake(True)
        test()
    elif sys.argv[1] == 'unix':
        PY = 'python3'
        unixmake(False)
        test()
        unixmake(True)
        test()
    print('test passed.')


if __name__ == '__main__':
    main()
