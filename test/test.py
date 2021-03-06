import os
import sys
import platform

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
    if os.system('msbuild stap++.vcxproj /m'):
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
    run('bar-6')
    run('test_truss_22')
    run('truss')

    os.chdir(DataDir)
    if os.system(PY + ' run-patch.py ' + STAP):
        quit(2)


def main():
    global PY
    if platform.system() == 'Windows':
        PY = 'python'
        winmake(False)
        test()
        winmake(True)
        test()
    elif platform.system() == 'Linux':
        PY = 'python3'
        unixmake(False)
        test()
        unixmake(True)
        test()
    else:
        print('unsupported platform')
    print('test passed.')


if __name__ == '__main__':
    main()
