import os
import re
import shutil

PATH = 'convergence rate data'

if os.path.exists(PATH):
    shutil.rmtree(PATH)

import genDat
genDat.main()

ns = []
for f in os.listdir(PATH):
    if re.match(r'\d+\.dat', f):
        ns.append(int(f[:-4]))

ns = sorted(ns)
for i in range(len(ns)):
    os.system('..' + os.sep + '..' + os.sep + 'build' + os.sep +
              'stap++ ' + '"%s%s%d.dat" > nul' % (PATH, os.sep, ns[i]))
    print('n = %d complete. < %d / %d >' % (ns[i], i + 1, len(ns)))

import calcErr
calcErr.main()
