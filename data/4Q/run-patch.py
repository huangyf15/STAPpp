import os
import re
import verify
import sys

os.system('..' + os.sep + '..' + os.sep + 'build' + os.sep + 'stap++ patch.dat > nul')
res = verify.main()
if res:
    print('patch test passed.')
    sys.exit(0)
else:
    print('patch test failed.')
    sys.exit(1)
