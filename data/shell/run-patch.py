import os
import re
import verify
import sys

res = 0 == os.system('..' + os.sep + '..' + os.sep + 'build' + os.sep + 'stap++ patch.dat > nul')
res = res and verify.main()
if res:
    print('patch test passed.')
    sys.exit(0)
else:
    print('patch test failed.')
    sys.exit(1)