import sys
from optparse import OptionParser
import os

import ABAQUSparser
import outputter
import ABAQUSCalc


def convert(fin, fout):
    data = ABAQUSparser.Parser(fin).parse()
    outputter.Outputter(data, fout).print()
    del data
    ABAQUSCalc.Calculator(fout).run()

    with open(fout, 'w') as f:
        for i in range(3):
            finname = fout + '.' + str(i+1)
            with open(finname) as fin_:
                for line in fin_:
                    f.write(line)
            os.remove(finname)



def main():
    parser = OptionParser('Usage: main.py [options] inputFile')
    parser.add_option(
        '-o', '--output',
        dest='dest', default=None,
        help='output dat file to dest'
    )

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        exit()

    fin = args[0]
    if fin[-4:] != '.inp':
        fin += '.inp'

    fout = options.dest
    if not fout:
        fout = fin[:-4] + '.dat'
    elif fout[-4:] != '.dat':
        fout += '.dat'

    convert(fin, fout)


if __name__ == '__main__':
    main()
