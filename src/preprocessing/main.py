import sys
from optparse import OptionParser

import ABAQUSparser
import outputter


def convert(fin, fout):
    data = ABAQUSparser.Parser(fin).parse()
    outputter.Outputter(data, fout).print()


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
