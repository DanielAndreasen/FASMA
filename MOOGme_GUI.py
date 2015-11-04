#!/usr/bin/env python
# -*- coding: utf8 -*-

from MOOGme import moogme
import argparse
from gooey import Gooey, GooeyParser


def ew(args):
    '''Create the StarMe.cfg file from the GUI'''
    fout = ''
    linelist = args.linelist.rpartition('/')[-1]
    if args.spectralType:
        if not args.temperature and not args.surfacegravity:
            fout += '%s spt:%s,' % (linelist, args.spectralType)
    else:
        if not args.temperature:
            print 'Warning: Temperature was set to 5777K'
            args.temperature = 5777
        if not args.surfacegravity:
            print 'Warning: Surface gravity was set to 4.44'
            args.surfacegravity = 4.44
        if not args.FeH:
            args.FeH = 0.0
        if not args.microturbulence:
            print 'Warning: Microturbulence was set to 1.00'
            args.microturbulence = 1.00
        fout += '%s %s %s %s %s ' % (linelist, args.temperature, args.surfacegravity, args.FeH, args.microturbulence)

    fout += 'model:%s,iterations:%s,weights:%s,rwslope:%s,epslope:%s,abdiff:%s,MOOGv:%s' % (args.model, args.Iterations, args.weights, args.RWslope, args.EPslope, args.Fedifference, args.MOOGv)
    if args.Fixteff:
        fout += ',teff'
    if args.Fixgravity:
        fout += ',logg'
    if args.FixFeH:
        fout += ',feh'
    if args.Fixmicroturbulence:
        fout += ',vt'
    with open('StarMe.cfg', 'w') as f:
        f.writelines(fout)
    moogme()


def synth():
    pass


def abund():
    pass


@Gooey(program_name='MOOG Made Easy - for deriving stellar parameters', default_size=(700, 1000))
def main():
    '''Take care of all the argparse stuff.
    :returns: the args
    '''
    parser = GooeyParser(description='Set parameters')

    subparsers = parser.add_subparsers()

    # Common to all
    parent_parser = GooeyParser(add_help=False)
    parent_parser.add_argument('--temperature',     help='Input initial temperature',      default=5777,  type=int)
    parent_parser.add_argument('--surfacegravity',  help='Input initial gravity',          default=4.44,  type=float, metavar='logg')
    parent_parser.add_argument('--FeH',             help='Input initial metallicity',      default='0.00',type=float, metavar='[Fe/H]')
    parent_parser.add_argument('--microturbulence', help='Input initial microturbulence',  default=1.0,   type=float)
    parent_parser.add_argument('--MOOGv',           help='Version of MOOG', default='2013', choices=['2013', '2014'], type=str, metavar='MOOG version')

    # For the EW method
    ew_parser = subparsers.add_parser('ew', parents=[parent_parser], help='EW method')
    ew_parser.add_argument('linelist',             help='Input linelist file',widget='FileChooser')
    ew_parser.add_argument('--spectralType',       help='Input spectral type (optional)', default=False, metavar='Spectral type')
    ew_parser.add_argument('--model',              help='Model atmosphere',    default='kurucz95', choices=['kurucz95', 'kurucz08', 'marcs', 'PHOENIX'])
    ew_parser.add_argument('--Fixteff',            help='Fix temperature',     action='store_true', metavar='Fix temperature')
    ew_parser.add_argument('--Fixgravity',         help='Fix gravity',         action='store_true', metavar='Fix gravity')
    ew_parser.add_argument('--FixFeH',             help='Fix metallicity',     action='store_true', metavar='Fix [Fe/H]')
    ew_parser.add_argument('--Fixmicroturbulence', help='Fix microturbulence', action='store_true', metavar='Fix microturbulence')
    ew_parser.add_argument('--Iterations',         help='Maximum number of iterations', default=160, type=int)
    ew_parser.add_argument('--weights',            help='Calculate the slopes of EP and RW with weights', type=str, default='null', choices=['null', 'median', 'sigma', 'mad'])
    ew_parser.add_argument('--EPslope',            help='EP slope to converge', default=0.001, type=float, metavar='EP slope')
    ew_parser.add_argument('--RWslope',            help='RW slope to converge', default=0.003, type=float, metavar='RW slope')
    ew_parser.add_argument('--Fedifference',       help='Difference between FeI and FeII', default='0.000',   type=float, metavar='|<Fel>-<Fell>|')
    ew_parser.set_defaults(func=ew)

    # For the synhtesis method
    synth_parser = subparsers.add_parser('synth', parents=[parent_parser], help='Synthesis method')
    synth_parser.add_argument('--test', help='this is test')
    synth_parser.set_defaults(func=synth)

    # For calculating the abundances
    abund_parser = subparsers.add_parser('abund', parents=[parent_parser], help='Abundances')
    abund_parser.set_defaults(func=abund)

    args = parser.parse_args()
    return args.func(args)


if __name__ == '__main__':
    main()
