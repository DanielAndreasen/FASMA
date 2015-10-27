#!/usr/bin/env python
# -*- coding: utf8 -*-

from MOOGme import moogme
import argparse
from gooey import Gooey, GooeyParser

@Gooey(program_name='MOOG Made Easy - for deriving stellar parameters', default_size=(610, 930))
def _parser():
    '''Take care of all the argparse stuff.
    :returns: the args
    '''
    parser = GooeyParser(description='Set parameters')

    parser.add_argument('linelist',             help='Input linelist file',widget='FileChooser')
    parser.add_argument('--temperature',        help='Input initial temperature',      default=5777,  type=int)
    parser.add_argument('--surfacegravity',     help='Input initial gravity',          default=4.44,  type=float)
    parser.add_argument('--FeH',                help='Input initial metallicity',      default='0.0', type=float)
    parser.add_argument('--microturbulence',    help='Input initial microturbulence',  default=1.0,   type=float)
    parser.add_argument('--spectralType',       help='Input spectral type (optional)', default=False)
    parser.add_argument('--model',              help='Model atmosphere', default='kurucz95', choices=['kurucz95', 'kurucz08', 'marcs', 'PHOENIX'])
    parser.add_argument('--Fixteff',            help='Fix temperature',               action='store_true')
    parser.add_argument('--Fixgravity',         help='Fix gravity',         action='store_true')
    parser.add_argument('--FixFeH',             help='Fix metallicity',     action='store_true')
    parser.add_argument('--Fixmicroturbulence', help='Fix microturbulence', action='store_true')
    parser.add_argument('--Iterations',         help='Maximum number of iterations', default=160, type=int)
    parser.add_argument('--weights',            help='Calculate the slopes of EP and RW with weights', type=str, default='null', choices=['null', 'median', 'sigma', 'mad'])
    parser.add_argument('--EPslope',            help='EP slope to converge', default=0.003, type=float)
    parser.add_argument('--RWslope',            help='RW slope to converge', default=0.003, type=float)
    parser.add_argument('--Fedifference',       help='Difference between FeI and FeII', default='0.000',   type=float)
    parser.add_argument('--MOOGv',              help='Version of MOOG', default='2013', choices=['2013', '2014'], type=str)
    return parser.parse_args()


def create_starme(args):
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


def main():
    args = _parser()
    create_starme(args)
    moogme()


if __name__ == '__main__':
    main()
