#!/usr/bin/env python
# -*- coding: utf8 -*-

from __future__ import print_function
from ewDriver import ewdriver
from synthDriver import synthdriver
from abundanceDriver import abundancedriver
from aresDriver import aresdriver
from gooey import Gooey, GooeyParser


def ew(args):
    '''Driver for the EW method'''
    fout = ''
    linelist = args.linelist.rpartition('/')[-1]
    args.outlier = False if args.outlier == 'False' else args.outlier
    if args.spectralType:
        if not args.temperature and not args.surfacegravity:
            fout += '%s spt:%s,' % (linelist, args.spectralType)
        else:
            print('Temperature and/or surface gravity set. Will not use spectral type.')
    else:
        if not args.temperature:
            print('Warning: Temperature was set to 5777K')
            args.temperature = 5777
        if not args.surfacegravity:
            print('Warning: Surface gravity was set to 4.44')
            args.surfacegravity = 4.44
        if not args.FeH:
            args.FeH = 0.0
        if not args.microturbulence:
            print('Warning: Microturbulence was set to 1.00')
            args.microturbulence = 1.00
        fout += '%s %s %s %s %s ' % (linelist, args.temperature, args.surfacegravity, args.FeH, args.microturbulence)

    fout += 'model:%s,iterations:%s,weights:%s,RWcrit:%s,EPcrit:%s,ABdiffcrit:%s,MOOGv:%s,sigma:%s' % (args.model, args.Iterations, args.weights, args.RWslope, args.EPslope, args.Fedifference, args.MOOGv, args.sigma)
    if args.outlier:
        fout += ',outlier:%s' % args.outlier
    if args.refine:
        fout += ',refine'
    if args.Fixteff:
        fout += ',teff'
    if args.Fixgravity:
        fout += ',logg'
    if args.FixFeH:
        fout += ',feh'
    if args.Fixmicroturbulence:
        fout += ',vt'
    if args.teffrange:
        fout += ',teffrange'
    if args.autofixvt:
        fout += ',autofixvt'
    if args.tmcalc:
        fout += ',tmcalc'
    with open('StarMe_ew.cfg', 'w') as f:
        f.writelines(fout)
    ewdriver(overwrite=args.overwrite)


def synth(args):
    """Driver for the synthesis method"""
    fout = ''
    linelist = args.linelist.rpartition('/')[-1]
    if args.spectralType:
        if not args.temperature and not args.surfacegravity:
            fout += '%s spt:%s,' % (linelist, args.spectralType)
        else:
            print('Temperature and/or surface gravity set. Will not use spectral type.')
    else:
        if not args.temperature:
            print('Warning: Temperature was set to 5777K')
            args.temperature = 5777
        if not args.surfacegravity:
            print('Warning: Surface gravity was set to 4.44')
            args.surfacegravity = 4.44
        if not args.FeH:
            args.FeH = 0.0
        if not args.microturbulence:
            print('Warning: Microturbulence was set to 1.00')
            args.microturbulence = 1.00
        if not args.macroturbulence:
            print('Warning: Macroturbulence was set to 3.21')
            args.macroturbulence = 3.21
        if not args.rotation:
            print('Warning: Rotation was set to 1.90')
            args.microturbulence = 1.90
        fout += '%s %s %s %s %s %s %s ' % (linelist, args.temperature, args.surfacegravity, args.FeH, args.microturbulence, args.macroturbulence, args.rotation)

    fout += 'model:%s,MOOGv:%s,step_wave:%s,step_flux:%s,inter_file:%s,limb:%s,damping:%s' % (args.model, args.MOOGv, args.step_wave, args.step_flux, args.inter_file, args.limb, args.damping)
    if args.observations:
        fout += ',observations:%s' % args.observations
    if args.resolution:
        fout += ',resolution:%s' % args.resolution
    if args.snr:
        fout += ',snr:%s' % args.snr
    if args.plot:
        fout += ',plot'
    if args.plot_res:
        fout += ',plot_res'
    if args.save:
        fout += ',save'
    if args.minimize:
        fout += ',minimize'
    if args.Fixteff:
        fout += ',teff'
    if args.Fixgravity:
        fout += ',logg'
    if args.FixFeH:
        fout += ',feh'
    if args.Fixmicroturbulence:
        fout += ',vt'
    if args.Fixmacroturbulence:
        fout += ',vmac'
    if args.Fixrotation:
        fout += ',vsini'
    if args.flag_vt:
        fout += ',flag_vt'
    if args.flag_vmac:
        fout += ',flag_vmac'
    with open('StarMe_synth.cfg', 'w') as f:
        f.writelines(fout)
    synthdriver(overwrite=args.overwrite)


def abund(args):
    """Driver for abundances"""

    fout = ''
    linelist = args.linelist.rpartition('/')[-1]

    if not args.temperature:
        print('Warning: Temperature was set to 5777K')
        args.temperature = 5777
    if not args.surfacegravity:
        print('Warning: Surface gravity was set to 4.44')
        args.surfacegravity = 4.44
    if not args.FeH:
        args.FeH = 0.0
    if not args.microturbulence:
        print('Warning: Microturbulence was set to 1.00')
        args.microturbulence = 1.00
    fout += '%s %s %s %s %s ' % (linelist, args.temperature, args.surfacegravity, args.FeH, args.microturbulence)

    fout += 'model:%s,MOOGv:%s\n' % (args.model, args.MOOGv)

    with open('StarMe_abund.cfg', 'w') as f:
        f.writelines(fout)

    abundancedriver(overwrite=args.overwrite)


def ares(args):
    """Driver for ARES"""
    def rejt_from_snr(snr):
        """Calculate rejt from SNR"""
        return 1-1/snr

    if args.SNR:
        rejt = rejt_from_snr(args.SNR)
    else:
        rejt = args.rejt
    rejt = 0.999 if rejt > 0.999 else rejt
    plot = 1 if args.plots else 0
    rvmask = False if args.RVmask == 0 else args.RVmask
    out = args.spectrum + '.ares' if not args.output else args.output

    # Make the StarMe_ares.cfg file from Gooey
    fout = args.linelist.rpartition('/')[2]
    fout += ' %s' % args.spectrum
    fout += ' lambdai:%s,lambdaf:%s,smoothder:%s' % (args.lambdai, args.lambdaf, args.smoothder)
    fout += ',space:%s,lineresol:%s' % (args.space, args.lineresol)
    fout += ',miniline:%s,EWcut:%s' % (args.miniline, args.EWcut)
    if args.SNR:
        fout += ',rejt:%s' % rejt
    else:
        fout += ',rejt:%s' % args.rejt
    if args.plots:
        fout += ',plots_flag:1'
    if args.output:
        fout += ',output:%s' % args.output
    if rvmask:
        fout += ',rvmask:%s' % rvmask
    if args.force:
        fout += ',force'

    with open('StarMe_ares.cfg', 'w') as f:
        f.writelines(fout)

    aresdriver()
    if args.output:
        print('Congratulations! The final line list are here: linelist/%s.moog' % args.output)
    else:
        linelist_out = 'linelist/%s.moog' % args.spectrum.rpartition('/')[2].rpartition('.')[0]
        print('Congratulations! The final line list are here: %s' % linelist_out)


@Gooey(program_name='FASMA - spectral analysis',
       default_size=(900, 1000),
       image_dir='./img')
def main():
    '''Take care of all the argparse stuff.
    :returns: the args
    '''
    parser = GooeyParser(description='Set parameters')

    subparsers = parser.add_subparsers()

    # Common to all
    parent_parser = GooeyParser(add_help=False)
    parent_parser.add_argument('--temperature',     help='(in K)',    default=5777,  type=int,   metavar='Temperature')
    parent_parser.add_argument('--surfacegravity',  help='(in dex)',  default=4.44,  type=float, metavar='logg')
    parent_parser.add_argument('--FeH',             help='(in dex)',  default='0.00',type=float, metavar='[Fe/H]')
    parent_parser.add_argument('--microturbulence', help='(in km/s)', default=1.0,   type=float, metavar='Microturbulence')
    parent_parser.add_argument('--MOOGv',           default='2014', choices=['2013', '2014', '2016'], type=str, metavar='MOOG version')
    parent_parser.add_argument('--model',           help='Grid of models', default='kurucz95', choices=['kurucz95', 'apogee_kurucz', 'kurucz08', 'marcs', 'PHOENIX'], metavar='Model atmosphere')

    # For the EW method
    ew_parser = subparsers.add_parser('ew', parents=[parent_parser], help='EW method')
    ew_parser.add_argument('linelist',             help='Input linelist file', widget='FileChooser')
    ew_parser.add_argument('--spectralType',       help='(optional)', default=False, metavar='Spectral type')
    ew_parser.add_argument('--Fixteff',            help='Fix Teff',   action='store_true', metavar='Fix temperature')
    ew_parser.add_argument('--Fixgravity',         help='Fix logg',   action='store_true', metavar='Fix gravity')
    ew_parser.add_argument('--FixFeH',             help='Fix [Fe/H]', action='store_true', metavar='Fix metallicity')
    ew_parser.add_argument('--Fixmicroturbulence', help='Fix vt',     action='store_true', metavar='Fix microturbulence')
    ew_parser.add_argument('--tmcalc',             help='Better guess on initial conditions',     action='store_true', metavar='Set initial conditions')
    ew_parser.add_argument('--refine',             help='Refine parameters',   action='store_true', metavar='Refine parameters')
    ew_parser.add_argument('--Iterations',         help='Maximum number of iterations', default=160, type=int)
    ew_parser.add_argument('--outlier',            help='Remove outliers', default='False', choices=['False', '1Iter', '1Once', 'allIter', 'allOnce'])
    ew_parser.add_argument('--weights',            help='Calculate the slopes of EP and RW with weights', type=str, default='null', choices=['null', 'sigma', 'mad'])
    ew_parser.add_argument('--EPslope',            help='EP slope to converge', default=0.001, type=float, metavar='EP slope')
    ew_parser.add_argument('--RWslope',            help='RW slope to converge', default=0.003, type=float, metavar='RW slope')
    ew_parser.add_argument('--Fedifference',       help='Difference between FeI and FeII', default='0.000', type=float, metavar='|Fel-Fell|')
    ew_parser.add_argument('--overwrite',          help='Overwrite results.csv', action='store_true', default=False)
    ew_parser.add_argument('--teffrange',          help='Give warning at high Teff, and remove lines at low Teff', action='store_true', default=False, metavar='Teff range')
    ew_parser.add_argument('--autofixvt',          help='Auto fix vt if it goes too low or too high', action='store_true', default=False)
    ew_parser.add_argument('--sigma',              help='Number of sigma, for sigma clipping outliers (default: 3)', type=float, default=3)
    ew_parser.set_defaults(driver=ew)

    # For the synhtesis method
    synth_parser = subparsers.add_parser('synth', parents=[parent_parser], help='Synthesis method')
    synth_parser.add_argument('linelist',                help='Input linelist file', widget='FileChooser')
    synth_parser.add_argument('--macroturbulence',    help='(in km/s)',  default=3.21, type=float, metavar='Macroturbulence')
    synth_parser.add_argument('--rotation',           help='(in km/s)',  default=1.90, type=float, metavar='Rotational velocity')
    synth_parser.add_argument('--spectralType',       help='(optional)', default=False, metavar='Spectral type')
    synth_parser.add_argument('--minimize',           help='Start minimization', action='store_true', metavar='Minimization procedure')
    synth_parser.add_argument('--Fixteff',            help='Fix Teff',   action='store_true', metavar='Fix temperature')
    synth_parser.add_argument('--Fixgravity',         help='Fix logg',   action='store_true', metavar='Fix gravity')
    synth_parser.add_argument('--FixFeH',             help='Fix [Fe/H]', action='store_true', metavar='Fix metallicity')
    synth_parser.add_argument('--Fixmicroturbulence', help='Fix vt',     action='store_true', metavar='Fix microturbulence')
    synth_parser.add_argument('--Fixmacroturbulence', help='Fix vmac',   action='store_true', metavar='Fix macroturbulence')
    synth_parser.add_argument('--Fixrotation',        help='Fix vsini',  action='store_true', metavar='Fix rotation')
    synth_parser.add_argument('--flag_vt',            help='For each iteration', action='store_true', metavar='Change vt in minimization')
    synth_parser.add_argument('--flag_vmac',          help='For each iteration', action='store_true', metavar='Change vmac in minimization')
    synth_parser.add_argument('--step_wave',          help='(in Angstroms)', default=0.01, type=float, metavar='Wavelength step for synthesis')
    synth_parser.add_argument('--step_flux',          help='(in Angstroms)', default=10.0, type=float, metavar='Flux step for synthesis')
    synth_parser.add_argument('--inter_file',         help='File with the intervals', metavar='Intervals', widget='FileChooser')
    synth_parser.add_argument('--limb',               help='Coefficient for vsini broadening', default=0.6, type=float, metavar='Limb darkening')
    synth_parser.add_argument('--damping',            help='Van der waals damping', default='1', metavar='Damping option',  choices=['0', '1', '2'])
    synth_parser.add_argument('--observations',       help='File with the observations', widget='FileChooser', metavar='Observations')
    synth_parser.add_argument('--resolution',         help='Instrumental resolution', default=None, metavar='Resolution')
    synth_parser.add_argument('--snr',                help='Signal-to-noise ratio',   default=None, metavar='SNR')
    synth_parser.add_argument('--save',               help='Save spectrum', action='store_true', metavar='Save output')
    synth_parser.add_argument('--plot',               help='Plot option',   action='store_true', metavar='Plot spectrum')
    synth_parser.add_argument('--plot_res',           help='Plot option',   action='store_true', metavar='Plot residuals')
    synth_parser.add_argument('--overwrite',          help='Overwrite results.csv', action='store_true', default=False)
    synth_parser.set_defaults(driver=synth)

    # For calculating the abundances
    abund_parser = subparsers.add_parser('abund', parents=[parent_parser], help='Abundances')
    abund_parser.add_argument('linelist',    help='Input linelist file', widget='FileChooser')
    abund_parser.add_argument('--overwrite', help='Overwrite abundances.csv', action='store_true', default=False)
    abund_parser.set_defaults(driver=abund)

    # Driver for ARES
    ares_parser = subparsers.add_parser('ares', help='ARES')
    ares_parser.add_argument('linelist',    help='Input linelist file', widget='FileChooser')
    ares_parser.add_argument('spectrum',    help='1D spectrum',         widget='FileChooser')
    ares_parser.add_argument('--output',    help='Name of the output file (optional)', metavar='Output')
    ares_parser.add_argument('--lambdai',   help='Start of wavelength interval',  default=3900,  type=int)
    ares_parser.add_argument('--lambdaf',   help='End of wavelength interval',        default=25000, type=int)
    ares_parser.add_argument('--smoothder', help='Noise smoother',                    default=4,     type=int)
    ares_parser.add_argument('--space',     help='Interval for the line computation', default=2.0,   type=float)
    ares_parser.add_argument('--rejt',      help='Continuum position',                default=False, type=float)
    ares_parser.add_argument('--lineresol', help='Line resolution',                   default=0.07,  type=float)
    ares_parser.add_argument('--miniline',  help='Weaker line to be printed out',     default=5.0,     type=int)
    ares_parser.add_argument('--SNR',       help='If specified, the rejt is calculated', type=float)
    ares_parser.add_argument('--EWcut',     help='Cut for the maximum EW value',      default=200.0, type=float)
    ares_parser.add_argument('--RVmask',    help='RV of spectrum (km/s)',             default='0.0', type=float)
    ares_parser.add_argument('--force',     help='Remove troublesome lines automatically', default=False, action='store_true', metavar='Force finish')
    ares_parser.set_defaults(driver=ares)

    args = parser.parse_args()
    return args.driver(args)


if __name__ == '__main__':
    main()
