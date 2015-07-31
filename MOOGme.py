#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import logging
import os
import argparse
try:
    import seaborn as sns
    sns.set_style('dark')
    sns.set_context('talk')
except ImportError:
    print('install seaborn for nice plots (pip install seaborn)')
import yaml

from utils import _get_model, _update_par
from model_interpolation import interpolator
from model_interpolation import save_model
from utils import fun_moog, fun_moog_fortran
from moog_minimization import minimize


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Plot fits file for ARES. Be'
                                     ' careful with large files')

    parser.add_argument('linelist', help='Input line list')
    parser.add_argument('-p', '--parfile',
                        help='The parameter file (default: batch.par)',
                        default='batch.par')
    parser.add_argument('-m', '--model',
                        help='Model atmosphere',
                        default='Kurucz95',
                        choices=['Kurucz95', 'Kn', 'Marcs', 'PHOENIX'])
    parser.add_argument('-i', '--initial',
                        help='Initial conditions (Teff, logg, [Fe/H], vt)',
                        nargs='+',
                        type=float,
                        default=False)
    parser.add_argument('-pl', '--plot',
                        help='Plot the slopes',
                        default=False)
    parser.add_argument('-ol', '--outliers',
                        help='Remove n*sigma outliers',
                        type=float,
                        default=False)
    parser.add_argument('-spt', '--spectralType',
                        help='Input spectral type (e.g. F4V) and get initial parameters',
                        default=False)
    parser.add_argument('-v', '--verbose',
                        help='Print information to the screen along the way',
                        default=False)
    return parser.parse_args()


def _getSpt(spt):
    """Get the spectral type from a string like 'F5V'."""
    if len(spt) > 4:
        raise ValueError('Spectral type most be of the form: F8V')
    if '.' in spt:
        raise ValueError('Do not use half spectral types as %s' % spt)
    with open('SpectralTypes.yml', 'r') as f:
        d = yaml.safe_load(f)
    temp = spt[0:2]
    lum = spt[2:]
    try:
        line = d[lum][temp]
    except KeyError:
        print('Was not able to find the spectral type: %s' % spt)
        print('Teff=5777 and logg=4.44')
        return 5777, 4.44
    try:
        line = line.split()
        teff = int(line[0])
        logg = float(line[1])
    except AttributeError:
        teff = line
        logg = 4.44
    return teff, logg


def _getMic(teff, logg):
    """Calculate micro turbulence. REF? Doyle 2014"""
    if logg >= 3.95:  # Dwarfs
        mic = 6.932 * teff * (10 ** (-4)) - 0.348 * logg - 1.437
        return round(mic, 2)
    else:  # Giants
        mic = 3.7 - (5.1 * teff * (10 ** (-4)))
        return round(mic, 2)


def _renaming(linelist, converged):
    """Save the output in a file related to the linelist"""
    if converged:
        cmd = 'cp summary.out %s.out' % linelist
    else:
        cmd = 'cp summary.out %s.NC.out' % linelist
        # os.system('cp minimization_profile.dat %s.profile.dat' % linelist)

    os.system(cmd)


def moogme(starLines, parfile='batch.par', model='kurucz95',
           initial=False, plot=False, outlier=False, spt=False):
    """Some doc"""
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler('captain.log')
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # TODO: Let us just insist on having a batch.par in the directory.
    # Create one if necessary
    # Preparing the batch file for MOOGSILENT
    if not os.path.isfile(parfile):
        raise IOError('%s does not exist' % parfile)
    if parfile != 'batch.par':
        os.system('cp %s batch.par' % parfile)
        rm_batch = True

    with open(starLines, 'r') as lines:
        for line in lines:
            if not line[0].isalpha():
                logger.debug('Skipping header: %s' % line.strip())
                continue
            logger.info('Line list: %s' % line.strip())
            fix_teff = False
            fix_logg = False
            fix_feh = False
            fix_vt = False
            line = line.strip()
            line = line.split(' ')
            if len(line) == 1:
                initial = (5777, 4.44, 0.00, 1.00)
                logger.info('Setting solar values {0}, {1}, {2}, {3}'.format(*initial))
                # Update batch.par
                _update_par(line_list=line[0])
            elif len(line) == 2:
                logger.info('Spectral type given: %s' % line[1])
                spt = line[1]
                Teff, logg = _getSpt(spt)
                mic = _getMic(Teff, logg)
                initial = (Teff, logg, 0.00, mic)
                logger.info('Initial parameters: {0}, {1}, {2}, {3}'.format(*initial))
                _update_par(line_list=line[0])
            elif len(line) == 5:
                logger.info('Initial parameters given by the user.')
                initial = map(float, line[1::])
                logger.info('Initial parameters: {0}, {1}, {2}, {3}'.format(*initial))
                _update_par(line_list=line[0])
            elif len(line) == 6:
                logger.info('Initial parameters given by user and some parameters fixed.')
                initial = map(float, line[1:-1])
                logger.info('Initial parameters: {0}, {1}, {2}, {3}'.format(*initial))
                fix = line[-1].lower()
                fix_teff = True if 'teff' in fix else False
                fix_logg = True if 'logg' in fix else False
                fix_feh  = True if 'feh'  in fix else False
                fix_vt   = True if 'vt'   in fix else False
                if fix_teff:
                    logger.info('Effective temperature fixed at: %i' % initial[0])
                elif fix_logg:
                    logger.info('Surface gravity fixed at: %s' % initial[1])
                elif fix_feh:
                    logger.info('Metallicity fixed at: %s' % initial[2])
                elif fix_vt:
                    logger.info('Micro turbulence fixed at: %s' % initial[3])
                _update_par(line_list=line[0])
            else:
                logger.error('Could not process information for this line list.')
                continue

            # Setting the models to use
            if model != 'Kurucz95':
                raise NotImplementedError('Your request for type: %s is not available' % model)

            # Get the initial grid models
            logger.info('Getting initial model grid')
            models, nt, nl, nf = _get_model(teff=initial[0], logg=initial[1], feh=initial[2])
            logger.info('Initial interpolation of model...')
            inter_model = interpolator(models,
                                       teff=(initial[0], nt),
                                       logg=(initial[1], nl),
                                       feh=(initial[2], nf))
            save_model(inter_model, params=initial)
            logger.info('Interpolation successful.')

            logger.info('Starting the minimization procedure...')
            # parameters, converged = minimize(initial, fun_moog, bounds=model,
            #                                  fix_teff=fix_teff, fix_logg=fix_logg,
            #                                  fix_feh=fix_feh, fix_vt=fix_vt)
            parameters, converged = minimize(initial, fun_moog_fortran, bounds=model,
                                             fix_teff=fix_teff, fix_logg=fix_logg,
                                             fix_feh=fix_feh, fix_vt=fix_vt)
            logger.info('Finished minimization procedure')
            _renaming(line[0], converged)

            print('\nCongratulation, you have won! Your final parameters are\n' + ', '.join(map(str,parameters)))
            print(line[0])
    return parameters


if __name__ == '__main__':
    args = _parser()
    parameters = moogme(args.linelist, args.parfile, args.model,
                        args.initial, args.plot,
                        args.outliers, args.spectralType)
