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
# from utils import error
from minimization import minimize


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Plot fits file for ARES. Be'
                                     ' careful with large files')

    parser.add_argument('-i', '--configfile',
                        help='Input config file',
                        default='StarMe.cfg')
    parser.add_argument('-p', '--parfile',
                        help='The parameter file (default: batch.par)',
                        default='batch.par')
    parser.add_argument('-m', '--model',
                        help='Model atmosphere',
                        default='Kurucz95',
                        choices=['Kurucz95', 'Kn', 'Marcs', 'PHOENIX'])
    parser.add_argument('-pl', '--plot',
                        help='Plot the slopes',
                        default=False)
    parser.add_argument('-ol', '--outliers',
                        help='Remove n*sigma outliers',
                        type=float,
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
        mic = 6.932 * teff * (10**(-4)) - 0.348 * logg - 1.437
        return round(mic, 2)
    else:  # Giants
        mic = 3.7 - (5.1 * teff * (10**(-4)))
        return round(mic, 2)


def _renaming(linelist, converged):
    """Save the output in a file related to the linelist"""
    if converged:
        cmd = 'cp summary.out ./results/%s.out' % linelist
    else:
        cmd = 'cp summary.out ./results/%s.NC.out' % linelist
        # os.system('cp minimization_profile.dat %s.profile.dat' % linelist)

    os.system(cmd)


def _options(options=False):
    '''Reads the options inside the config file'''
    defaults = {'spt': False,
                'outlier': False,
                'plot': False,
                'models':'K95',
                'teff': False,
                'logg': False,
                'feh': False,
                'vt': False,
                'iterations': 25,
                'epslope': 0.001,
                'rwslope': 0.001,
                'abdiff': 0.01
                }
    if not options:
        return defaults
    else:
        options = options.split(',')
        for option in options:
            if ':' in option:
                option = option.split(':')
                defaults[option[0]] = option[1]
            else:
                defaults[option] = True
        return defaults


def moogme(starLines, parfile='batch.par', model='kurucz95',
           plot=False, outlier=False):
    """The function that glues everything together

    Input:
    starLines   -   Configuration file (default: StarMe.cfg)
    parfile     -   The configuration file for MOOG
    model       -   Type of model atmospheres
    plot        -   Plot results (currently not implemented)
    outlier     -   Remove outliers (currently not implemented)

    Output:
    <linelist>.(NC).out     -   NC=not converget.
    results.csv             -   Easy readable table with results from many linelists
    """
    try:  # Cleaning from previous runs
        os.remove('captain.log')
    except OSError:
        pass
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
        logger.error('%s does not exist' % parfile)
        raise IOError('%s does not exist' % parfile)
    if parfile != 'batch.par':
        logger.info('Copying %s to batch.par' % parfile)
        os.system('cp %s batch.par' % parfile)
        rm_batch = True

    if plot:
        logger.debug('plot keyword not implemented yet')
    if outlier:
        logger.debug('outlier keyword not implemented yet')

    with open('results.csv', 'w') as output:
        tmp = ['linelist', 'teff', 'tefferr', 'logg', 'loggerr', 'feh', 'feherr', 'vt', 'vterr', 'convergence']
        output.write('\t'.join(tmp)+'\n')

    #Check if there is a directory called linelist, if not create it and ask the user to put files there
    if not os.path.isdir('linelist'):
        logger.error('Error: The directory linelist does not exist!')
        os.mkdir('linelist')
        logger.info('linelist directory was created')
        raise IOError('linelist directory did not exist! Put the linelists inside that directory, please.')

    #Create results directory
    if not os.path.isdir('results'):
        os.mkdir('results')
        logger.info('results directory was created')

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

            #Check if the linelist is inside the directory if not log it and pass to next linelist
            if not os.path.isfile('./linelist/'+line[0]):
                logger.error('Error: The linelist has to be inside the directory linelist')
                parameters = None
                continue
            else:
                _update_par(line_list='./linelist/'+line[0])
            if len(line) == 1:
                initial = (5777, 4.44, 0.00, 1.00)
                options = _options()
                logger.info('Setting solar values {0}, {1}, {2}, {3}'.format(*initial))
                # Update batch.par

            elif len(line) == 5:
                logger.info('Initial parameters given by the user.')
                initial = map(float, line[1::])
                options = _options()
                logger.info('Initial parameters: {0}, {1}, {2}, {3}'.format(*initial))

            elif len(line) == 2:
                logger.info('Spectral type given: %s' % line[1])
                options = _options(line[1])
                if options['spt']:
                    Teff, logg = _getSpt(options['spt'])
                mic = _getMic(Teff, logg)
                initial = (Teff, logg, 0.00, mic)
                logger.info('Initial parameters: {0}, {1}, {2}, {3}'.format(*initial))
                fix_teff = options['teff']
                fix_logg = options['logg']
                fix_feh  = options['feh']
                fix_vt   = options['vt']
                if fix_teff:
                    logger.info('Effective temperature fixed at: %i' % initial[0])
                elif fix_logg:
                    logger.info('Surface gravity fixed at: %s' % initial[1])
                elif fix_feh:
                    logger.info('Metallicity fixed at: %s' % initial[2])
                elif fix_vt:
                    logger.info('Micro turbulence fixed at: %s' % initial[3])

            elif len(line) == 6:
                logger.info('Initial parameters given by user and some parameters fixed.')
                initial = map(float, line[1:-1])
                logger.info('Initial parameters: {0}, {1}, {2}, {3}'.format(*initial))
                options = _options(line[-1])
                fix_teff = options['teff']
                fix_logg = options['logg']
                fix_feh  = options['feh']
                fix_vt   = options['vt']
                if fix_teff:
                    logger.info('Effective temperature fixed at: %i' % initial[0])
                elif fix_logg:
                    logger.info('Surface gravity fixed at: %s' % initial[1])
                elif fix_feh:
                    logger.info('Metallicity fixed at: %s' % initial[2])
                elif fix_vt:
                    logger.info('Micro turbulence fixed at: %s' % initial[3])

            else:
                logger.error('Could not process information for this line list: %s' % line)
                continue

            # Setting the models to use
            if options['models'] != 'K95':
                raise NotImplementedError('Your request for type: %s is not available' % model)

            # Get the initial grid models
            logger.info('Getting initial model grid')
            if initial[1] > 4.99:  # quick fix
                initial[1] = 4.99
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

            with open('results.csv', 'a') as output:
                tmp = [line[0], parameters[0], 0, parameters[1], 0.00, parameters[2], 0.00, parameters[3], 0.00, converged]
                output.write('\t'.join(map(str, tmp))+'\n')
            logger.info('Saved results to: results.csv')

            print('\nCongratulation, you have won! Your final parameters are\n' + ' '.join(map(str,parameters)))
            print(line[0])
            # error(line[0])
    return parameters


if __name__ == '__main__':
    args = _parser()
    parameters = moogme(args.configfile, args.parfile, args.model,
                        args.plot, args.outliers)
