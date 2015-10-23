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
from utils import error
from minimization import minimize


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
        cmd = 'cp summary.out results/%s.out' % linelist
    else:
        cmd = 'cp summary.out results/%s.NC.out' % linelist
        # os.system('cp minimization_profile.dat %s.profile.dat' % linelist)

    os.system(cmd)


def _options(options=False):
    '''Reads the options inside the config file'''
    defaults = {'spt': False,
                'weights': 'null',
                'plot': False,
                'model':'kurucz95',
                'teff': False,
                'logg': False,
                'feh': False,
                'vt': False,
                'iterations': 160,
                'EPslope': 0.001,
                'RWslope': 0.001,
                'Fedifference': 0.01
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
        defaults['model'] = defaults['model'].lower()
        defaults['iterations'] = int(defaults['iterations'])
        defaults['EPslope'] = float(defaults['EPslope'])
        defaults['RWslope'] = float(defaults['RWslope'])
        defaults['Fedifference'] = float(defaults['Fedifference'])
        return defaults


def moogme(starLines='StarMe.cfg'):
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
            if not os.path.isfile('linelist/%s' % line[0]):
                logger.error('Error: linelist/%s not found.' % line[0])
                parameters = None
                continue
            else:
                _update_par(line_list='linelist/%s' % line[0])

            if len(line) == 1:
                initial = (5777, 4.44, 0.00, 1.00)
                options = _options()
                logger.info('Setting solar values {0}, {1}, {2}, {3}'.format(*initial))

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
                else:
                    initial = (5777, 4.44, 0.00, 1.00)
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
            if options['model'] != 'kurucz95' and options['model'] != 'kurucz08':
                logger.error('Your request for type: %s is not available' % model)
                continue

            # Get the initial grid models
            logger.info('Getting initial model grid')
            # TODO: Fix the interpolation please!
            if initial[1] > 4.99:  # quick fix
                initial[1] = 4.99
            models, nt, nl, nf = _get_model(teff=initial[0], logg=initial[1], feh=initial[2], atmtype=options['model'])
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
            #                                  fix_feh=fix_feh, fix_vt=fix_vt,
            #                                  weights=options['weights'],
            #                                  iteration=options['iterations'],
            #                                  EPcrit=options['EPslope'],
            #                                  RWcrit=options['RWslope'],
            #                                  ABdiffcrit=options['Fedifference'])
            parameters, converged = minimize(initial, fun_moog_fortran, bounds='kurucz95',
                                             fix_teff=fix_teff, fix_logg=fix_logg,
                                             fix_feh=fix_feh, fix_vt=fix_vt,
                                             weights=options['weights'],
                                             iteration=options['iterations'],
                                             EPcrit=options['EPslope'],
                                             RWcrit=options['RWslope'],
                                             ABdiffcrit=options['Fedifference'])
            logger.info('Finished minimization procedure')
            _renaming(line[0], converged)
            parameters = error(line[0])
            with open('results.csv', 'a') as output:
                tmp = [line[0]] + list(parameters) + [converged]
                output.write('\t'.join(map(str, tmp))+'\n')
            logger.info('Saved results to: results.csv')



            if __name__ == '__main__':
                if converged:
                    print('\nCongratulation, you have won! Your final parameters are:')
                else:
                    print('\nSorry, you did not win. However, your final parameters are:')
                print(u' Teff:    %i\u00B1%i\n logg:    %.2f\u00B1%.2f\n [Fe/H]: %.2f\u00B1%.2f\n vt:      %.2f\u00B1%.2f\n' %
                     (parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5],parameters[6],parameters[7]))
            elif __name__ == 'MOOGme':
                if converged:
                    print('\nCongratulation, you have won! Your final parameters are:')
                else:
                    print('\nSorry, you did not win. However, your final parameters are:')
                print('Teff:      %i+/-%i\nlogg:    %.2f+/-%.2f\n[Fe/H]:  %.2f+/-%.2f\nvt:        %.2f+/-%.2f\n' %
                     (parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5],parameters[6],parameters[7]))
    return parameters


if __name__ == '__main__':
    parameters = moogme()
