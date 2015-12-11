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
from interpolation import interpolator
from interpolation import save_model
from utils import fun_moog, fun_moog_fortran, _run_moog, _read_moog
from utils import error


def _renaming(linelist):
    """Save the output in a file related to the linelist"""
    cmd = 'cp summary.out results/%s.out' % linelist
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
        return defaults

def moogme_ab(starLines='StarMe.cfg'):
    """The function that glues everything together

    Input:
    starLines   -   Configuration file (default: StarMe.cfg)
    parfile     -   The configuration file for MOOG
    model       -   Type of model atmospheres
    plot        -   Plot results (currently not implemented)
    outlier     -   Remove outliers (currently not implemented)

    Output:
    <linelist>.out          -   Output file
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
            _run_moog()
            element, EP_slopes, RW_slopes, abundances = _read_moog()
            for e, a in zip(element, abundances):
                print('Element: ', e, 'Abundance:', a)
            _renaming(line[0])
    return 

if __name__ == '__main__':
    moogme_ab()
