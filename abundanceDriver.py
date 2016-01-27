#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import logging
import os
import numpy as np

from utils import _get_model, _update_par
from interpolation import interpolator
from interpolation import save_model
from utils import _run_moog, _read_moog


def save(dic):
    """Write results"""
    linelists = dic.pop('linelist')
    teff = dic.pop('Temperature')
    logg = dic.pop('Gravity')
    feh = dic.pop('[Fe/H]')
    vt = dic.pop('microturbulence')
    elements = dic.keys()
    header = 'linelist,temperature,logg,[Fe/H],vt,' + ','.join(elements)

    with open('abundances.csv', 'w') as fout:
        fout.writelines(header + '\n')
        for i, linelist in enumerate(linelists):
            line = '%s,%i,%.2f,%.2f,%.2f' % (linelist, teff[i], logg[i], feh[i], vt[i])
            for element in elements:
                line += ',%s' % dic[element][i]
            fout.writelines(line + '\n')


def _options(options=False):
    '''Reads the options inside the config file'''
    defaults = {'model':'kurucz95',
                'MOOGv': 2013
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


def abundancedriver(starLines='StarMe.cfg'):
    """The function that glues everything together

    Input:
    starLines   -   Configuration file (default: StarMe.cfg)

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
        # TODO: We should not remove this file from previous runs. New standard is to append.
        if os.path.isfile('results.csv'):
            os.remove('results.csv')
        counter = 0
        abundance_dict = {}
        for line in lines:
            if not line[0].isalpha():
                logger.debug('Skipping header: %s' % line.strip())
                continue
            logger.info('Line list: %s' % line.strip())
            line = line.strip()
            line = line.split(' ')

            # Check if the linelist is inside the directory if not log it and pass to next linelist
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

            elif len(line) == 6:
                logger.info('Initial parameters given by the user.')
                initial = map(float, line[1:-1])
                options = _options(line[-1])
                logger.info('Initial parameters: {0}, {1}, {2}, {3}'.format(*initial))
            else:
                logger.error('Could not process information for this line list: %s' % line)
                continue

            # Setting the models to use
            if (options['model'] != 'kurucz95') and (options['model'] != 'kurucz08'):
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


            elements, EP_slopes, RW_slopes, abundances = _read_moog()
            if counter == 0:
                abundance_dict={'linelist':[line[0]],
                                'Temperature':[initial[0]],
                                'Gravity': [initial[1]],
                                '[Fe/H]':[initial[2]],
                                'microturbulence':[initial[3]]}
            else:
                abundance_dict['linelist'].append(line[0])
                abundance_dict['Temperature'].append(initial[0])
                abundance_dict['Gravity'].append(initial[1])
                abundance_dict['[Fe/H]'].append(initial[2])
                abundance_dict['microturbulence'].append(initial[3])

            N = len(abundance_dict['linelist'])
            for element, abundance in zip(elements, abundances):
                print('Element: ', element, 'Abundance:', abundance)
                if element in abundance_dict.keys():
                    abundance_dict[element].append(abundance)
                else:
                    abundance_dict[element]=[np.nan]*counter+[abundance]

            for key in abundance_dict.keys():
                if len(abundance_dict[key]) < N:
                    abundance_dict[key].append(np.nan)
            counter += 1

    save(abundance_dict)


if __name__ == '__main__':
    abundancedriver()
