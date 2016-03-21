#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import logging
import os
import numpy as np
from utils import GetModels, _update_par
from interpolation import interpolator
from interpolation import save_model
from utils import _run_moog, Readmoog
import pandas as pd


def save(dic, overwrite):
    """Write results"""
    linelist = dic.pop('linelist')
    teff = dic.pop('Temperature')
    logg = dic.pop('Gravity')
    feh = dic.pop('[Fe/H]')
    vt = dic.pop('microturbulence')
    elements = dic.keys()


    def _get_header(elements):
        """Get the header and append new elements to it"""
        with open('abundances.csv', 'r') as f:
            header = f.readline()
        header = header.strip('\n').split(',')
        for element in elements:
            if element not in header:
                header.append(element)
        return ','.join(header)


    if not os.path.isfile('abundances.csv'):
        # The file does not exists, so create it and return from here
        header = 'linelist,temperature,logg,[Fe/H],vt,' + ','.join(elements)
        with open('abundances.csv', 'w') as fout:
            fout.writelines(header + '\n')
            line = '%s,%i,%.2f,%.2f,%.2f' % (linelist, teff, logg, feh, vt)
            for element in elements:
                line += ',%s' % dic[element]
            fout.writelines(line + '\n')
        return
    else:
        header = _get_header(elements)

    if overwrite:
        with open('abundances.csv', 'w') as fout:
            fout.writelines(header + '\n')
    else:
        try:
            df = pd.read_csv('abundances.csv', na_values='...')
            header = header.split(',')
            for key in header:
                if key not in df.columns:
                    df[key] = np.nan
            df.to_csv(path_or_buf='abundances.csv', header=header, index=False, na_rep='...')
        except IOError:
            # It does not exists yet, so create the file from scratch
            with open('abundances.csv', 'w') as fout:
                fout.writelines(header + '\n')


    df = pd.read_csv('abundances.csv', na_values='...')
    rows = df.shape[0]
    for element in header[5:]:
        if element in elements:
            df[element] = dic[element]
        else:
            df[element] = np.nan
    df['linelist'] = linelist
    df['temperature'] = teff
    df['logg'] = logg
    df['[Fe/H]'] = feh
    df['vt'] = vt

    if rows:
        df.drop(df.index[range(rows-1)], inplace=True)

    df.to_csv(path_or_buf='abundances.csv', header=False, index=False, mode='a', na_rep='...')


def _options(options=False):
    '''Reads the options inside the config file'''
    defaults = {'model':'kurucz95',
                'MOOGv': 2014
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


def abundancedriver(starLines='StarMe.cfg', overwrite=False):
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
        for line in lines:
            abundance_dict = {}
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
            grid = GetModels(teff=initial[0], logg=initial[1], feh=initial[2], atmtype=options['model'])
            models, nt, nl, nf = grid.getmodels()
            logger.info('Initial interpolation of model...')
            inter_model = interpolator(models, teff=(initial[0], nt), logg=(initial[1], nl), feh=(initial[2], nf))
            save_model(inter_model, params=initial)
            logger.info('Interpolation successful.')
            _run_moog()

            m = Readmoog()
            elements, abundances = m.elements()
            abundance_dict={'linelist': line[0],
                            'Temperature': initial[0],
                            'Gravity': initial[1],
                            '[Fe/H]': initial[2],
                            'microturbulence': initial[3]}

            for element, abundance in zip(elements, abundances):
                abundance_dict[element] = abundance

            save(abundance_dict, overwrite=overwrite)


if __name__ == '__main__':
    abundancedriver()
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)
    df = pd.read_csv('abundances.csv')
    print df.to_string(index=False, justify='right')
    # print '\n\n ---------------'
    # print '|Some statistics|'
    # print ' ---------------\n'
    # print df.describe()
