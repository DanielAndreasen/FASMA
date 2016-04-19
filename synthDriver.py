#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import logging
import os
from shutil import copyfile
import yaml
import numpy as np
from utils import GetModels, _update_par_synth
from interpolation import interpolator
from interpolation import save_model
from utils import _run_moog
from observations import read_observations, plot_synth, plot_synth_obs, chi2
import seaborn

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


def _getMic(teff, logg, feh):
    """Calculate micro turbulence."""
    if logg >= 3.95:  # Dwarfs Tsantaki 2013
        mic = 6.932 * teff * (10**(-4)) - 0.348 * logg - 1.437
        return round(mic, 2)
    else:  # Giants Adibekyan 2015
        mic = 2.72 - (0.457 * logg) + (0.072 * feh)
        return round(mic, 2)


def _options(options=None):
    '''Reads the options inside the config file'''
    defaults = {'spt': False,
                'weights': 'null',
                'model': 'kurucz95',
                'MOOGv': 2014,
                'plotpars': 1,
                'plot': False,  #This is irrelevant with the batch.par value
                'start_wave': False, #if this is set by the user
                'end_wave': False,
                'step_wave': 0.01,
                'step_flux': 5.0,
                'iterations': 160,
                'observations': False,
                'resolution': 0.06,
                'vmac': 0.0,
                'vsini': 0.0,
                'limb': 0.0,
                'lorentz': 0.0
                }
    if not options:
        return defaults
    else:
        for option in options.split(','):
            if ':' in option:
                option = option.split(':')
                defaults[option[0]] = option[1]
            else:
                # Clever way to change the boolean
                if option in ['teff', 'logg', 'feh', 'vt']:
                    option = 'fix_%s' % option
                defaults[option] = False if defaults[option] else True
        defaults['model'] = defaults['model'].lower()
        defaults['iterations'] = int(defaults['iterations'])
        defaults['step_wave'] = float(defaults['step_wave'])
        defaults['step_flux'] = float(defaults['step_flux'])
        defaults['plotpars'] = int(defaults['plotpars'])
        #defaults['plot'] = int(defaults['plot'])
        #defaults['observations'] = str(defaults['observations'])
        defaults['resolution'] = float(defaults['resolution'])
        defaults['vmac'] = float(defaults['vmac'])
        defaults['vsini'] = float(defaults['vsini'])
        defaults['limb'] = float(defaults['limb'])
        defaults['lorentz'] = float(defaults['lorentz'])
        defaults['MOOGv'] = int(defaults['MOOGv'])
        return defaults

def read_wave(linelist): 
    """Read the wavelenth intervals of the line list"""

    with open(linelist, 'r') as f:

        lines = f.readlines()
    first_line = lines[0].split()

    if len(first_line) == 1: 
        start_wave = first_line[0].split('-')[0]
        end_wave = first_line[0].split('-')[1]
    else:
        start_wave = first_line[0]
        end_wave = lines[-1].split()[0]
    return start_wave, end_wave  


def synthdriver(starLines='StarMe.cfg', overwrite=False):
    """The function that glues everything together

    Input:
    starLines   -   Configuration file (default: StarMe.cfg)
    parfile     -   The configuration file for MOOG
    model       -   Type of model atmospheres
    plot        -   Plot results (currently not implemented)

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

    # Check if there is a directory called linelist, if not create it and ask the user to put files there
    if not os.path.isdir('linelist'):
        logger.error('Error: The directory linelist does not exist!')
        os.mkdir('linelist')
        logger.info('linelist directory was created')
        raise IOError('linelist directory did not exist! Put the linelists inside that directory, please.')

    # Create results directory
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

            # Check if the linelist is inside the directory if not log it and pass to next linelist
            if not os.path.isfile('linelist/%s' % line[0]):
                logger.error('Error: linelist/%s not found.' % line[0])
                parameters = None
                continue

            if len(line) == 1:
                initial = [5777, 4.44, 0.00, 1.00]
                options = _options()
                plot_flag = False
                x_obs, y_obs = (None, None)
                _update_par_synth('linelist/%s' % line[0], start_wave=start_wave, end_wave=end_wave, options=options)
                logger.info('Setting solar values {0}, {1}, {2}, {3}'.format(*initial))

            elif len(line) == 2:
                logger.info('Spectral type given: %s' % line[1])
                options = _options(line[1])
                if options['spt']:
                    Teff, logg = _getSpt(options['spt'])
                    mic = _getMic(Teff, logg)
                    initial = (Teff, logg, 0.00, mic)
                else:
                    initial = [5777, 4.44, 0.00, 1.00]

                if options['start_wave'] and options['end_wave']:
                    start_wave = options['start_wave']
                    end_wave = options['end_wave']
                else:
                    start_wave, end_wave = read_wave('linelist/%s' % line[0])

                if options['observations']:
                    plot_flag = True 
                    x_obs, y_obs = read_observations(fname=options['observations'], start_synth=start_wave, end_synth=end_wave)
                elif options['plot']:
                    plot_flag = True
                    x_obs, y_obs = (None, None)
                else:
                    plot_flag = False
                    x_obs, y_obs = (None, None)
                _update_par_synth('linelist/%s' % line[0], start_wave=start_wave, end_wave=end_wave, options=options)
                logger.info('Initial parameters: {0}, {1}, {2}, {3}'.format(*initial))

            elif len(line) == 5:
                logger.info('Initial parameters given by the user.')
                initial = map(float, line[1::])
                initial[0] = int(initial[0])
                options = _options()
                plot_flag = False
                x_obs, y_obs = (None, None)
                logger.info('Initial parameters: {0}, {1}, {2}, {3}'.format(*initial))
                _update_par_synth('linelist/%s' % line[0], start_wave=start_wave, end_wave=end_wave, options=options)

            elif len(line) == 6:
                logger.info('Initial parameters given by user.')
                initial = map(float, line[1:-1])
                initial[0] = int(initial[0])
                logger.info('Initial parameters: {0}, {1}, {2}, {3}'.format(*initial))
                options = _options(line[-1])
                if options['start_wave'] and options['end_wave']:
                    start_wave = options['start_wave']
                    end_wave = options['end_wave']
                else:
                    start_wave, end_wave = read_wave('linelist/%s' % line[0])
                if options['observations']:
                    plot_flag = True
                    x_obs, y_obs = read_observations(fname=options['observations'], start_synth=start_wave, end_synth=end_wave)
                elif options['plot']:
                    plot_flag = True
                    x_obs, y_obs = (None, None)
                else:
                    plot_flag = False
                    x_obs, y_obs = (None, None)
                _update_par_synth('linelist/%s' % line[0], start_wave=start_wave, end_wave=end_wave, options=options)

            else:
                logger.error('Could not process information for this line list: %s' % line)
                continue

            if options['model'] != 'kurucz95' and options['model'] != 'apogee_kurucz' and options['model'] != 'marcs':
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
            inter_model = interpolator(models,
                                       teff=(initial[0], nt),
                                       logg=(initial[1], nl),
                                       feh=(initial[2], nf))
            save_model(inter_model, params=initial)
            logger.info('Interpolation successful.')

            logger.info('Starting the minimization procedure...')

            # Options not in use will be removed
            if __name__ == '__main__':
                options['GUI'] = False  # Running batch mode
            else:
                options['GUI'] = True  # Running GUI mode
            options.pop('spt')

            # Setting the models to use
            if options['model'] != 'kurucz95' and options['model'] != 'kurucz08':
                logger.error('Your request for type: %s is not available' % model)
                continue

            # Get the initial grid models
            logger.info('Getting initial model grid')
            print('This is your synthetic spectrum: results/%s' % line[0] + '.spec')
    synth_out = 'results/%s' % line[0] + '.spec'

    return plot_flag, x_obs, y_obs, synth_out 

if __name__ == '__main__':
    plot_flag, x_obs, y_obs, synth_out = synthdriver(starLines='StarMe.cfg')
    _run_moog(driver='synth')
    if x_obs is not None:
        plot_synth_obs(x_obs, y_obs, fname=synth_out) #print obs and synth
        chi2(x_obs, y_obs, fname=synth_out)
    if plot_flag and x_obs is None: #print only synth
        plot_synth(fname=synth_out)    
      

