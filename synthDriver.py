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
from utils import fun_moog as func
from interpolation import interpolator, save_model
from observations import read_obs_intervals, plot, chi2
from synthetic import save_synth_spec,read_linelist
import seaborn
from minimization import minimize_synth

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
                'model': 'kurucz95',
                'MOOGv': 2014,
                'plotpars': 1,
                'plot': False,  #This is irrelevant with the batch.par value
                'step_wave': 0.01,
                'step_flux': 5.0,
                'minimize' : False,
                'observations' : False,
                'resolution': None,
                'vmac': 0.0,
                'vsini': 0.0,
                'limb': 0.6,
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
        defaults['step_wave'] = float(defaults['step_wave'])
        defaults['step_flux'] = float(defaults['step_flux'])
        defaults['plotpars'] = int(defaults['plotpars'])
        #defaults['plot'] = int(defaults['plot'])
        #defaults['observations'] = str(defaults['observations'])
        #defaults['resolution'] = int(defaults['resolution'])
        defaults['vmac'] = float(defaults['vmac'])
        defaults['vsini'] = float(defaults['vsini'])
        defaults['limb'] = float(defaults['limb'])
        defaults['lorentz'] = float(defaults['lorentz'])
        defaults['MOOGv'] = int(defaults['MOOGv'])
        return defaults


def synthdriver(starLines='StarMe_synth.cfg', overwrite=False):
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

            #Check if configuration parameters are correct
            if len(line) not in [1, 2, 5, 6]:
                logger.error('Could not process this information: %s' % line)
                continue

            # Check if the linelist is inside the directory, if not log it and pass to next linelist
            if not os.path.isfile('linelist/%s' % line[0]):
                logger.error('Error: linelist/%s not found.' % line[0])
                continue

            if len(line) == 1:
                initial = [5777, 4.44, 0.00, 1.00]
                options = _options()
                x_obs, y_obs = (None, None)
                #Create initial synthetic model
                N, r, fout = read_linelist(line[0])
                logger.info('Getting initial model grid')
                x_initial, y_initial = func(initial, atmtype=options['model'], driver='synth', version=options['MOOGv'], N=N, r=r, fout=fout, options=options)
                logger.info('Save initial synthetic spectrum')
                save_synth_spec(x_initial, y_initial, fname='initial.spec')
                logger.info('Interpolation successful.')
                logger.info('Setting solar values {0}, {1}, {2}, {3}'.format(*initial))

            elif len(line) == 2:
                options = _options(line[1])
                if options['spt']:
                    logger.info('Spectral type given: %s' % line[1])
                    Teff, logg = _getSpt(options['spt'])
                    mic = _getMic(Teff, logg)
                    initial = (Teff, logg, 0.00, mic)
                else:
                    initial = [5777, 4.44, 0.00, 1.00]

                #Create initial synthetic model
                N, r, fout = read_linelist(line[0])
                logger.info('Getting initial model grid')
                x_initial, y_initial = func(initial, atmtype=options['model'], driver='synth', version=options['MOOGv'], N=N, r=r, fout=fout, options=options)
                logger.info('Save initial synthetic spectrum')
                save_synth_spec(x_initial, y_initial, fname='initial.spec')
                logger.info('Interpolation successful.')
                logger.info('Setting solar values {0}, {1}, {2}, {3}'.format(*initial))

                if options['observations']:
                    # Check if observations exit, if not pass another line
                    if (not (os.path.isfile('spectra/%s' % options['observations'])) and (not os.path.isfile(options['observations']))):
                        logger.error('Error: %s not found.' % options['observations'])
                        continue

                    print('This is your observed spectrum: %s' % options['observations'])
                    x_obs, y_obs = read_obs_intervals('spectra/%s' % options['observations'], N, r)
                    if options['minimize']:
                        print('Starting minimization...')
                        logger.info('Starting the minimization procedure...')
                        params, x_final, y_final = minimize_synth(initial, x_obs, y_obs, N, r, fout, options)
                        logger.info('Minimization done.')

                        #Some statistics. Here synthetic spectrum is interpolated to the observed
                        chi2(x_obs, y_obs, x_initial, y_initial)
                        chi2(x_obs, y_obs, x_final, y_final)

                else:
                    x_obs, y_obs = (None, None)
                    x_final, y_final = (None, None)

                if options['plot']: #if there in no observed only the synthetic will be plotted
                    plot(x_obs, y_obs, x_initial, y_initial)
                    if options['minimize']:
                        plot(x_obs, y_obs, x_initial, y_initial)
                        plot(x_obs, y_obs, x_final, y_final)

            elif len(line) == 5:
                logger.info('Initial parameters given by the user.')
                initial = map(float, line[1::])
                initial[0] = int(initial[0])
                options = _options()
                x_obs, y_obs = (None, None) #No observed spectrum, no plots

                #Create initial synthetic model
                N, r, fout = read_linelist(line[0])
                logger.info('Getting initial model grid')
                x_initial, y_initial = func(initial, atmtype=options['model'], driver='synth', version=options['MOOGv'], N=N, r=r, fout=fout, options=options)
                logger.info('Save initial synthetic spectrum')
                save_synth_spec(x_initial, y_initial, fname='initial.spec')
                logger.info('Interpolation successful.')
                logger.info('Setting solar values {0}, {1}, {2}, {3}'.format(*initial))

            elif len(line) == 6:
                logger.info('Initial parameters given by user.')
                initial = map(float, line[1:-1])
                initial[0] = int(initial[0])
                logger.info('Initial parameters: {0}, {1}, {2}, {3}'.format(*initial))
                options = _options(line[-1])

                #Create initial synthetic model
                N, r, fout = read_linelist(line[0])
                logger.info('Getting initial model grid')
                x_initial, y_initial = func(initial, atmtype=options['model'], driver='synth', version=options['MOOGv'], N=N, r=r, fout=fout, options=options)
                logger.info('Save initial synthetic spectrum')
                save_synth_spec(x_initial, y_initial, fname='initial.spec')
                logger.info('Interpolation successful.')
                logger.info('Setting solar values {0}, {1}, {2}, {3}'.format(*initial))

                if options['observations']:
                    # Check if observations exit, if not pass another line
                    if (not (os.path.isfile('spectra/%s' % options['observations'])) and (not os.path.isfile(options['observations']))):
                        logger.error('Error: %s not found.' % options['observations'])
                        continue

                    print('This is your observed spectrum: %s' % options['observations'])
                    x_obs, y_obs = read_obs_intervals('spectra/%s' % options['observations'], N, r)
                    if options['minimize']:
                        print('Starting minimization...')
                        logger.info('Starting the minimization procedure...')
                        params, x_final, y_final = minimize_synth(initial, x_obs, y_obs, N, r, fout, options)
                        logger.info('Minimization done.')

                        #Some statistics. Here synthetic spectrum is interpolated to the observed
                        chi2(x_obs, y_obs, x_initial, y_initial)
                        chi2(x_obs, y_obs, x_final, y_final)

                else:
                    x_obs, y_obs = (None, None)
                    x_final, y_final = (None, None)

                if options['plot']: #if there in no observed only the synthetic will be plotted
                    plot(x_obs, y_obs, x_initial, y_initial)
                    if options['minimize']:
                        plot(x_obs, y_obs, x_initial, y_initial)
                        plot(x_obs, y_obs, x_final, y_final)

            else:
                logger.error('Could not process information for this line list: %s' % line)
                continue

            if options['model'] != 'kurucz95' and options['model'] != 'apogee_kurucz' and options['model'] != 'marcs':
                logger.error('Your request for type: %s is not available' % model)
                continue

            # Options not in use will be removed
            if __name__ == '__main__':
                options['GUI'] = False  # Running batch mode
            else:
                options['GUI'] = True  # Running GUI mode
            options.pop('spt')

    return

if __name__ == '__main__':
    synthdriver()
