#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import logging
import os
import yaml
from utils import fun_moog_synth as func
from observations import read_obs_intervals, plot
from minimization import minimize_synth
from synthetic import read_linelist, save_synth_spec


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


def _getMac(teff, logg):
    """Calculate macro turbulence (Doyle et al. 2014)."""
    # 5200 < teff < 6400
    # 4.0 < logg < 4.6
    mac = 3.21 + (2.33 * (teff - 5777.) * (10**(-3)))
    + (2.00 * ((teff - 5777.)**2) * (10**(-6))) - (2.00 * (logg - 4.44))
    return round(mac, 2)


def _options(options=None):
    '''Reads the options inside the config file'''
    defaults = {'spt': False,
                'model': 'kurucz95',
                'MOOGv': 2014,
                'plotpars': 0,
                'fix_teff': False,
                'fix_logg': False,
                'fix_feh': False,
                'fix_vt': False,
                'fix_vmac': False,
                'fix_vsini': False,
                'plot': False,  # This is irrelevant with the batch.par value
                'damping': 1,
                'step_wave': 0.01,
                'step_flux': 10.0,
                'minimize': False,
                'observations': False,
                'snr': 100.0,
                'resolution': None,
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
                if option in ['teff', 'logg', 'feh', 'vt', 'vmac', 'vsini']:
                    option = 'fix_%s' % option
                defaults[option] = False if defaults[option] else True
        defaults['model'] = defaults['model'].lower()
        defaults['step_wave'] = float(defaults['step_wave'])
        defaults['step_flux'] = float(defaults['step_flux'])
        defaults['snr'] = float(defaults['snr'])
        defaults['plotpars'] = int(defaults['plotpars'])
        defaults['limb'] = float(defaults['limb'])
        defaults['lorentz'] = float(defaults['lorentz'])
        defaults['MOOGv'] = int(defaults['MOOGv'])
        return defaults


def wave_step(dl_obs, step_wave=0.01):
    '''Find the step of synthesis in wavelength depending the observations'''

    if dl_obs < step_wave:
        step_wave = dl_obs
    elif dl_obs > step_wave:
        step_wave = dl_obs
    else:
        step_wave
    return round(step_wave,3)


def synthdriver(starLines='StarMe_synth.cfg', overwrite=False):
    """The function that glues everything together

    Input
    -----
    starLines   -   Configuration file (default: StarMe.cfg)

    Output
    -----
    to be decided
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

            # Check if configuration parameters are correct
            if len(line) not in [1, 2, 7, 8]:
                logger.error('Could not process this information: %s' % line)
                continue

            # Check if the linelist is inside the directory, if not log it and pass to next linelist
            if not os.path.isfile('rawLinelist/%s' % line[0]):
                logger.error('Error: rawLinelist/%s not found.' % line[0])
                raise IOError('Linelist file does not exist in the folder!')

            ranges, atomic_data = read_linelist(line[0])
            # Create synthetic spectrum with solar values and the default options
            if len(line) == 1:
                initial = [5777, 4.44, 0.00, 1.00, 3.21, 1.90]
                options = _options()
                x_obs, y_obs = (None, None)
                # Create initial synthetic model
                logger.info('Getting initial model grid')
                x_initial, y_initial = func(initial, atmtype=options['model'],
                ranges=ranges, driver='synth', version=options['MOOGv'], **options)
                logger.info('Interpolation successful.')
                logger.info('Setting solar values {0}, {1}, {2}, {3}'.format(*initial))
                save_synth_spec(x_initial, y_initial, fname='initial.spec')
                logger.info('Save initial synthetic spectrum')
                print('Synthetic spectrum contains %s points' % len(x_initial))

            # Create synthetic spectrum with solar values and the options defined by the user
            elif len(line) == 2:
                options = _options(line[1])
                if options['spt']:
                    logger.info('Spectral type given: %s' % line[1])
                    Teff, logg = _getSpt(options['spt'])
                    mic = _getMic(Teff, logg)
                    mac = _getMac(Teff, logg)
                    initial = (Teff, logg, 0.00, mic, mac, 1.90)
                else:
                    initial = [5777, 4.44, 0.00, 1.00, 3.21, 1.90]

                # Create initial synthetic model
                logger.info('Getting initial model grid')
                x_initial, y_initial = func(initial, atmtype=options['model'],
                ranges=ranges, driver='synth', version=options['MOOGv'], **options)
                logger.info('Interpolation successful.')
                logger.info('Setting solar values {0}, {1}, {2}, {3}'.format(*initial))
                save_synth_spec(x_initial, y_initial, fname='initial.spec')
                logger.info('Save initial synthetic spectrum')
                print('Synthetic spectrum contains %s points' % len(x_initial))

                if options['observations']:
                    # Check if observations exit, if not pass another line
                    if (not (os.path.isfile('spectra/%s' % options['observations'])) and (not os.path.isfile(options['observations']))):
                        logger.error('Error: %s not found.' % options['observations'])
                        continue
                    print('This is your observed spectrum: %s' % options['observations'])
                    x_obs, y_obs = read_obs_intervals('spectra/%s' % options['observations'], ranges, snr=options['snr'])
                    dl_obs = x_obs[1] - x_obs[0]
                    print('Observed spectrum contains %s points' % len(x_obs))
                    # If observations exists, then why not, find best fit?
                    if options['minimize']:
                        print('Starting minimization...')
                        logger.info('Starting the minimization procedure...')
                        options['step_wave'] = wave_step(dl_obs)
                        params, x_final, y_final = minimize_synth(initial, x_obs, y_obs, ranges=ranges, **options)
                        logger.info('Minimization done.')

                else:
                    x_obs, y_obs = (None, None)
                    x_final, y_final = (None, None)

                if options['plot']:  # if there in no observed only the synthetic will be plotted
                    plot(x_obs, y_obs, x_initial, y_initial)
                    if options['minimize']:
                        # Plot also final spectra
                        plot(x_obs, y_obs, x_final, y_final)

            # Create spectra with parameters defined by the user but the rest are set to default
            elif len(line) == 7:
                logger.info('Initial parameters given by the user.')
                initial = map(float, line[1::])
                initial[0] = int(initial[0])
                options = _options()
                x_obs, y_obs = (None, None)  # No observed spectrum, no plots

                # Create initial synthetic model
                logger.info('Getting initial model grid')
                x_initial, y_initial = func(initial, atmtype=options['model'],
                ranges=ranges, driver='synth', version=options['MOOGv'], **options)
                logger.info('Save initial synthetic spectrum')
                save_synth_spec(x_initial, y_initial, fname='initial.spec')
                logger.info('Interpolation successful.')
                logger.info('Setting solar values {0}, {1}, {2}, {3}'.format(*initial))
                print('Synthetic spectrum contains %s points' % len(x_initial))

            # Create synthetic spectra with values set by the user and options altered.
            elif len(line) == 8:
                logger.info('Initial parameters given by user.')
                initial = map(float, line[1:-1])
                initial[0] = int(initial[0])
                logger.info('Initial parameters: {0}, {1}, {2}, {3}'.format(*initial))
                options = _options(line[-1])

                # Create initial synthetic model
                logger.info('Getting initial model grid')
                x_initial, y_initial = func(initial, atmtype=options['model'],
                ranges=ranges, driver='synth', version=options['MOOGv'], **options)
                logger.info('Save initial synthetic spectrum')
                save_synth_spec(x_initial, y_initial, fname='initial.spec')
                logger.info('Interpolation successful.')
                logger.info('Setting solar values {0}, {1}, {2}, {3}'.format(*initial))
                print('Synthetic spectrum contains %s points' % len(x_initial))

                if options['observations']:
                    # Check if observations exit, if not pass another line
                    if (not (os.path.isfile('spectra/%s' % options['observations'])) and (not os.path.isfile(options['observations']))):
                        logger.error('Error: %s not found.' % options['observations'])
                        continue
                    print('This is your observed spectrum: %s' % options['observations'])
                    x_obs, y_obs = read_obs_intervals('spectra/%s' % options['observations'], ranges, snr=options['snr'])
                    dl_obs = x_obs[1] - x_obs[0]
                    print('Observed spectrum contains %s points' % len(x_obs))

                    # If observations exists, then why not, find best fit?
                    if options['minimize']:
                        print('Starting minimization...')
                        logger.info('Starting the minimization procedure...')
                        options['step_wave'] = wave_step(dl_obs)
                        params, x_final, y_final = minimize_synth(initial, x_obs, y_obs, ranges=ranges, **options)
                        logger.info('Minimization done.')

                else:
                    x_obs, y_obs = (None, None)
                    x_final, y_final = (None, None)

                if options['plot']:  # if there in no observed only the synthetic will be plotted
                    plot(x_obs, y_obs, x_initial, y_initial)
                    if options['minimize']:
                        # Plot also final spectra
                        plot(x_obs, y_obs, x_final, y_final)

            else:
                logger.error('Could not process information for this line list: %s' % line)
                continue

            if options['model'] != 'kurucz95' and options['model'] != 'apogee_kurucz' and options['model'] != 'marcs':
                logger.error('Your request for type: %s is not available' % options['model'])
                continue

            # Options not in use will be removed
            if __name__ == '__main__':
                options['GUI'] = False  # Running batch mode
            else:
                options['GUI'] = True  # Running GUI mode
            options.pop('spt')

    return

if __name__ == '__main__':
    import time
    start_time = time.time()
    synthdriver()
    end_time = time.time()-start_time
    print('Calculations finished in %s sec' % int(end_time))
