#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import logging
import os

from utils import GetModels, _update_par_synth
from interpolation import interpolator
from interpolation import save_model
from utils import read_observations, _run_moog, plot_synthetic
from utils import interpol_synthetic


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
                'plotpars': 0,
                'plot': 0,
                'step_wave': 0.01,
                'step_flux': 5.0,
                'model': 'kurucz95',
                'observations': 'observations',
                'resolution': 0.6,
                'vmac': 0.0,
                'vsini': 0.0,
                'limb': 0.0,
                'lorentz': 0.0,
                'teff': False,
                'logg': False,
                'feh': False,
                'vt': False,
                'iterations': 160,
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
        defaults['step_wave'] = float(defaults['step_wave'])
        defaults['step_flux'] = float(defaults['step_flux'])
        defaults['plotpars'] = int(defaults['plotpars'])
        defaults['plot'] = int(defaults['plot'])
        defaults['observations'] = str(defaults['observations'])
        defaults['resolution'] = float(defaults['resolution'])
        defaults['vmac'] = float(defaults['vmac'])
        defaults['vsini'] = float(defaults['vsini'])
        defaults['limb'] = float(defaults['limb'])
        defaults['lorentz'] = float(defaults['lorentz'])
        return defaults


def synthdriver(starLines='StarMe.cfg'):
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
            fix_teff = False
            fix_logg = False
            fix_feh = False
            fix_vt = False
            fix_vmacro = False
            fix_vsini = False
            line = line.strip()
            line = line.split(' ')

            # Check if the linelist is inside the directory if not log it and pass to next linelist
            if not os.path.isfile('linelist/%s' % line[0]):
                logger.error('Error: linelist/%s not found.' % line[0])
                parameters = None
                continue

            if len(line) == 7:
                logger.info('Initial parameters given by the user.')
                initial = map(float, line[1:5])
                options = _options()
                _update_par_synth(line_list='linelist/%s' % line[0], start_wave=line[5], end_wave=line[6])
                logger.info('Initial parameters: {0}, {1}, {2}, {3}'.format(*initial))

            elif len(line) == 3:
                initial = (5777, 4.44, 0.00, 1.00)
                options = _options()
                _update_par_synth(line_list='linelist/%s' % line[0], start_wave=line[5], end_wave=line[6])
                logger.info('Setting solar values {0}, {1}, {2}, {3}'.format(*initial))

            elif len(line) == 4:
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

            elif len(line) == 8:
                logger.info('Initial parameters given by user and some parameters fixed.')
                initial = map(float, line[1:5])
                options = _options(line[-1])
                if options['observations']:
                    _update_par_synth(line_list='linelist/%s' % line[0], start_wave=line[5],
                                      end_wave=line[6], step_wave=options['step_wave'],
                                      step_flux=options['step_flux'], vsini=options['vsini'],
                                      plot=options['plot'], plotpars=options['plotpars'],
                                      vmac=options['vmac'], resolution=options['resolution'],
                                      lorentz=options['lorentz'], limb=options['limb'], obfile=options['observations'])
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
                else:
                    _update_par_synth(line_list='linelist/%s' % line[0], start_wave=line[5],
                                      end_wave=line[6], step_wave=options['step_wave'],
                                      step_flux=options['step_flux'], vsini=options['vsini'],
                                      plot=options['plot'], plotpars=options['plotpars'],
                                      vmac=options['vmac'], resolution=options['resolution'],
                                      lorentz=options['lorentz'], limb=options['limb'])
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
            # if initial[1] > 4.99:  # quick fix
            #     initial[1] = 4.99
            # models, nt, nl, nf = _get_model(teff=initial[0], logg=initial[1], feh=initial[2], atmtype=options['model'])
            # logger.info('Initial interpolation of model...')
            # inter_model = interpolator(models,
            #                            teff=(initial[0], nt),
            #                            logg=(initial[1], nl),
            #                            feh=(initial[2], nf))
            # save_model(inter_model, params=initial)
            # logger.info('Interpolation successful.')
    return 0


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import seaborn
    import numpy as np
    synthdriver()
    _run_moog(driver='synth')
    # plot_synthetic()
    wavelength_obs, flux_obs = np.loadtxt('sun_harps_ganymede.txt', unpack=True, usecols=(0, 1))
    wavelength_obs, flux_obs, flux_inter_synth = interpol_synthetic(wavelength_obs, flux_obs, 6444.672, 6447.340)

    flux_obs /= np.median(flux_obs)
    # Normalization (use first 50 points below 1.2 as constant continuum)
    maxes = flux_obs[(flux_obs < 1.2)].argsort()[-50:][::-1]
    flux_obs /= np.median(flux_obs[maxes])

    # flux_obs *= np.median(flux_inter_synth)/np.median(flux_obs)

    plt.plot(wavelength_obs, flux_obs, '-k', wavelength_obs, flux_inter_synth, '-r')
    plt.plot(wavelength_obs, flux_obs - flux_inter_synth + 1, '--g')
    plt.show()
    # read_observations('sun_6000-7000.txt', 6444.0, 6448.0)
