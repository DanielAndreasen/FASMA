#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import logging
import os
from shutil import copyfile
import yaml
import numpy as np
from utils import GetModels, _update_par
from interpolation import interpolator
from interpolation import save_model
from utils import fun_moog, Readmoog
from utils import error
from minimization import Minimize


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
        copyfile('summary.out', 'results/%s.out' % linelist)
    else:
        copyfile('summary.out', 'results/%s.NC.out' % linelist)


def _options(options=None):
    """Reads the options inside the config file"""
    defaults = {'spt': False,
                'weights': 'null',
                'model':'kurucz95',
                'fix_teff': False,
                'fix_logg': False,
                'fix_feh': False,
                'fix_vt': False,
                'refine': False,
                'iterations': 160,
                'EPcrit': 0.001,
                'RWcrit': 0.001,
                'ABdiffcrit': 0.01,
                'MOOGv': 2014,
                'loggLC': False,
                'outlier': False,
                'teffrange': False
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
        defaults['EPcrit'] = float(defaults['EPcrit'])
        defaults['RWcrit'] = float(defaults['RWcrit'])
        defaults['ABdiffcrit'] = float(defaults['ABdiffcrit'])
        defaults['MOOGv'] = int(defaults['MOOGv'])
        if defaults['outlier'] not in [False, '1Iter', '1Once', 'allIter', 'allOnce']:
            print('Invalid option set for option "outlier"')
            defaults['outlier'] = False
        return defaults


def _output(overwrite=None, header=None, parameters=None):
    """Create the output file 'results.csv'

    Input
    -----
    overwrite - Overwrite the file
    header    - Only use True if this is for the file to be created
    """
    if header:
        hdr = ['linelist', 'teff', 'tefferr', 'logg', 'loggerr', 'feh', 'feherr',
               'vt', 'vterr', 'convergence', 'fixteff', 'fixlogg', 'fixfeh', 'fixvt',
               'loggLC', 'outlier', 'weights', 'model', 'refine', 'EPcrit', 'RWcrit',
               'ABdiffcrit']
        if overwrite:
            with open('results.csv', 'w') as output:
                output.write('\t'.join(hdr)+'\n')
        else:
            if not os.path.isfile('results.csv'):
                with open('results.csv', 'w') as output:
                    output.write('\t'.join(hdr)+'\n')
    else:
        with open('results.csv', 'a') as output:
            output.write('\t'.join(map(str, parameters))+'\n')


def _setup(line):
    """Do the setup with initial parameters and options"""
    if len(line) == 1:
        initial = [5777, 4.44, 0.00, 1.00]
        options = _options()
    elif len(line) == 5:
        initial = map(float, line[1::])
        initial[0] = int(initial[0])
        options = _options()
    elif len(line) == 2:
        options = _options(line[1])
        if options['spt']:
            Teff, logg = _getSpt(options['spt'])
            mic = _getMic(Teff, logg)
            initial = (Teff, logg, 0.00, mic)
        else:
            initial = [5777, 4.44, 0.00, 1.00]
    elif len(line) == 6:
        initial = map(float, line[1:-1])
        initial[0] = int(initial[0])
        options = _options(line[-1])
    return initial, options


def _outlierRunner(type, linelist, parameters, options):
    """Remove the potential outliers based on a given type

    Input
    -----
    type - Can be '1Iter', '1Once', 'allIter', or 'allOnce'
    """
    tmpll = 'linelist/tmplinelist.moog'
    copyfile('linelist/'+linelist, tmpll)
    _update_par(line_list=tmpll)
    newLineList = False
    Noutlier = 0
    outliers = hasOutlier()
    if type == '1Iter':
        # Remove one outlier above 3 sigma iteratively
        while outliers:
            Noutlier += 1
            newLineList = True  # At the end, create a new linelist
            wavelength = outliers[max(outliers.keys())]
            removeOutlier(tmpll, wavelength)
            print('Removing line: %.2f. Outliers removed: %d' % (wavelength, Noutlier))
            print('Restarting the minimization routine...\n')
            function = Minimize(parameters, fun_moog, **options)
            parameters, converged = function.minimize()
            outliers = hasOutlier()

    elif type == '1Once':
        # Remove one outlier above 3 sigma once
        if outliers:
            Noutlier += 1
            newLineList = True  # At the end, create a new linelist
            wavelength = outliers[max(outliers.keys())]
            removeOutlier(tmpll, wavelength)
            print('Removing line: %.2f. Outliers removed: %d' % (wavelength, Noutlier))
            print('Restarting the minimization routine...')
            function = Minimize(parameters, fun_moog, **options)
            parameters, converged = function.minimize()
            outliers = hasOutlier()

    elif type == 'allIter':
        # Remove all outliers above 3 sigma iteratively
        while outliers:
            newLineList = True  # At the end, create a new linelist
            for wavelength in outliers.itervalues():
                removeOutlier(tmpll, wavelength)
                Noutlier += 1
                print('Removing line: %.2f. Outliers removed: %d' % (wavelength, Noutlier))
            print('Restarting the minimization routine...')
            function = Minimize(parameters, fun_moog, **options)
            parameters, converged = function.minimize()
            outliers = hasOutlier()

    elif type == 'allOnce':
        # Remove all outliers above 3 sigma once
        if outliers:
            newLineList = True  # At the end, create a new linelist
            for wavelength in outliers.itervalues():
                removeOutlier(tmpll, wavelength)
                Noutlier += 1
                print('Removing line: %.2f. Outliers removed: %d' % (wavelength, Noutlier))
            print('Restarting the minimization routine...')
            function = Minimize(parameters, fun_moog, **options)
            parameters, converged = function.minimize()
            outliers = hasOutlier()

    if newLineList:
        newName = linelist.replace('.moog', '_outlier.moog')
        copyfile(tmpll, 'linelist/'+newName)
        os.remove(tmpll)
        _update_par(line_list='linelist/'+newName)
        return newLineList, newName
    _update_par(line_list=linelist)
    os.remove(tmpll)
    return newLineList, linelist


def hasOutlier(MOOGv=2014):
    """Function that reads the summary.out file and return a dictionary
    with key being the deviation (above 3 sigma), and value the wavelength"""
    idx = 1 if MOOGv > 2013 else 0
    s = Readmoog(version=MOOGv)
    d = s.fe_statistics()
    fe1 = d[-2]  # All the FeI lines
    fe2 = d[-1]  # All the FeII lines
    m1, m2 = np.mean(fe1[:, 5+idx]), np.mean(fe2[:, 5+idx])
    s1, s2 = 3*np.std(fe1[:, 5+idx]), 3*np.std(fe2[:, 5+idx])

    d = {}
    for i, fe1i in enumerate(fe1[:, 5+idx]):
        dev = abs(fe1i-m1)
        if dev >= s1:
            d[dev] = fe1[i, 0]
    for i, fe2i in enumerate(fe2[:, 5+idx]):
        dev = abs(fe2i-m2)
        if dev >= s2:
            d[dev] = fe2[i, 0]

    if len(d.keys()):
        return d
    else:
        return False


def removeOutlier(fname, wavelength):
    """Remove an outlier from the linelist fname, and save it in the same name

    Input:
    fname      -- Name of the linelist
    wavelength -- The wavelength of the line to remove
    """
    wavelength = str(round(wavelength, 2))
    with open(fname, 'r') as lines:
        fout = ''
        for line in lines:
            if line.replace(' ','').startswith(wavelength):
                continue
            fout += line
    with open(fname, 'w') as f:
        f.writelines(fout)


def ewdriver(starLines='StarMe.cfg', overwrite=None):
    """The function that glues everything together

    Input:
    starLines   -   Configuration file (default: StarMe.cfg)
    overwrite   -   Overwrite the results.csv file (default: False)

    Output:
    <linelist>.(NC).out     -   NC=not converged.
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

    # Creating the output file
    _output(overwrite=overwrite, header=True)

    with open(starLines, 'r') as lines:
        for line in lines:
            if not line[0].isalnum():
                logger.debug('Skipping header: %s' % line.strip())
                continue
            logger.info('Line list: %s' % line.strip())
            line = line.strip()
            line = line.split(' ')
            if len(line) not in [1, 2, 5, 6]:
                logger.error('Could not process information for this line list: %s' % line)
                continue
            # Check if the linelist is inside the directory if not log it and pass to next linelist
            if not os.path.isfile('linelist/%s' % line[0]):
                logger.error('Error: linelist/%s not found.' % line[0])
                parameters = None
                continue
            else:
                _update_par(line_list='linelist/%s' % line[0])

            initial, options = _setup(line)
            logger.info('Initial parameters: {0}, {1}, {2}, {3}'.format(*initial))
            if options['fix_teff']:
                logger.info('Effective temperature fixed at: %i' % initial[0])
            elif options['fix_logg']:
                logger.info('Surface gravity fixed at: %s' % initial[1])
            elif options['fix_feh']:
                logger.info('Metallicity fixed at: %s' % initial[2])
            elif options['fix_vt']:
                logger.info('Micro turbulence fixed according to an emperical relation')

            # Setting the models to use
            if options['model'] not in ['kurucz95', 'apogee_kurucz']:
                logger.error('Your request for type: %s is not available' % model)
                continue

            # Get the initial grid models
            logger.info('Getting initial model grid')
            grid = GetModels(teff=initial[0], logg=initial[1], feh=initial[2], atmtype=options['model'])
            models, nt, nl, nf = grid.getmodels()
            logger.info('Initial interpolation of model...')
            inter_model = interpolator(models,
                                       teff=(initial[0], nt),
                                       logg=(initial[1], nl),
                                       feh=(initial[2], nf))
            save_model(inter_model, params=initial)
            logger.info('Interpolation successful.')

            # Adjusting the options for the minimization routine
            if __name__ == '__main__':
                options['GUI'] = False  # Running batch mode
            else:
                options['GUI'] = True  # Running GUI mode

            logger.info('Starting the minimization procedure...')
            function = Minimize(initial, fun_moog, **options)
            try:
                parameters, converged = function.minimize()
            except ValueError:
                print('No FeII lines were measured.')
                print('Skipping to next linelist..\n')
                logger.error('No FeII lines found for %s. Skipping to next linelist' % line[0])
                continue

            if options['outlier']:
                newLineList, newName = _outlierRunner(options['outlier'], line[0], parameters, options)
                line[0] = newName
            else:
                newLineList = False

            if options['teffrange']:
                d = np.loadtxt('rawLinelist/coolNormalDiff.lines')
                ll = np.loadtxt('linelist/%s' % line[0], skiprows=1, usecols=(0,))
                normalLL = np.in1d(ll, d)
                if np.any(normalLL) and (parameters[0] > 7000):
                    logger.warning('Effective temperature probably to high for this line list')
                elif np.any(normalLL) and (parameters[0] < 5200):
                    logger.info('Removing lines from the line list to compensate for the low Teff')
                    print('Removing lines to compensate for low Teff\n')
                    for li in ll[normalLL]:
                        removeOutlier('linelist/%s' % line[0], li)

                    # Restart the minimization procedure from the last best point
                    function = Minimize(parameters, fun_moog, **options)
                    try:
                        parameters, converged = function.minimize()
                    except ValueError:
                        print('No FeII lines were measured.')
                        print('Skipping to next linelist..\n')
                        logger.error('No FeII lines found for %s. Skipping to next linelist' % line[0])
                    if options['outlier']:
                        newLineList, newName = _outlierRunner(options['outlier'], line[0], parameters, options)
                        line[0] = newName
                    else:
                        newLineList = False

            # Refine the parameters
            if converged and options['refine']:
                logger.info('Refining the parameters')
                print('\nRefining the parameters')
                print('This might take some time...')
                options['EPcrit'] = 0.001
                options['RWcrit'] = 0.001
                options['ABdiffcrit'] = 0.01
                function = Minimize(parameters, fun_moog, **options)
                p1, converged = function.minimize()
                if converged:
                    print('Adjusting the final parameters')
                    parameters = p1  # overwrite with new best results
            logger.info('Finished minimization procedure')

            _renaming(line[0], converged)
            parameters = error(line[0], converged, atmtype=options['model'], version=options['MOOGv'], weights=options['weights'])
            parameters = list(parameters)
            if options['loggLC']:
                parameters[2] = round(parameters[2] - 4.57E-4*parameters[0] + 2.59, 2)

            tmp = [line[0]] + parameters +\
                  [converged, options['fix_teff'], options['fix_logg'],
                  options['fix_feh'], options['fix_vt'], options['loggLC'],
                  options['outlier']]+[options['weights'], options['model'],
                  options['refine'], options['EPcrit'], options['RWcrit'], options['ABdiffcrit']]
            _output(parameters=tmp)
            logger.info('Saved results to: results.csv')

            if __name__ == '__main__':
                if converged:
                    print('\nCongratulation, you have won! Your final parameters are:')
                else:
                    print('\nSorry, you did not win. However, your final parameters are:')
                try:
                    print(u' Teff:{:>8d}\u00B1{:d}\n logg:{:>8.2f}\u00B1{:1.2f}\n [Fe/H]:{:>+6.2f}\u00B1{:1.2f}\n vt:{:>10.2f}\u00B1{:1.2f}\n\n\n\n'.format(*parameters))
                except UnicodeEncodeError:
                    print(' Teff:{:>8d}({:d})\n logg:{:>8.2f}({:1.2f})\n [Fe/H]:{:>+6.2f}({:1.2f})\n vt:{:>10.2f}({:1.2f})\n\n\n\n'.format(*parameters))
            elif __name__ == 'ewDriver':
                if converged:
                    print('\nCongratulation, you have won! Your final parameters are:')
                else:
                    print('\nSorry, you did not win. However, your final parameters are:')
                print(u' Teff:{:>8d}+/-{:d}\n logg:{:>8.2f}+/-{:1.2f}\n [Fe/H]:{:>+6.2f}+/-{:1.2f}\n vt:{:>10.2f}+/-{:1.2f}\n\n\n\n'.format(*parameters))
    return parameters


if __name__ == '__main__':
    parameters = ewdriver()
