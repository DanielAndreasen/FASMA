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
        cmd = 'cp summary.out results/%s.out' % linelist
    else:
        cmd = 'cp summary.out results/%s.NC.out' % linelist
    os.system(cmd)


def _options(options=False):
    '''Reads the options inside the config file'''
    defaults = {'spt': False,
                'weights': 'null',
                'model':'kurucz95',
                'teff': False,
                'logg': False,
                'feh': False,
                'vt': False,
                'refine': False,
                'iterations': 160,
                'EPcrit': 0.001,
                'RWcrit': 0.001,
                'ABdiffcrit': 0.01,
                'MOOGv': 2014,
                'loggLC': False,
                'outlier': False
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
                # Clever way to change the boolean
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
            # print(wavelength, line.strip())
            if line.replace(' ','').startswith(wavelength):
                continue
            fout += line
    with open(fname, 'w') as f:
        f.writelines(fout)


def ewdriver(starLines='StarMe.cfg', overwrite=False):
    """The function that glues everything together

    Input:
    starLines   -   Configuration file (default: StarMe.cfg)
    overwrite   -   Overwrite the results.csv file (default: False)

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
    
    if overwrite:
        with open('results.csv', 'w') as output:
            tmp = ['linelist', 'teff', 'tefferr', 'logg', 'loggerr', 'feh', 'feherr', 'vt', 'vterr', 'convergence',  'fixteff',  'fixlogg', 'fixfeh',  'fixvt','loggLC','outlier', 'weights',  'model',  'refine', 'EPcrit', 'RWcrit',  'ABdiffcrit'] 
            output.write('\t'.join(tmp)+'\n')
    else:
        if not os.path.isfile('results.csv'):
            with open('results.csv', 'w') as output:
                tmp = ['linelist', 'teff', 'tefferr', 'logg', 'loggerr', 'feh', 'feherr', 'vt', 'vterr', 'convergence',  'fixteff',  'fixlogg', 'fixfeh',  'fixvt','loggLC','outlier', 'weights',  'model',  'refine', 'EPcrit', 'RWcrit',  'ABdiffcrit'] 
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
            if not line[0].isalnum():
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
                initial = [5777, 4.44, 0.00, 1.00]
                options = _options()
                logger.info('Setting solar values {0}, {1}, {2}, {3}'.format(*initial))

            elif len(line) == 5:
                logger.info('Initial parameters given by the user.')
                initial = map(float, line[1::])
                initial[0] = int(initial[0])
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
                    initial = [5777, 4.44, 0.00, 1.00]
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
                initial[0] = int(initial[0])
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
            if options['model'] != 'kurucz95' and options['model'] != 'apogee_kurucz':
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
            refine = options.pop('refine')
            # Fixing parameters
            fix_teff = options.pop('teff')
            fix_logg = options.pop('logg')
            fix_feh = options.pop('feh')
            fix_vt = options.pop('vt')
            loggLC = options.pop('loggLC')
            outlier = options.pop('outlier')

            fff = Minimize(initial, fun_moog,
                           fix_teff=fix_teff, fix_logg=fix_logg,
                           fix_feh=fix_feh, fix_vt=fix_vt, **options)
            parameters, converged = fff.minimize()

            newLineList = False
            if outlier:
                tmpll = 'linelist/tmplinelist.moog'
                Noutlier = 0
                outliers = hasOutlier()
                if outlier == '1Iter':
                    # Remove one outlier above 3 sigma iteratively
                    while outliers:
                        Noutlier += 1
                        newLineList = True  # At the end, create a new linelist
                        if not os.path.isfile(tmpll):
                            copyfile('linelist/'+line[0], tmpll)
                            _update_par(line_list=tmpll)
                        wavelength = outliers[max(outliers.keys())]
                        removeOutlier(tmpll, wavelength)
                        print('Removing line: %.2f. Outliers removed: %d' % (wavelength, Noutlier))
                        print('Restarting the minimization routine...')
                        fff = Minimize(parameters, fun_moog,
                                       fix_teff=fix_teff, fix_logg=fix_logg,
                                       fix_feh=fix_feh, fix_vt=fix_vt, **options)
                        parameters, converged = fff.minimize()
                        outliers = hasOutlier()

                if outlier == '1Once':
                    # Remove one outlier above 3 sigma once
                    if outliers:
                        Noutlier += 1
                        newLineList = True  # At the end, create a new linelist
                        if not os.path.isfile(tmpll):
                            copyfile('linelist/'+line[0], tmpll)
                            _update_par(line_list=tmpll)
                        wavelength = outliers[max(outliers.keys())]
                        removeOutlier(tmpll, wavelength)
                        print('Removing line: %.2f. Outliers removed: %d' % (wavelength, Noutlier))
                        print('Restarting the minimization routine...')
                        fff = Minimize(parameters, fun_moog,
                                       fix_teff=fix_teff, fix_logg=fix_logg,
                                       fix_feh=fix_feh, fix_vt=fix_vt, **options)
                        parameters, converged = fff.minimize()
                        outliers = hasOutlier()

                if outlier == 'allIter':
                    # Remove all outliers above 3 sigma iteratively
                    while outliers:
                        newLineList = True  # At the end, create a new linelist
                        if not os.path.isfile(tmpll):
                            copyfile('linelist/'+line[0], tmpll)
                            _update_par(line_list=tmpll)
                        for wavelength in outliers.itervalues():
                            removeOutlier(tmpll, wavelength)
                            Noutlier += 1
                            print('Removing line: %.2f. Outliers removed: %d' % (wavelength, Noutlier))
                        print('Restarting the minimization routine...')
                        fff = Minimize(parameters, fun_moog,
                                       fix_teff=fix_teff, fix_logg=fix_logg,
                                       fix_feh=fix_feh, fix_vt=fix_vt, **options)
                        parameters, converged = fff.minimize()
                        outliers = hasOutlier()

                if outlier == 'allOnce':
                    # Remove all outliers above 3 sigma once
                    if outliers:
                        newLineList = True  # At the end, create a new linelist
                        if not os.path.isfile(tmpll):
                            copyfile('linelist/'+line[0], tmpll)
                            _update_par(line_list=tmpll)
                        for wavelength in outliers.itervalues():
                            removeOutlier(tmpll, wavelength)
                            Noutlier += 1
                            print('Removing line: %.2f. Outliers removed: %d' % (wavelength, Noutlier))
                        print('Restarting the minimization routine...')
                        fff = Minimize(parameters, fun_moog,
                                       fix_teff=fix_teff, fix_logg=fix_logg,
                                       fix_feh=fix_feh, fix_vt=fix_vt, **options)
                        parameters, converged = fff.minimize()
                        outliers = hasOutlier()

                if newLineList:
                    newName = line[0].replace('.moog', '_outlier.moog')
                    copyfile(tmpll, 'linelist/'+newName)
                    _update_par(line_list='linelist/'+newName)

                if os.path.isfile(tmpll):
                    os.remove(tmpll)

            if converged and refine:
                logger.info('Refining the parameters')
                print('\nRefining the parameters')
                print('This might take some time...')
                options['EPcrit'] = 0.001
                options['RWcrit'] = 0.001
                options['ABdiffcrit'] = 0.01
                fff = Minimize(parameters, fun_moog,
                               fix_teff=fix_teff, fix_logg=fix_logg,
                               fix_feh=fix_feh, fix_vt=fix_vt, **options)
                p1, converged = fff.minimize()
                if converged:
                    print('reseting the parameters')
                    parameters = p1  # overwrite with new best results
            logger.info('Finished minimization procedure')

            if newLineList:
                _renaming(newName, converged)
                parameters = error(newName, converged, atmtype=options['model'], version=options['MOOGv'], weights=options['weights'])
            else:
                _renaming(line[0], converged)
                parameters = error(line[0], converged, atmtype=options['model'], version=options['MOOGv'], weights=options['weights'])
            parameters = list(parameters)
            if loggLC:
                parameters[2] = round(parameters[2] - 4.57E-4*parameters[0] + 2.59, 2)

            _ = options.pop('MOOGv')
            _ = options.pop('iterations')
            with open('results.csv', 'a') as output:
                if newLineList:
                    tmp = [newName] + parameters + [converged, fix_teff, fix_logg,fix_feh,fix_vt,loggLC,outlier] +[options['weights'], options['model'], refine, options['EPcrit'], options['RWcrit'], options['ABdiffcrit']]
                else:
                    tmp = [line[0]] + parameters +  [converged, fix_teff, fix_logg,fix_feh,fix_vt,loggLC,outlier] + [options['weights'], options['model'], refine, options['EPcrit'], options['RWcrit'], options['ABdiffcrit']]
                output.write('\t'.join(map(str, tmp))+'\n')
            logger.info('Saved results to: results.csv')

            if __name__ == '__main__':
                if converged:
                    print('\nCongratulation, you have won! Your final parameters are:')
                else:
                    print('\nSorry, you did not win. However, your final parameters are:')
                try:
                    print(u' Teff:    %i\u00B1%i\n logg:    %.2f\u00B1%.2f\n [Fe/H]: %.2f\u00B1%.2f\n vt:      %.2f\u00B1%.2f\n\n\n\n' %
                        (parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5],parameters[6],parameters[7]))
                except UnicodeEncodeError:
                    print('Teff:      %i(%i)\nlogg:    %.2f(%.2f)\n[Fe/H]:  %.2f(%.2f)\nvt:        %.2f(%.2f)\n\n\n\n' %
                         (parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5],parameters[6],parameters[7]))
            elif __name__ == 'ewDriver':
                if converged:
                    print('\nCongratulation, you have won! Your final parameters are:')
                else:
                    print('\nSorry, you did not win. However, your final parameters are:')
                print('Teff:      %i+/-%i\nlogg:    %.2f+/-%.2f\n[Fe/H]:  %.2f+/-%.2f\nvt:        %.2f+/-%.2f\n\n\n\n' %
                     (parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],parameters[5],parameters[6],parameters[7]))
    return parameters


if __name__ == '__main__':
    parameters = ewdriver()
