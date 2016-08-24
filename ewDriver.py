#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import logging
import os
from shutil import copyfile
import yaml
import numpy as np
from utils import _update_par
from interpolation import interpolator
from utils import fun_moog, Readmoog
from utils import error
from minimization import Minimize


def _getSpt(spt):
    """Get the spectral type from a string like 'F5V'.

    Input
    -----
    spt : str
      The spectral type given in the form: F5V

    Output
    ------
    teff : int
      The effective temperature corresponding to the spectral type
    logg : float
      The surface gravity corresponding to the spectral type
    """
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
    """Calculate micro turbulence based on emperical relations.
    For dwarfs (logg>=3.95) we use the relation by Tsantaki 2013,
    and for giants (logg<3.95) we use the relation by Adibekyan 2015

    Input
    -----
    teff : int
      The effective temperature
    logg : float
      The surface gravity
    feh : float
      The metallicity ([Fe/H] notation)

    Output
    ------
    vt : float
      The microturbulence
    """
    if logg >= 3.95:  # Dwarfs Tsantaki 2013
        mic = 6.932 * teff * (10**(-4)) - 0.348 * logg - 1.437
        return round(mic, 2)
    else:  # Giants Adibekyan 2015
        mic = 2.72 - (0.457 * logg) + (0.072 * feh)
        return round(mic, 2)


def _tmcalc(linelist):
    """Initial guess on atmospheric parameters. Estimate based on TMCalc

    Input
    -----
    linelist : str
      The raw output from ARES

    Output
    ------
    teff : int
      The effective temperature
    logg : float
      Solar surface gravity, 4.44. TMCalc can not give an estimate on this
    feh : float
      The metallicity ([Fe/H] notation)
    vt : float
      The microturbulence calculated from the emperical relation in _getMic
    """
    import sys
    sys.path.append('TMCALC/tmcalc_cython')
    from tmcalc_module import get_temperature_py as get_teff
    from tmcalc_module import get_feh_py as get_feh

    data = np.loadtxt('linelist/%s' % linelist, skiprows=1, usecols=(0, 4))
    X = np.zeros((data.shape[0], 9))
    X[:, 0] = data[:, 0]
    X[:, 4] = data[:, 1]
    np.savetxt('tmp.ares', X, '%.2f')

    teff = get_teff('TMCALC/tmcalc_cython/gteixeira_teff_cal.dat', 'tmp.ares')
    feh = get_feh('TMCALC/tmcalc_cython/gteixeira_feh_cal.dat', 'tmp.ares', teff[0], teff[1], teff[2], teff[3])[0]
    teff = teff[0]
    vt = _getMic(teff, 4.44, feh)
    os.remove('tmp.ares')
    return [teff, 4.44, feh, vt]


def _renaming(linelist, converged):
    """Save the output in a file related to the linelist.

    Input
    -----
    linelist : str
      The PATH for the line list
    converged : bool
      True if the minimization converged, False otherwise

    Output
    ------
    results/<linelist>(.NC).out : file
      Copy the final summary.out to this file
    """
    if converged:
        copyfile('summary.out', 'results/%s.out' % linelist)
    else:
        copyfile('summary.out', 'results/%s.NC.out' % linelist)


def _options(options=None):
    """Reads the options inside the config file.

    Input
    -----
    options : str (optional)
      The line from the configuration file with the user options

    Output
    ------
    defaults : dict
      A dictionary with all the options for the EW method.
    """
    defaults = {'spt': False,
                'weights': 'null',
                'model': 'kurucz95',
                'fix_teff': False,
                'fix_logg': False,
                'fix_feh': False,
                'fix_vt': False,
                'refine': False,
                'iterations': 160,
                'EPcrit': 0.001,
                'RWcrit': 0.003,
                'ABdiffcrit': 0.01,
                'MOOGv': 2014,
                'outlier': False,
                'teffrange': False,
                'autofixvt': False,
                'tmcalc': False,
                'sigma': 3
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


def _output(linelist=None, parameters=None, converged=None, options=None,
            overwrite=None, header=None):
    """Create the output file 'results.csv'

    Input
    -----
    linelist : str
      The line list for the star
    parameters : list
      The parameters for the star to be saved
    converged : bool
      True if the minimization routine converged, False otherwise
    options : dict
      The options dictionary
    overwrite : bool
      Overwrite the results.csv file
    header : bool
      Only use True if this is for the file to be created

    Output
    ------
    results.csv : file
      If overwrite is True, then make a new file, otherwise append the results
      to this file
    """
    if header:
        hdr = ['linelist', 'teff', 'tefferr', 'logg', 'loggerr', 'feh', 'feherr',
               'vt', 'vterr', 'loggastero', 'dloggastero', 'loggLC', 'dloggLC',
               'convergence', 'fixteff', 'fixlogg', 'fixfeh', 'fixvt', 'outlier',
               'weights', 'model', 'refine', 'EPcrit', 'RWcrit', 'ABdiffcrit']
        if overwrite:
            with open('results.csv', 'w') as output:
                output.write('\t'.join(hdr)+'\n')
        else:
            if not os.path.isfile('results.csv'):
                with open('results.csv', 'w') as output:
                    output.write('\t'.join(hdr)+'\n')
    else:
        tmp = [linelist] + parameters + [converged, options['fix_teff'],
                                         options['fix_logg'], options['fix_feh'],
                                         options['fix_vt'], options['outlier']] +\
            [options['weights'], options['model'], options['refine'],
             options['EPcrit'], options['RWcrit'], options['ABdiffcrit']]
        with open('results.csv', 'a') as output:
            output.write('\t'.join(map(str, tmp))+'\n')


def _setup(line):
    """Do the setup with initial parameters and options.

    Input
    -----
    line : list
      A line from the configuration file after being split at spaces

    Output
    ------
    initial : list
      The initial parameters for a given star
    options : dict
      The options to use for a given star
    """
    if len(line) == 1:
        initial = [5777, 4.44, 0.00, 1.00]
        options = _options()
    elif len(line) == 5:
        initial = map(float, line[1::])
        initial[0] = int(initial[0])
        options = _options()
    elif len(line) == 2:
        options = _options(line[1])
        initial = [5777, 4.44, 0.00, 1.00]
        if options['spt']:
            Teff, logg = _getSpt(options['spt'])
            mic = _getMic(Teff, logg, 0.00)
            initial = (Teff, logg, 0.00, mic)
        if options['tmcalc']:
            initial = _tmcalc(line[0])
    elif len(line) == 6:
        initial = map(float, line[1:-1])
        initial[0] = int(initial[0])
        options = _options(line[-1])
    return initial, options


def _outlierRunner(type, linelist, parameters, options):
    """Remove the potential outliers based on a given method. After outliers
    are removed, then restarts the minimization routine at the previous best
    found parameters.

    Input
    -----
    type : str
      The method to remove outliers (above n sigma).
       '1Iter': Remove 1 outlier iteratively.
       '1Once': Remove 1 outlier once.
       'allIter': Remove all outliers iteratively.
       'allOnce': Remove all outliers once.
    linelist : str
      The name of the line list
    parameters : list
      The parameters (used as a starting point for the minimization when it
      restarts)
    options : dict
      The options dictionary

    Output
    ------
    linelist : str
      The new name of the line list contains a _outlier.moog ending. Only
      applies if the were removed outliers.
    parameters : list
      The new parameters after outlier removal
    """
    tmpll = 'linelist/tmplinelist.moog'
    copyfile('linelist/'+linelist, tmpll)
    _update_par(line_list=tmpll)
    newLineList = False
    Noutlier = 0
    outliers = hasOutlier(MOOGv=options['MOOGv'], n=options['sigma'])
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
            outliers = hasOutlier(MOOGv=options['MOOGv'], n=options['sigma'])

    elif type == '1Once':
        # Remove one outlier above 3 sigma once
        if outliers:
            Noutlier += 1
            newLineList = True  # At the end, create a new linelist
            wavelength = outliers[max(outliers.keys())]
            removeOutlier(tmpll, wavelength)
            print('Removing line: %.2f. Outliers removed: %d' % (wavelength, Noutlier))
            print('Restarting the minimization routine...\n')
            function = Minimize(parameters, fun_moog, **options)
            parameters, converged = function.minimize()
            outliers = hasOutlier(MOOGv=options['MOOGv'], n=options['sigma'])

    elif type == 'allIter':
        # Remove all outliers above 3 sigma iteratively
        while outliers:
            newLineList = True  # At the end, create a new linelist
            for wavelength in outliers.itervalues():
                removeOutlier(tmpll, wavelength)
                Noutlier += 1
                print('Removing line: %.2f. Outliers removed: %d' % (wavelength, Noutlier))
            print('Restarting the minimization routine...\n')
            function = Minimize(parameters, fun_moog, **options)
            parameters, converged = function.minimize()
            outliers = hasOutlier(MOOGv=options['MOOGv'], n=options['sigma'])

    elif type == 'allOnce':
        # Remove all outliers above 3 sigma once
        if outliers:
            newLineList = True  # At the end, create a new linelist
            for wavelength in outliers.itervalues():
                removeOutlier(tmpll, wavelength)
                Noutlier += 1
                print('Removing line: %.2f. Outliers removed: %d' % (wavelength, Noutlier))
            print('Restarting the minimization routine...\n')
            function = Minimize(parameters, fun_moog, **options)
            parameters, converged = function.minimize()
            outliers = hasOutlier(MOOGv=options['MOOGv'], n=options['sigma'])

    if newLineList:
        newName = linelist.replace('.moog', '_outlier.moog')
        copyfile(tmpll, 'linelist/'+newName)
        os.remove(tmpll)
        _update_par(line_list='linelist/'+newName)
        return newName, parameters
    _update_par(line_list='linelist/'+linelist)
    os.remove(tmpll)
    return linelist, parameters


def hasOutlier(MOOGv=2014, n=3):
    """Function that reads the summary.out file and return a dictionary
    with key being the deviation (above n sigma), and value the wavelength.

    Input
    -----
    MOOGv : int
      The version of MOOG (default: 2014)
    n : float
      The number of sigma for the outlier identification (default: 3)

    Output
    ------
    d : dict
      A dictionary with {n*deviation: wavelength}
    """
    idx = 1 if MOOGv > 2013 else 0
    s = Readmoog(version=MOOGv)
    d = s.fe_statistics()
    fe1 = d[-2]  # All the FeI lines
    fe2 = d[-1]  # All the FeII lines
    m1, m2 = np.mean(fe1[:, 5+idx]), np.mean(fe2[:, 5+idx])
    s1, s2 = n*np.std(fe1[:, 5+idx]), n*np.std(fe2[:, 5+idx])

    d = {}
    for i, fe1i in enumerate(fe1[:, 5+idx]):
        dev = abs(fe1i-m1)
        if dev >= s1:
            d[dev] = fe1[i, 0]
    if fe2.shape[0] > 10:
        for i, fe2i in enumerate(fe2[:, 5+idx]):
            dev = abs(fe2i-m2)
            if dev >= s2:
                d[dev] = fe2[i, 0]

    if len(d.keys()):
        return d
    else:
        return False


def removeOutlier(fname, wavelength):
    """Remove an outlier from a line list, and save it in the same name

    Input
    -----
    fname : str
      Name of the line list
    wavelength : float
      The wavelength of the line to remove

    Output
    ------
    fname : file
      Remove the line from the line list and save it in the same name
    """
    wavelength = str(round(wavelength, 2))
    with open(fname, 'r') as lines:
        fout = ''
        for line in lines:
            if line.replace(' ', '').startswith(wavelength):
                continue
            fout += line
    with open(fname, 'w') as f:
        f.writelines(fout)


def genStar(starLines):
    """A generator for the configuration file.

    Input
    -----
    starLines : str
      The name of the configuration file

    Output
    ------
    initial : list
      The initial parameters for the minimization
    options : dict
      The options dictionary
    """

    lines = open(starLines, 'r')
    for line in lines:
        if not line[0].isalnum():
            # Header
            continue
        line = line.strip()
        line = line.split(' ')
        if len(line) not in [1, 2, 5, 6]:
            # Not the expected format
            continue
        initial, options = _setup(line)
        yield initial, options, line


def _prepare(linelist, initial, options):
    """Prepare the run with setup and first interpolation.

    Input
    -----
    linelist : str
      The line list to be used
    initial : list
      The initial parameters for the minimization
    options : dict
      The options dictionary

    Output
    ------
    options : dict
      The updated options dictionary with the GUI keyword

    If the line list does not exists the return None
    """

    if not os.path.isfile('linelist/%s' % linelist):
        return None
    else:
        _update_par(line_list='linelist/%s' % linelist)

    # Make the initial interpolation
    interpolator(params=initial, atmtype=options['model'])

    # Adjusting the options for the minimization routine
    if __name__ == '__main__':
        options['GUI'] = False  # Running batch mode
    else:
        options['GUI'] = True  # Running GUI mode

    return options


def _minizationRunner(initial, options):
    """A function to run the minimization routine

    Input
    -----
    initial : list
      The initial parameters for the minimization
    options : dict
      The options dictionary

    Output
    ------
    parameters : list
      The best found parameters
    converged : bool
      True if the minimization routine converged, False otherwise
    """
    # Run the minimization routine first time
    function = Minimize(initial, fun_moog, **options)
    try:
        parameters, converged = function.minimize()
        return parameters, converged
    except ValueError:
        print('No FeII lines were measured.')
        print('Skipping to next linelist..\n')
        return None


def _teffrangeRunner(linelist, parameters, converged, options):
    """Adjust the line list if the temperature is too low for the normal
    line list be Sousa+ 2008, to represent that of Tsantaki+ 2013.

    Input
    -----
    linelist : str
      The line list
    parameters : list
      The initial parameters for the minimization routine
    converged : bool
      Whether the minimization converged before. Useful is no lines are removed
    options : dict
      The options dictionary

    Output
    ------
    linelist : str
      The line list will have a new name if the outlier was removed
    parameters : list
      The initial parameters for the minimization routine
    converged : bool
      True if the minimization routine converged, False otherwise
    """
    d = np.loadtxt('rawLinelist/coolNormalDiff.lines')
    ll = np.loadtxt('linelist/%s' % linelist, skiprows=1, usecols=(0,))
    normalLL = np.in1d(ll, d)
    if np.any(normalLL) and (parameters[0] < 5200):
        print('Removing lines to compensate for low Teff\n')
        for li in ll[normalLL]:
            removeOutlier('linelist/%s' % linelist, li)

        # Restart the minimization procedure from the last best point
        parameters, converged = _minizationRunner(parameters, options)
        if options['outlier']:
            newName, parameters = _outlierRunner(options['outlier'], linelist, parameters, options)
            linelist = newName
        return linelist, parameters, converged
    else:
        return linelist, parameters, converged


def _autofixvtRunner(parameters, options):
    """
    Check and fix the microturbulence if it is close to the boundaries of the
    allowed range, i.e. 0.05 < vt < 9.95, and with a big slope of abundance vs.
    RW.

    Input
    -----
    parameters : list
      The initial parameters for the minimization routine
    options : dict
      The options dictionary

    Output
    ------
    parameters : list
      The initial parameters for the minimization routine now with fixed vt
    converged : bool
      True if the minimization routine converged, False otherwise
    """
    _, _, RWs, _, _ = fun_moog(parameters, options['model'], weight=options['weights'], version=options['MOOGv'])
    vt = parameters[-1]
    if ((vt < 0.05) and (abs(RWs) > 0.050)) or (vt > 5.0):
        options['fix_vt'] = True
        print('Running minimization with vt fixed...\n')
        parameters, converged = _minizationRunner(parameters, options)
        return parameters, converged
    return parameters, True


def _refineRunner(parameters, options):
    """
    Refine the parameters using stricter convergence criteria(66 percent stricter)

    Input
    -----
    parameters : list
      The initial parameters for the minimization routine
    options : dict
      The options dictionary

    Output
    ------
    parameters : list
      Final parameters obtained with the new convergence criteria
    options : dict
      The updated options dictionary with new convergence criteria
    """

    print('\nRefining the parameters with stricter convergence criteria...\n')
    options['EPcrit'] = round(options['EPcrit']/3, 4)
    options['RWcrit'] = round(options['RWcrit']/3, 4)
    options['ABdiffcrit'] = round(options['ABdiffcrit']/3, 4)
    p1, c = _minizationRunner(parameters, options)
    if c:
        print('Adjusting the final parameters...')
        parameters = p1  # overwrite with new best results
    return parameters, options


def _loggCorrections(parameters):
    """
    Function that corrects the logg values using the light curve and
    asteroseismic data. For additional documentation check Mortier et alu. 2014
    This correction is valid for stars in the Teff interval:
    5000 < Teff < 6500 K

    Input
    -----
    parameters : list
      The initial parameters for the minimization routine
    Output
    ------
    parameters : list
      The parameters with the corrected logg
    """
    parameters = list(parameters)

    # Lightcurve corrected logg
    loggLC = round(parameters[2] - 4.57E-4*parameters[0] + 2.59, 2)
    error_loggLC = round(np.sqrt((4.57e-4*parameters[1])**2 + (parameters[3])**2), 2)
    # Asteroseismic corrected logg
    loggastero = round(parameters[2] - 3.89E-4*parameters[0] + 2.10, 2)
    error_loggastero = round(np.sqrt((3.89e-4*parameters[1])**2 + (parameters[3])**2), 2)

    parameters.append(loggastero)
    parameters.append(error_loggastero)
    parameters.append(loggLC)
    parameters.append(error_loggLC)
    return parameters


def printToScreen(parameters, converged):
    """
    Function which prints the current parameters of the minimization routine

    Input
    -----
    parameters : list
      The current parameters of the minimization routine
    converged : bool
      True if the minimization routine converged, False otherwise

    Output
    ------
    Print to screen Teff, logg, [Fe/H], vt and whether it converged
    """
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


def sanityCheck(logger):
    """Check if all folders exists.

    Input
    -----
    logger : obj
      The logger object for the log file

    Output
    ------
    logger : obj
      The updated logger object if some tests did not pass
    """
    if not os.path.isdir('linelist'):
        logger.error('Error: The directory linelist does not exist!')
        os.mkdir('linelist')
        logger.info('linelist directory was created\n')
        raise IOError('linelist directory did not exist! Put the linelists inside that directory, please.')

    # Create results directory
    if not os.path.isdir('results'):
        os.mkdir('results')
        logger.info('results directory was created')


def _logging():
    """Create the logger object for the log file

    Output
    ------
    logger : obj
      The logger object for the log file
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
    return logger


def ewdriver(starLines='StarMe_ew.cfg', overwrite=None):
    """The function that glues everything together for the EW method

    Input
    -----
    starLines : str
      Configuration file (default: StarMe_ew.cfg)
    overwrite : bool
      Overwrite the results.csv file (default: False)

    Output
    ------
    <linelist>.(NC).out : file
      The output line list; NC=not converged.
    results.csv : file
      Easy readable table with results from many linelists
    """

    logger = _logging()
    try:
        sanityCheck(logger)
    except IOError, e:
        raise IOError(e)

    # Creating the output file
    _output(overwrite=overwrite, header=True)

    for (initial, options, line) in genStar(starLines):
        logger.info('Start with line list: %s' % line[0])
        logger.info('Initial parameters: {:d}, {:.2f}, {:.2f}, {:.2f}'.format(*initial))
        options = _prepare(line[0], initial, options)
        if options is None:
            continue  # The line list does not exists

        logger.info('Starting the initial minimization routine...')
        status = _minizationRunner(initial, options)
        if status is None:
            logger.error('The minimization routine did not finish succesfully.')
            continue  # Problem with the minimization routine
        else:
            parameters, converged = status
            logger.info('The minimization routine finished succesfully.')

        if options['outlier']:
            logger.info('Removing outliers.')
            newName, parameters = _outlierRunner(options['outlier'], line[0], parameters, options)
            line[0] = newName

        if options['teffrange']:
            logger.info('Correcting the line list, if necessary, for low Teff.')
            linelist, parameters, converged = _teffrangeRunner(line[0], parameters, converged, options)
            line[0] = linelist

        if options['autofixvt']:
            logger.info('Fixing vt if necessary.')
            parameters, converged = _autofixvtRunner(parameters, options)

        if options['refine'] and converged:
            logger.info('Refining the parameters.')
            parameters, options = _refineRunner(parameters, options)

        logger.info('Final parameters: {:.0f}, {:.2f}, {:.2f}, {:.2f}\n'.format(*parameters))
        _renaming(line[0], converged)
        parameters = error(line[0], converged, parameters,
                           atmtype=options['model'], version=options['MOOGv'],
                           weights=options['weights'])

        parameters = _loggCorrections(parameters)

        _output(linelist=line[0], parameters=parameters, converged=converged, options=options)

        printToScreen(parameters, converged)
    return parameters


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        cfgfile = sys.argv[1]
    else:
        cfgfile = 'StarMe_ew.cfg'
    parameters = ewdriver(starLines=cfgfile)
