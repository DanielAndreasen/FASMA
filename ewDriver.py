#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import os
import yaml
import logging
import numpy as np
from shutil import copyfile
from minimization import Minimize
from loggf_update import update_loggf
from interpolation import interpolator
from utils import fun_moog, Readmoog, _update_par, error


class EWmethod:

    def __init__(self, cfgfile='StarMe_ew.cfg', overwrite=None):
        """The function that glues everything together for the EW method

        Input
        -----
        cfgfile : str
          Configuration file (default: StarMe_ew.cfg)
        overwrite : bool
          Overwrite the EWresults.dat file (default: False)

        Output
        ------
        <linelist>.(NC).out : file
          The output line list; NC=not converged.
        EWresults.dat : file
          Easy readable table with results from many linelists
        """
        self.cfgfile = cfgfile
        self.overwrite = overwrite

        # Setup of logger
        if os.path.isfile('captain.log'):  # Cleaning from previous runs
            os.remove('captain.log')
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)
        handler = logging.FileHandler('captain.log')
        handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)
        self._sanityCheck()

    def _sanityCheck(self):
        """Check if all folders exists."""
        if not os.path.isdir('linelist'):
            self.logger.error('Error: The directory linelist does not exist!')
            os.mkdir('linelist')
            self.logger.info('linelist directory was created\n')
            raise IOError('linelist directory did not exist! Put the linelists inside that directory, please.')

        # Create results directory
        if not os.path.isdir('results'):
            os.mkdir('results')
            self.logger.info('results directory was created')

    def _setup(self, line):
        """Do the setup with initial parameters and options.

        Input
        -----
        line : list
          A line from the configuration file after being split at spaces
        """
        self.linelist = line[0]
        if len(line) == 1:
            self.initial = [5777, 4.44, 0.00, 1.00]
            self._options()
        elif len(line) == 5:
            self.initial = map(float, line[1::])
            self.initial[0] = int(self.initial[0])
            self._options()
        elif len(line) == 2:
            self._options(line[1])
            self.initial = [5777, 4.44, 0.00, 1.00]
            if self.options['spt']:
                self.teff, self.logg = self._getSpt(self.options['spt'])
                self.vt = self._getMic(self.teff, self.logg, 0.00)
                self.initial = (self.teff, self.logg, 0.00, self.vt)
            if self.options['tmcalc']:
                self._tmcalc()
                self.initial = [self.teff, self.logg, self.feh, self.vt]
        elif len(line) == 6:
            self.initial = map(float, line[1:-1])
            self.initial[0] = int(self.initial[0])
            self._options(line[-1])
            if self.options['tmcalc']:
                self._tmcalc()
                self.initial = [self.teff, self.logg, self.feh, self.vt]

    def _getSpt(self):
        """Get the spectral type from a string like 'F5V'."""
        spt = self.options['spt']
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
            self.teff = 5777
            self.logg = 4.44
        try:
            line = line.split()
            self.teff = int(line[0])
            self.logg = float(line[1])
        except AttributeError:
            self.teff = line
            self.logg = 4.44

    def _getMic(self):
        """Calculate micro turbulence based on emperical relations.
        For dwarfs (logg>=3.95) we use the relation by Tsantaki+ 2013,
        and for giants (logg<3.95) we use the relation by Adibekyan+ 2015."""
        if self.logg >= 3.95:  # Dwarfs Tsantaki 2013
            self.vt = 6.932 * self.teff * (10**(-4)) - 0.348 * self.logg - 1.437
            self.vt = round(self.vt, 2)
        else:  # Giants Adibekyan 2015
            self.vt = 2.72 - (0.457 * self.logg) + (0.072 * self.feh)
            self.vt = round(self.vt, 2)

    def _tmcalc(self):
        """Initial guess on atmospheric parameters. Estimate based on TMCalc."""
        import sys
        sys.path.append('TMCALC/tmcalc_cython')
        from tmcalc_module import get_temperature_py as get_teff
        from tmcalc_module import get_feh_py as get_feh

        data = np.loadtxt('linelist/%s' % self.linelist, skiprows=1, usecols=(0, 4))
        X = np.zeros((data.shape[0], 9))
        X[:, 0] = data[:, 0]
        X[:, 4] = data[:, 1]
        np.savetxt('tmp.ares', X, '%.2f')

        teff = get_teff('TMCALC/tmcalc_cython/gteixeira_teff_cal.dat', 'tmp.ares')
        self.feh = get_feh('TMCALC/tmcalc_cython/gteixeira_feh_cal.dat', 'tmp.ares', teff[0], teff[1], teff[2], teff[3])[0]
        self.teff = teff[0]
        self.logg = 4.44
        self._getMic()
        os.remove('tmp.ares')

    def _renaming(self):
        """Save the output in a file related to the linelist."""
        if self.converged:
            copyfile('summary.out', 'results/%s.out' % self.linelist)
        else:
            copyfile('summary.out', 'results/%s.NC.out' % self.linelist)

    def _options(self, options=None):
        """Reads the options inside the config file."""
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
            self.options = defaults
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
            self.options = defaults

    def _genStar(self):
        """A generator for the configuration file."""
        lines = open(self.cfgfile, 'r')
        for line in lines:
            if not line[0].isalnum():
                # Header
                continue
            line = line.strip()
            line = line.split(' ')
            if len(line) not in [1, 2, 5, 6]:
                # Not the expected format
                continue
            self._setup(line)
            yield self.initial, self.options, line

    def _prepare(self):
        """Prepare the run with setup and first interpolation."""
        if not os.path.isfile('linelist/%s' % self.linelist):
            return None
        else:
            _update_par(line_list='linelist/%s' % self.linelist)

        # Make the initial interpolation
        interpolator(params=self.initial, atmtype=self.options['model'])

        # Adjusting the options for the minimization routine
        if __name__ in ('__main__', 'ewDriver'):
            self.options['GUI'] = False  # Running batch mode
        else:
            self.options['GUI'] = True  # Running GUI mode

        # Update loggf values
        with open('linelist/%s' % self.linelist, 'r') as lines:
            line = lines.readlines()[-1]
        line = filter(None, line.split())
        w = float(line[0])
        if w < 7000:
            region = 'EWoptical'
        else:
            region = 'EWNIR'
        update_loggf(self.options['model'], 'linelist/%s' % self.linelist, region=region)

    def _output(self, header=None):
        """Create the output file 'EWresults.dat'."""

        hdr = ['linelist', 'teff', 'tefferr', 'logg', 'loggerr', 'feh', 'feherr',
               'vt', 'vterr', 'loggastero', 'dloggastero', 'loggLC', 'dloggLC',
               'convergence', 'fixteff', 'fixlogg', 'fixfeh', 'fixvt', 'outlier',
               'weights', 'model', 'refine', 'EPcrit', 'RWcrit', 'ABdiffcrit']
        if header is not None:
            if self.overwrite:
                with open('EWresults.dat', 'w') as output:
                    output.write('\t'.join(hdr)+'\n')
            else:
                if not os.path.isfile('EWresults.dat'):
                    with open('EWresults.dat', 'w') as output:
                        output.write('\t'.join(hdr)+'\n')
        else:
            tmp = [self.linelist] + self.parameters +\
                  [self.converged, self.options['fix_teff'],
                   self.options['fix_logg'], self.options['fix_feh'],
                   self.options['fix_vt'], self.options['outlier']] +\
                  [self.options['weights'], self.options['model'],
                   self.options['refine'], self.options['EPcrit'],
                   self.options['RWcrit'], self.options['ABdiffcrit']]
            with open('EWresults.dat', 'a') as output:
                output.write('\t'.join(map(str, tmp))+'\n')

    def _printToScreen(self):
        """
        Function which prints the current parameters of the minimization routine."""

        if __name__ == '__main__':
            if self.converged:
                print('\nCongratulation, you have won! Your final parameters are:')
            else:
                print('\nSorry, you did not win. However, your final parameters are:')
            try:
                print(u' Teff:{:>8d}\u00B1{:.0f}\n logg:{:>8.2f}\u00B1{:1.2f}\n [Fe/H]:{:>+6.2f}\u00B1{:1.2f}\n vt:{:>10.2f}\u00B1{:1.2f}\n\n\n\n'.format(*self.parameters))
            except UnicodeEncodeError:
                print(' Teff:{:>8d}({:.0f})\n logg:{:>8.2f}({:1.2f})\n [Fe/H]:{:>+6.2f}({:1.2f})\n vt:{:>10.2f}({:1.2f})\n\n\n\n'.format(*self.parameters))
        elif __name__ == 'ewDriver':
            if self.converged:
                print('\nCongratulation, you have won! Your final parameters are:')
            else:
                print('\nSorry, you did not win. However, your final parameters are:')
            print(u' Teff:{:>8d}+/-{:.0f}\n logg:{:>8.2f}+/-{:1.2f}\n [Fe/H]:{:>+6.2f}+/-{:1.2f}\n vt:{:>10.2f}+/-{:1.2f}\n\n\n\n'.format(*self.parameters))

    def minizationRunner(self, p=None):
        """A function to run the minimization routine

        Output
        ------
        _ : bool
          True if the minimization run succesfully
        """
        # Run the minimization routine first time
        if p is not None:
            function = Minimize(p, fun_moog, **self.options)
        else:
            function = Minimize(self.initial, fun_moog, **self.options)
        try:
            self.parameters, self.converged = function.minimize()
            return True
        except ValueError:
            print('No FeII lines were measured.')
            print('Skipping to next linelist..\n')
            return None

    def outlierRunner(self):
        """Remove the potential outliers based on a given method. After outliers
        are removed, then restarts the minimization routine at the previous best
        found parameters."""
        type = self.options['outlier']
        tmpll = 'linelist/tmplinelist.moog'
        copyfile('linelist/'+self.linelist, tmpll)
        _update_par(line_list=tmpll)
        newLineList = False
        Noutlier = 0
        outliers = self._hasOutlier()
        if type == '1Iter':
            # Remove one outlier above 3 sigma iteratively
            while outliers:
                Noutlier += 1
                newLineList = True  # At the end, create a new linelist
                wavelength = outliers[max(outliers.keys())]
                self.removeOutlier(tmpll, wavelength)
                print('Removing line: %.2f. Outliers removed: %d' % (wavelength, Noutlier))
                print('Restarting the minimization routine...\n')
                function = Minimize(self.parameters, fun_moog, **self.options)
                self.parameters, self.converged = function.minimize()
                outliers = self._hasOutlier()

        elif type == '1Once':
            # Remove one outlier above 3 sigma once
            if outliers:
                Noutlier += 1
                newLineList = True  # At the end, create a new linelist
                wavelength = outliers[max(outliers.keys())]
                self.removeOutlier(tmpll, wavelength)
                print('Removing line: %.2f. Outliers removed: %d' % (wavelength, Noutlier))
                print('Restarting the minimization routine...\n')
                function = Minimize(self.parameters, fun_moog, **self.options)
                self.parameters, self.converged = function.minimize()
                outliers = self._hasOutlier()

        elif type == 'allIter':
            # Remove all outliers above 3 sigma iteratively
            while outliers:
                newLineList = True  # At the end, create a new linelist
                for wavelength in outliers.itervalues():
                    self.removeOutlier(tmpll, wavelength)
                    Noutlier += 1
                    print('Removing line: %.2f. Outliers removed: %d' % (wavelength, Noutlier))
                print('Restarting the minimization routine...\n')
                function = Minimize(self.parameters, fun_moog, **self.options)
                self.parameters, self.converged = function.minimize()
                outliers = self._hasOutlier()

        elif type == 'allOnce':
            # Remove all outliers above 3 sigma once
            if outliers:
                newLineList = True  # At the end, create a new linelist
                for wavelength in outliers.itervalues():
                    self.removeOutlier(tmpll, wavelength)
                    Noutlier += 1
                    print('Removing line: %.2f. Outliers removed: %d' % (wavelength, Noutlier))
                print('Restarting the minimization routine...\n')
                function = Minimize(self.parameters, fun_moog, **self.options)
                self.parameters, self.converged = function.minimize()
                outliers = self._hasOutlier()

        if newLineList:
            self.linelist = self.linelist.replace('.moog', '_outlier.moog')
            copyfile(tmpll, 'linelist/'+self.linelist)
        _update_par(line_list='linelist/'+self.linelist)
        os.remove(tmpll)

    def _hasOutlier(self):
        """Function that reads the summary.out file and return a dictionary
        with key being the deviation (above n sigma), and value the wavelength.

        Output
        ------
        d : dict
          A dictionary with {n*deviation: wavelength}
        """
        MOOGv = self.options['MOOGv']
        n = self.options['sigma']
        idx = 1 if MOOGv > 2013 else 0
        n = float(n)
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

    def removeOutlier(self, fname, wavelength):
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
        wavelength = str(wavelength)[:-1]
        with open(fname, 'r') as lines:
            fout = ''
            for line in lines:
                if line.replace(' ', '').startswith(wavelength):
                    continue
                fout += line
        with open(fname, 'w') as f:
            f.writelines(fout)

    def teffrangeRunner(self):
        """Adjust the line list if the temperature is too low for the normal
        line list be Sousa+ 2008, to represent that of Tsantaki+ 2013."""

        d = np.loadtxt('rawLinelist/coolNormalDiff.lines')
        wavelengths = np.loadtxt('linelist/%s' % self.linelist, skiprows=1, usecols=(0,))
        normalLL = np.in1d(wavelengths, d)
        if np.any(normalLL) and (self.parameters[0] < 5200):
            print('Removing lines to compensate for low Teff\n')
            for wavelength in wavelengths[normalLL]:
                self.removeOutlier('linelist/%s' % self.linelist, wavelength)

            # Restart the minimization procedure from the last best point
            _ = self.minizationRunner()
            if self.options['outlier']:
                self.outlierRunner()

    def autofixvtRunner(self):
        """Check and fix the microturbulence if it is close to the boundaries of
        the allowed range, i.e. 0.05 < vt < 9.95, and with a big slope of
        abundance vs. RW."""
        _, _, RWs, _, _ = fun_moog(self.parameters,
                                   self.options['model'],
                                   weight=self.options['weights'],
                                   version=self.options['MOOGv'])
        vt = self.parameters[-1]
        if ((vt < 0.05) and (abs(RWs) > 0.050)) or (vt > 5.0):
            self.options['fix_vt'] = True
            print('Running minimization with vt fixed...\n')
            _ = self.minizationRunner()

    def refineRunner(self):
        """Refine the parameters using stricter convergence criteria:
        66% stricter."""

        print('\nRefining the parameters with stricter convergence criteria...\n')
        self.options['EPcrit'] = round(self.options['EPcrit']/3, 4)
        self.options['RWcrit'] = round(self.options['RWcrit']/3, 4)
        self.options['ABdiffcrit'] = round(self.options['ABdiffcrit']/3, 4)
        p = tuple(self.parameters)
        _ = self.minizationRunner(p=self.parameters)
        if self.converged:
            print('Adjusted the final parameters...')
        else:
            # Use old parameters and set converged to True, since this is a
            # criteria to even use this functionality.
            self.parameters = list(p)
            self.converged = True

    def loggCorrections(self):
        """
        Function that corrects the logg values using the light curve and
        asteroseismic data. For additional documentation check Mortier+ 2014.
        This correction is valid for stars in the Teff interval:
        5000 < Teff < 6500 K.
        """
        self.parameters = list(self.parameters)

        # Lightcurve corrected logg
        loggLC = round(self.parameters[2] - 4.57E-4*self.parameters[0] + 2.59, 2)
        error_loggLC = round(np.sqrt((4.57e-4*self.parameters[1])**2 + (self.parameters[3])**2), 2)
        # Asteroseismic corrected logg
        loggastero = round(self.parameters[2] - 3.89E-4*self.parameters[0] + 2.10, 2)
        error_loggastero = round(np.sqrt((3.89e-4*self.parameters[1])**2 + (self.parameters[3])**2), 2)

        self.parameters.append(loggastero)
        self.parameters.append(error_loggastero)
        self.parameters.append(loggLC)
        self.parameters.append(error_loggLC)

    def ewdriver(self):
        # Creating the output file
        self._output(header=True)

        for (self.initial, self.options, self.line) in self._genStar():
            self.logger.info('Start with line list: %s' % self.linelist)
            self.logger.info('Initial parameters: {:.0f}, {:.2f}, {:.2f}, {:.2f}'.format(*self.initial))
            self._prepare()
            if self.options is None:
                self.logger.error('The line list does not exists!\n')
                continue  # The line list does not exists

            self.logger.info('Starting the initial minimization routine...')
            status = self.minizationRunner()
            if status is None:
                self.logger.error('The minimization routine did not finish succesfully.')
                continue  # Problem with the minimization routine
            else:
                self.logger.info('The minimization routine finished succesfully.')

            if self.options['outlier']:
                self.logger.info('Removing outliers.')
                self.outlierRunner()

            if self.options['teffrange']:
                self.logger.info('Correcting the line list, if necessary, for low Teff.')
                self.teffrangeRunner()

            if self.options['autofixvt']:
                self.logger.info('Fixing vt if necessary.')
                self.autofixvtRunner()

            if self.options['refine'] and self.converged:
                self.logger.info('Refining the parameters.')
                self.refineRunner()

            self.logger.info('Final parameters: {:.0f}, {:.2f}, {:.2f}, {:.2f}\n'.format(*self.parameters))
            self._renaming()
            self.parameters = error(self.linelist, self.converged,
                                    self.parameters,
                                    atmtype=self.options['model'],
                                    version=self.options['MOOGv'],
                                    weights=self.options['weights'])

            self.loggCorrections()
            self._output()
            self._printToScreen()
        return self.parameters


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        cfgfile = sys.argv[1]
    else:
        cfgfile = 'StarMe_ew.cfg'
    driver = EWmethod(cfgfile=cfgfile, overwrite=None)
    parameters = driver.ewdriver()
