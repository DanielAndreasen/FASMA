#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import os
import pandas as pd
import numpy as np
from aresDriver import aresdriver
from ewDriver import ewdriver
from abundanceDriver import abundancedriver
from time import time

'''Get from the spectrum to parameters and abundances'''


class FullSpectralAnalysis:

    def __init__(self, cfgfile='StarMe_all.cfg'):
        self.cfgfile = cfgfile

    def _options(self, options=None):
        '''Reads the options inside the config file.

        Input
        -----
        options : str (optional)
          The line from the configuration file with the user options

        Output
        ------
        defaults : dict
          A dictionary with all the options for the EW method.
        '''

        # There a general (gen) options, specific for parameters (par), and
        # specific for EW measurements (ews)
        defaults = {'gen.model': 'kurucz95',
                    'gen.MOOGv': 2014,
                    'par.fix_teff': False,
                    'par.fix_logg': False,
                    'par.fix_feh': False,
                    'par.fix_vt': False,
                    'par.refine': False,
                    'par.iterations': 160,
                    'par.EPcrit': 0.001,
                    'par.RWcrit': 0.003,
                    'par.ABdiffcrit': 0.01,
                    'par.outlier': False,
                    'par.teffrange': False,
                    'par.autofixvt': False,
                    'par.tmcalc': False,
                    'par.sigma': 3,
                    'ews.lambdai': '3900.0',
                    'ews.lambdaf': '25000.0',
                    'ews.smoothder': '4',
                    'ews.space': '2.0',
                    'ews.rejt': False,
                    'ews.lineresol': '0.07',
                    'ews.miniline': '2.0',
                    'ews.plots_flag': False,
                    'ews.EWcut': '200.0',
                    'ews.snr': False,
                    'ews.output': False,
                    'ews.rvmask': '"0,0"',
                    'ews.force': False,
                    'ews.extra': None}

        if options is None:
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
            defaults['gen.model'] = defaults['gen.model'].lower()
            defaults['gen.MOOGv'] = int(defaults['gen.MOOGv'])
            defaults['par.iterations'] = int(defaults['par.iterations'])
            defaults['par.EPcrit'] = float(defaults['par.EPcrit'])
            defaults['par.RWcrit'] = float(defaults['par.RWcrit'])
            defaults['par.ABdiffcrit'] = float(defaults['par.ABdiffcrit'])

            defaults['ews.lambdai'] = float(defaults['ews.lambdai'])
            defaults['ews.lambdaf'] = float(defaults['ews.lambdaf'])
            defaults['ews.smoothder'] = int(defaults['ews.smoothder'])
            defaults['ews.space'] = float(defaults['ews.space'])
            defaults['ews.lineresol'] = float(defaults['ews.lineresol'])
            defaults['ews.miniline'] = float(defaults['ews.miniline'])
            defaults['ews.rvmask'] = str(defaults['ews.rvmask'])
            defaults['ews.EWcut'] = float(defaults['ews.EWcut'])

            if defaults['par.outlier'] not in [False, '1Iter', '1Once', 'allIter', 'allOnce']:
                print('Invalid option set for option "outlier"')
                defaults['par.outlier'] = False

            if defaults['ews.plots_flag']:
                defaults['ews.plots_flag'] = '1'

            if not defaults['ews.rejt']:
                defaults['ews.rejt'] = '3;5764,5766,6047,6053,6068,6076'
            try:
                defaults['ews.rejt'] = float(defaults['ews.rejt'])
            except ValueError:
                defaults['ews.rejt'] = str(defaults['ews.rejt'])

        return defaults

    def ares(self):
        aresOptions = {}
        for option in self.options.iterkeys():
            if option.startswith('ews'):
                newoption = option.replace('ews.', '')
                aresOptions[newoption] = self.options[option]

        with open('StarMe_all1.cfg', 'w') as fout:
            opt1 = ''
            opt2 = []
            for ai in aresOptions.iterkeys():
                if aresOptions[ai] is True:
                    opt1 += '%s,' % ai
                elif (aresOptions[ai] is None) or (aresOptions[ai] is False):
                    continue
                else:
                    opt2.append('%s:%s' % (ai, aresOptions[ai]))

            opt2 = ','.join(opt2)
            opt = opt1 + opt2
            fout.write('Sousa2007_opt.lst {0} extra:Neves2009_elements.lst,{1}'.format(self.spectrum, opt))
        self.snr = aresdriver('StarMe_all1.cfg')

    def ewmethod(self):
        ewmethodOptions = {}
        for option in self.options.iterkeys():
            if option.startswith('par'):
                ewmethodOptions[option.replace('par.', '')] = self.options[option]
            elif option.startswith('gen'):
                ewmethodOptions[option.replace('gen.', '')] = self.options[option]

        with open('StarMe_all2.cfg', 'w') as fout:
            opt1 = ''
            opt2 = []
            for ai in ewmethodOptions.iterkeys():
                if ewmethodOptions[ai] is True:
                    opt1 += '%s,' % ai
                elif (ewmethodOptions[ai] is None) or (ewmethodOptions[ai] is False):
                    continue
                else:
                    opt2.append('%s:%s' % (ai, ewmethodOptions[ai]))

            opt2 = ','.join(opt2)
            opt = opt1 + opt2
            p = ' '.join(map(str, self.initial))
            linelist = self.spectrum.replace('.fits', '.moog')
            fout.write('{0} {1} {2}'.format(linelist, p, opt))
        self.params = ewdriver('StarMe_all2.cfg')[:-4]

    def abundances(self):
        abundanceOptions = {}
        for option in self.options.iterkeys():
            if option.startswith('gen'):
                abundanceOptions[option.replace('gen.', '')] = self.options[option]

        with open('StarMe_all3.cfg', 'w') as fout:
            opt1 = ''
            opt2 = []
            for ai in abundanceOptions.iterkeys():
                if abundanceOptions[ai] is True:
                    opt1 += '%s,' % ai
                elif (abundanceOptions[ai] is None) or (abundanceOptions[ai] is False):
                    continue
                else:
                    opt2.append('%s:%s' % (ai, abundanceOptions[ai]))

            opt2 = ','.join(opt2)
            opt = opt1 + opt2
            p = ' '.join(map(str, self.params[::2]))

            s = self.spectrum.rpartition('.')
            linelist = s[0] + '_sec.moog'
            fout.write('{0} {1} {2}'.format(linelist, p, opt))
        self.abundance = abundancedriver('StarMe_all3.cfg')

    def saveResults(self, dict):
        '''Would like to have:
        spectrum, SNR, parameters, abundances, options
        '''
        spectrum = dict.pop('spectrum')
        snr = dict.pop('SNR')
        teff = dict.pop('Teff')
        tefferr = dict.pop('Tefferr')
        logg = dict.pop('logg')
        loggerr = dict.pop('loggerr')
        feh = dict.pop('feh')
        feherr = dict.pop('feherr')
        vt = dict.pop('vt')
        vterr = dict.pop('vterr')
        options = dict.pop('options')
        dict = dict.pop('abundances')
        elements = dict.keys()

        def _get_header(elements):
            """Get the header and append new elements to it"""
            with open('FASMA_all.dat', 'r') as f:
                header = f.readline()
            header = header.strip('\n').split(',')
            for element in elements:
                if element not in header:
                    header.append(element)
            return ','.join(header)

        if not os.path.isfile('FASMA_all.dat'):
            # The file does not exists, so create it and return from here
            header = 'spectrum,SNR,Teff,dTeff,logg,dlogg,[Fe/H],d[Fe/H],vt,dvt,' + ','.join(elements)
            with open('FASMA_all.dat', 'w') as fout:
                fout.writelines(header + '\n')
                line = '%s,%i,%i,%i,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f' % (spectrum, snr, teff, tefferr, logg, loggerr, feh, feherr, vt, vterr)
                for element in elements:
                    line += ',%s' % dict[element]
                fout.writelines(line + '\n')
            return
        else:
            header = _get_header(elements)

        try:
            # Setting previous elements to nan if they do nos exists
            df = pd.read_csv('FASMA_all.dat', na_values='....')
            header = header.split(',')
            for key in header:
                if key not in df.columns:
                    df[key] = np.nan
            df.to_csv(path_or_buf='FASMA_all.dat', header=header, index=False, na_rep='....')
        except IOError:
            # It does not exists yet, so create the file from scratch
            with open('FASMA_all.dat', 'w') as fout:
                fout.writelines(header + '\n')

        df = pd.read_csv('FASMA_all.dat', na_values='....')
        rows = df.shape[0]
        for element in header[10:]:
            if element in elements:
                df[element] = dict[element]
            else:
                df[element] = np.nan
        df['spectrum'] = spectrum
        df['SNR'] = snr
        df['Teff'] = teff
        df['dTeff'] = tefferr
        df['logg'] = logg
        df['dlogg'] = loggerr
        df['[Fe/H]'] = feh
        df['d[Fe/H]'] = feherr
        df['vt'] = vt
        df['dvt'] = vterr

        if rows:
            df.drop(df.index[range(rows-1)], inplace=True)

        df.to_csv(path_or_buf='FASMA_all.dat', header=False, index=False, mode='a', na_rep='....')

    def get_all(self):
        '''Get EW measurements, parameters from EW method, and abundances'''
        with open(self.cfgfile, 'r') as stars:
            for star in stars:
                if not star[0].isalpha():  # Skip comments
                    continue
                star = star.strip().split(' ')
                self.spectrum = star[0]
                print '*' * 42
                s = ' Analyzing: %s ' % self.spectrum.rpartition('.')[0]
                print s.center(42, '*')
                print '*' * 42

                Nopt = len(star)
                if (Nopt == 1) or (Nopt == 5):  # Using pure defaults
                    self.options = self._options()
                else:  # Use user-defined options
                    self.options = self._options(star[-1])

                if Nopt >= 5:
                    self.initial = map(float, star[1:5])
                else:
                    self.initial = 5777, 4.44, 0.0, 1.0

                print '\nMeasuring EWs. Please wait...'
                t = time()
                self.ares()
                print 'Done in %.2fs!\n' % (time()-t)
                print 'Getting parameters. Please wait...'
                t = time()
                self.ewmethod()
                print 'Done in %.2fs!\n' % (time()-t)
                print 'Getting abundances. Please wait...'
                t = time()
                self.abundances()
                print 'Done in %.2fs!\n\n' % (time()-t)

                results = {'spectrum': self.spectrum,
                           'SNR': self.snr,
                           'Teff': self.params[0],
                           'Tefferr': self.params[1],
                           'logg': self.params[2],
                           'loggerr': self.params[3],
                           'feh': self.params[4],
                           'feherr': self.params[5],
                           'vt': self.params[6],
                           'vterr': self.params[7],
                           'abundances': self.abundance,
                           'options': self.options}

                self.saveResults(results)


if __name__ == '__main__':
    import sys
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)
    if len(sys.argv) > 1:
        cfgfile = sys.argv[1]
    else:
        cfgfile = 'StarMe_all.cfg'
    analysis = FullSpectralAnalysis(cfgfile)
    analysis.get_all()

    df = pd.read_csv('FASMA_all.dat')

    s = df.to_string(justify='right', formatters={'temperature': lambda x: '%d' % x,
                                                  'logg': lambda x: '%.2f' % x,
                                                  '[Fe/H]': lambda x: '%.2f' % x,
                                                  'vt': lambda x: '%.2f' % x})
    for i, line in enumerate(s.split('\n')):
        if i == 0:
            print line
            continue
        val = len(str(i-1))
        print ' '*val + line[val::]
