#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import os
import logging
import numpy as np
import pandas as pd
from loggf_update import update_loggf
from interpolation import interpolator
from utils import _update_par, _run_moog, Readmoog

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


class AbundanceDriver:

    def __init__(self, cfgfile='StarMe_abund.cfg', overwrite=None):
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
            logger.info('results directory was created')

    def _options(self, options=None):
        '''Reads the options inside the config file'''
        defaults = {'model': 'kurucz95',
                    'MOOGv': 2014
                    }
        if options is None:
            self.options = defaults
        else:
            options = options.split(',')
            for option in options:
                if ':' in option:
                    option = option.split(':')
                    defaults[option[0]] = option[1]
                else:
                    defaults[option] = True
            defaults['model'] = defaults['model'].lower()
            self.options = defaults

    def save(self):
        """Write results"""
        linelist = self.abundance_dict.pop('linelist')
        teff = self.abundance_dict.pop('Temperature')
        logg = self.abundance_dict.pop('Gravity')
        feh = self.abundance_dict.pop('[Fe/H]')
        vt = self.abundance_dict.pop('microturbulence')
        elements = self.abundance_dict.keys()

        def _get_header(elements):
            """Get the header and append new elements to it"""
            with open('abundresults.dat', 'r') as f:
                header = f.readline()
            header = header.strip('\n').split(',')
            for element in elements:
                if element not in header:
                    header.append(element)
            return ','.join(header)

        if not os.path.isfile('abundresults.dat'):
            # The file does not exists, so create it and return from here
            header = 'linelist,temperature,logg,[Fe/H],vt,' + ','.join(elements)
            with open('abundresults.dat', 'w') as fout:
                fout.writelines(header + '\n')
                line = '%s,%i,%.2f,%.2f,%.2f' % (linelist, teff, logg, feh, vt)
                for element in elements:
                    line += ',%s' % self.abundance_dict[element]
                fout.writelines(line + '\n')
            return
        else:
            header = _get_header(elements)

        if self.overwrite:
            with open('abundresults.dat', 'w') as fout:
                fout.writelines(header + '\n')
        else:
            try:
                df = pd.read_csv('abundresults.dat', na_values='...')
                header = header.split(',')
                for key in header:
                    if key not in df.columns:
                        df[key] = np.nan
                df.to_csv(path_or_buf='abundresults.dat', header=header, index=False, na_rep='...')
            except IOError:
                # It does not exists yet, so create the file from scratch
                with open('abundresults.dat', 'w') as fout:
                    fout.writelines(header + '\n')

        df = pd.read_csv('abundresults.dat', na_values='...')
        rows = df.shape[0]
        for element in header[5:]:
            if element in elements:
                df[element] = self.abundance_dict[element]
            else:
                df[element] = np.nan
        df['linelist'] = linelist
        df['temperature'] = teff
        df['logg'] = logg
        df['[Fe/H]'] = feh
        df['vt'] = vt

        if rows:
            df.drop(df.index[range(rows-1)], inplace=True)

        df.to_csv(path_or_buf='abundresults.dat', header=False, index=False, mode='a', na_rep='...')

    def weighted_avg_and_std(self, values):
        """Get the weighted average and standard deviation.

        Input
        -----
        values : ndarray
          The list of values from which to calculate the weighted
          average and standard deviation

        Output
        ------
        average : float
          The weighted average
        std : float
          The weighted standard deviation
        """
        values = np.array(values)
        weights =  (np.abs(values-np.median(values))/(np.std(values)+1E-13)+0.25)
        weights_rounded = (np.round(weights/0.5+1E-13)*0.5)**(-1)
        average = round(np.average(values, weights=weights_rounded), 3)
        std = np.sqrt(np.average((values-average)**2, weights=weights_rounded))
        return average, std

    def abundancedriver(self):
        """The function that glues everything together

        Output
        ------
        abundresults.dat : file
          Easy readable table with results from many linelists
        """

        with open(self.cfgfile, 'r') as lines:
            for line in lines:
                self.abundance_dict = {}
                if not line[0].isalpha():
                    self.logger.debug('Skipping header: %s' % line.strip())
                    continue
                self.logger.info('Line list: %s' % line.strip())
                line = line.strip()
                line = line.split(' ')

                # Check if the linelist is inside the directory if not log it and pass to next linelist
                if not os.path.isfile('linelist/%s' % line[0]):
                    self.logger.error('Error: linelist/%s not found.' % line[0])
                    continue
                else:
                    _update_par(line_list='linelist/%s' % line[0])

                if len(line) == 1:
                    self.initial = (5777, 4.44, 0.00, 1.00)
                    options = _options()
                    self.logger.info('Setting solar values {0}, {1}, {2}, {3}'.format(*self.initial))

                elif len(line) == 5:
                    self.logger.info('Initial parameters given by the user.')
                    self.initial = map(float, line[1::])
                    self._options()
                    self.logger.info('Initial parameters: {0}, {1}, {2}, {3}'.format(*self.initial))

                elif len(line) == 6:
                    self.logger.info('Initial parameters given by the user.')
                    self.initial = map(float, line[1:-1])
                    self._options(line[-1])
                    self.logger.info('Initial parameters: {0}, {1}, {2}, {3}'.format(*self.initial))
                else:
                    self.logger.error('Could not process information for this line list: %s' % line)
                    continue

                # Setting the models to use
                if self.options['model'] not in ['kurucz95', 'marcs', 'apogee_kurucz']:
                    self.logger.error('Your request for type: %s is not available' % self.options['model'])
                    continue

                update_loggf(self.options['model'], 'linelist/%s' % line[0], region='ABoptical')
                # Get the initial grid models
                self.logger.info('Interpolation of model...')
                interpolator(params=self.initial, atmtype=self.options['model'])
                self.logger.info('Interpolation successful.')
                _run_moog()

                table = Readmoog(version=self.options['MOOGv']).all_table()
                elements = table.atom.unique()
                self.abundance_dict = {'linelist': line[0],
                                       'Temperature': self.initial[0],
                                       'Gravity': self.initial[1],
                                       '[Fe/H]': self.initial[2],
                                       'microturbulence': self.initial[3]}

                for element in elements:
                    sub_table = table[table.atom == element]
                    if len(sub_table) > 1:
                        abundance, _ = self.weighted_avg_and_std(sub_table.abund.values)
                    else:
                        abundance = sub_table.abund.values[0]
                    self.abundance_dict[element] = abundance

                self.save()

        return self.abundance_dict

    def to_screen(self):
        df = pd.read_csv('abundresults.dat')

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


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        cfgfile = sys.argv[1]
    else:
        cfgfile = 'StarMe_abund.cfg'

    driver = AbundanceDriver(cfgfile=cfgfile)
    _ = driver.abundancedriver()
    driver.to_screen()
