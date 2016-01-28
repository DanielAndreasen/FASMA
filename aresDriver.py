#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import logging
import os
import numpy as np


def _run_ares():
    """Run ARES"""
    os.system('ARES > /dev/null')
    for tmp in ['tmp', 'tmp2', 'tmp3']:
        if os.path.isfile(tmp):
            os.remove(tmp)


def make_linelist(line_file, ares, cut=200):
    """This function creates a MOOG readable file from the ARES output using the
    line list and atomic data of make_linelist.dat file"""

    # Read the line list and check for multiple identical lines
    linelist = np.loadtxt(line_file, skiprows=2, usecols=range(4))
    assert (len(np.unique(linelist[:, 0])) == len(linelist[:, 0])), 'Check for multiple lines the linelist: %s' % line_file

    # Read the lines and ews in the ares data and check for identical lines
    data = np.loadtxt(ares, usecols=(0, 4))
    _, idx = np.unique(data[:, 0], return_index=True)
    data = data[idx]
    idx = np.argsort(data[:, 0])
    data = data[idx]

    # Cut high EW lines away
    idx = data[:, 1] > cut
    print('%s line(s) with EW higher than %s were deleted' % (len(data[idx, 0]), cut))
    data = data[~idx]
    # Wavelength and EW taken from the ares file.
    # Test whether each element of a 1D array is also present in a second array
    idx = np.in1d(data[:, 0], linelist[:, 0])
    data = data[idx]

    # Sort common elements from ares by wavelength
    idx = np.argsort(data[:, 0])
    data = data[idx]

    # Wavelength and atomic data taken from the make_linelist.dat file.
    # Test whether each element of a 1D array is also present in a second array
    linelist_index = np.in1d(linelist[:, 0], data[:, 0])
    index_lines_not_found = np.in1d(linelist[:, 0], data[:, 0])
    lines_not_found = linelist[~index_lines_not_found, 0]
    print('ARES did not find %i lines' % len(lines_not_found))

    linelist = linelist[linelist_index]
    print('Lines in the new line list: %i' % len(linelist[:, 0]))

    # Sort common elements from line list by wavelength
    idx = np.argsort(linelist[:, 0])
    linelist = linelist[idx]

    # Merge line list data with the EW from ARES
    # Sort the FeI and the FeII lines using the atomic number
    values = np.column_stack((data[:, 0], linelist[:, 1], linelist[:, 2], linelist[:, 3], data[:, 1]))
    idx = np.argsort(values[:, 1])
    sorted_values = values[idx]

    # Write results in MOOG readable format
    assert np.array_equal(data[:, 0], linelist[:, 0]), 'There is something wrong with the common elements of ARES and the line list'
    data = zip(sorted_values[:, 0], sorted_values[:, 1], sorted_values[:, 2], sorted_values[:, 3], sorted_values[:, 4])
    fout = '%s.moog' % ares.rpartition('.')[0]
    np.savetxt(fout, data, fmt=('%9.3f', '%10.1f', '%9.2f', '%9.3f', '%28.1f'), header=' %s' % ares)
    os.remove(ares)


def _options(options=False):
    '''Reads the options inside the config file'''
    defaults = {'lambdai':'3900.0',
                'lambdaf': '25000.0',
                'smoothder': '4',
                'space': '2.0',
                'rejt': '0.995',
                'lineresol' : '0.07',
                'miniline' : '2.0',
                'plots_flag': False,
                'EWcut': '200.0',
                'snr': False,
                'rvmask' : '"0,0"'
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
        defaults['lambdai'] = float(defaults['lambdai'])
        defaults['lambdaf'] = float(defaults['lambdaf'])
        defaults['smoothder'] = int(defaults['smoothder'])
        defaults['space'] = float(defaults['space'])
        defaults['rejt'] = float(defaults['rejt'])
        defaults['lineresol'] = float(defaults['lineresol'])
        defaults['miniline'] = float(defaults['miniline'])
        defaults['rvmask'] = str(defaults['rvmask'])
        defaults['EWcut'] = float(defaults['EWcut'])
        if defaults['plots_flag']:
            defaults['plots_flag'] = '1'
        return defaults


def update_ares(line_list, spectrum, out, options):
    """Driver for ARES"""

    default_options = options
    for key, value in default_options.iteritems():
        if key not in options.keys():
            options[key] = value

    def rejt_from_snr(snr):
        """Calculate rejt from SNR"""
        return 1.0-(1.0/snr)

    if options['snr']:
        rejt = rejt_from_snr(options['snr'])
    else:
        rejt = options['rejt']
    rejt = 0.999 if rejt > 0.999 else rejt
    plot = 1 if options['plots_flag'] else 0

    fout = 'specfits=\'spectra/%s\'\n' % spectrum
    fout += 'readlinedat=\'rawLinelist/%s\'\n' % line_list
    fout += 'fileout=\'linelist/%s\'\n' % out
    fout += 'lambdai=%s\n' % options['lambdai']
    fout += 'lambdaf=%s\n' % options['lambdaf']
    fout += 'smoothder=%s\n' % options['smoothder']
    fout += 'space=%s\n' % options['space']
    fout += 'rejt=%s\n' % rejt
    fout += 'lineresol=%s\n' % options['lineresol']
    fout += 'miniline=%s\n' % options['miniline']
    fout += 'plots_flag=%s\n' % plot

    with open('mine.opt', 'w') as f:
        f.writelines(fout)


def aresdriver(starLines='StarMe_ares.cfg'):
    """The function that glues everything together

    Input:
    starLines   -   Configuration file (default: StarMe_ares.cfg)

    Output:
    <linelist>.out          -   Output file
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

    #Check if there is a directory called linelist, if not create it and ask the user to put files there
    if not os.path.isdir('linelist'):
        os.mkdir('linelist')
        logger.info('linelist directory was created')

    if not os.path.isdir('rawLinelist'):
        os.mkdir('rawLinelist')
        logger.info('linelist directory was created')
        raise IOError('Please put linelists in rawLinelist folder')

    with open(starLines, 'r') as lines:
        for line in lines:
            if not line[0].isalpha():
                logger.debug('Skipping header: %s' % line.strip())
                continue
            logger.info('Processing: %s' % line.strip())
            line = line.strip()
            line = line.split(' ')

            if len(line) == 3:
                options = _options(line[-1])
                line_list = line[0]
                spectrum = line[1]
                out = '%s.ares' % spectrum.rpartition('.')[0]
                update_ares(line_list, spectrum, out, options)
            else:
                logger.error('Could not process information for this line: %s' % line)
                continue
            _run_ares()
            line_list = line[0]
            make_linelist('rawLinelist/'+line_list, 'linelist/'+out, cut=options['EWcut'])

if __name__ == '__main__':
    aresdriver(starLines='StarMe_ares.cfg')